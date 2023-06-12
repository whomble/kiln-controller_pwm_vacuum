import threading
import time
import random
import datetime
import logging
import json
import config
import os

log = logging.getLogger(__name__)

class DupFilter(object):
    def __init__(self):
        self.msgs = set()

    def filter(self, record):
        rv = record.msg not in self.msgs
        self.msgs.add(record.msg)
        return rv

class Duplogger():
    def __init__(self):
        self.log = logging.getLogger("%s.dupfree" % (__name__))
        dup_filter = DupFilter()
        self.log.addFilter(dup_filter)
    def logref(self):
        return self.log

duplog = Duplogger().logref()


class Output(object):
    def __init__(self):
        self.active = False
        self.load_libs()

    def load_libs(self):
        try:
            import RPi.GPIO as GPIO
            GPIO.setmode(GPIO.BCM)
            GPIO.setwarnings(False)
            GPIO.setup(config.gpio_heat_pwm, GPIO.OUT)
            GPIO.setup(config.gpio_primaire, GPIO.OUT)    # Sortie Primaire
            GPIO.setup(config.gpio_secondaire, GPIO.OUT)  # Sortie Secondaire
            GPIO.setup(config.gpio_argon_in, GPIO.OUT)    # Sortie Argon_In
            GPIO.setup(config.gpio_argon_out, GPIO.OUT)   # Sortie Argon_Out
            pwm_heat = GPIO.PWM(config.gpio_heat_pwm, 1000)  # Configuration du pin pwm_heat en mode pwm
            pwm_heat.start(0)
            self.active = True
            self.GPIO = GPIO
        except:
            msg = "Could not initialize GPIOs, oven operation will only be simulated!"
            log.warning(msg)
            self.active = False

    def heat(self, pwm_value):
        self.GPIO.output(config.gpio_heat, self.GPIO.HIGH)
        config.gpio_heat_pwm.ChangeDutyCycle(pwm_value)

    def Primaire(self, value):
        self.GPIO.output(config.gpio_primaire, value)

    def Secondaire(self, value):
        self.GPIO.output(config.gpio_secondaire, value)

    def Argon_In(self, value):
        self.GPIO.output(config.gpio_argon_in, value)

    def Argon_Out(self, value):
        self.GPIO.output(config.gpio_argon_out, value)

class Board(object):
    def __init__(self):
        self.name = None
        self.active = False
        self.temp_sensor = None
        self.gpio_active = False
        self.load_libs()
        self.create_temp_sensor()
        self.temp_sensor.start()

    def load_libs(self):
        if config.max31855:
            try:
                #from max31855 import MAX31855, MAX31855Error
                self.name='MAX31855'
                self.active = True
                log.info("import %s " % (self.name))
            except ImportError:
                msg = "max31855 config set, but import failed"
                log.warning(msg)

        if config.max31856:
            try:
                #from max31856 import MAX31856, MAX31856Error
                self.name='MAX31856'
                self.active = True
                log.info("import %s " % (self.name))
            except ImportError:
                msg = "max31856 config set, but import failed"
                log.warning(msg)

    def create_temp_sensor(self):
        if config.simulate == True:
            self.temp_sensor = TempSensorSimulated()
        else:
            self.temp_sensor = TempSensorReal()

class BoardSimulated(object):
    def __init__(self):
        self.temp_sensor = TempSensorSimulated()


class TempSensor(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.daemon = True
        self.temperature = 0
        self.bad_percent = 0
        self.time_step = config.sensor_time_wait
        self.noConnection = self.shortToGround = self.shortToVCC = self.unknownError = False

class TempSensorSimulated(TempSensor):
    '''not much here, just need to be able to set the temperature'''
    def __init__(self):
        TempSensor.__init__(self)

class TempSensorReal(TempSensor):
    '''real temperature sensor thread that takes N measurements
       during the time_step'''
    def __init__(self):
        TempSensor.__init__(self)
        self.sleeptime = self.time_step / float(config.temperature_average_samples)
        self.bad_count = 0
        self.ok_count = 0
        self.bad_stamp = 0

        if config.max31855:
            log.info("init MAX31855")
            from max31855 import MAX31855, MAX31855Error
            self.thermocouple = MAX31855(config.gpio_sensor_cs,
                                     config.gpio_sensor_clock,
                                     config.gpio_sensor_data,
                                     config.temp_scale)

        if config.max31856:
            log.info("init MAX31856")
            from max31856 import MAX31856
            software_spi = { 'cs': config.gpio_sensor_cs,
                             'clk': config.gpio_sensor_clock,
                             'do': config.gpio_sensor_data,
                             'di': config.gpio_sensor_di }
            self.thermocouple = MAX31856(tc_type=config.thermocouple_type,
                                         software_spi = software_spi,
                                         units = config.temp_scale,
                                         ac_freq_50hz = config.ac_freq_50hz,
                                         )

    def run(self):
        '''use a moving average of config.temperature_average_samples across the time_step'''
        temps = []
        while True:
            # reset error counter if time is up
            if (time.time() - self.bad_stamp) > (self.time_step * 2):
                if self.bad_count + self.ok_count:
                    self.bad_percent = (self.bad_count / (self.bad_count + self.ok_count)) * 100
                else:
                    self.bad_percent = 0
                self.bad_count = 0
                self.ok_count = 0
                self.bad_stamp = time.time()

            temp = self.thermocouple.get()
            self.noConnection = self.thermocouple.noConnection
            self.shortToGround = self.thermocouple.shortToGround
            self.shortToVCC = self.thermocouple.shortToVCC
            self.unknownError = self.thermocouple.unknownError

            is_bad_value = self.noConnection | self.unknownError
            if not config.ignore_tc_short_errors:
                is_bad_value |= self.shortToGround | self.shortToVCC

            if not is_bad_value:
                temps.append(temp)
                if len(temps) > config.temperature_average_samples:
                    del temps[0]
                self.ok_count += 1

            else:
                log.error("Problem reading temp N/C:%s GND:%s VCC:%s ???:%s" % (self.noConnection,self.shortToGround,self.shortToVCC,self.unknownError))
                self.bad_count += 1

            if len(temps):
                self.temperature = self.get_avg_temp(temps)
            time.sleep(self.sleeptime)

    def get_avg_temp(self, temps, chop=25):
        '''
        strip off chop percent from the beginning and end of the sorted temps
        then return the average of what is left
        '''
        chop = chop / 100
        temps = sorted(temps)
        total = len(temps)
        items = int(total*chop)
        temps = temps[items:total-items]
        return sum(temps) / len(temps)

class Oven(threading.Thread):
    '''parent oven class. this has all the common code
       for either a real or simulated oven'''
    def __init__(self):
        threading.Thread.__init__(self)
        self.daemon = True
        self.temperature = 0
        self.time_step = config.sensor_time_wait
        self.reset()

    def reset(self):
        self.cost = 0
        self.state = "IDLE"
        self.profile = None
        self.start_time = 0
        self.runtime = 0
        self.totaltime = 0
        self.target_temperature = 0
        self.target_pressure = 0
        self.target_power = 0
        self.target_vacuum = 0
        self.heat = 0
        self.pid = PID(ki=config.pid_ki, kd=config.pid_kd, kp=config.pid_kp)

    def run_profile(self, profile, startat=0):
        self.reset()

        if self.board.temp_sensor.noConnection:
            log.info("Refusing to start profile - thermocouple not connected")
            return
        if self.board.temp_sensor.shortToGround:
            log.info("Refusing to start profile - thermocouple short to ground")
            return
        if self.board.temp_sensor.shortToVCC:
            log.info("Refusing to start profile - thermocouple short to VCC")
            return
        if self.board.temp_sensor.unknownError:
            log.info("Refusing to start profile - thermocouple unknown error")
            return

        self.startat = startat * 60
        self.runtime = self.startat
        self.start_time = datetime.datetime.now() - datetime.timedelta(seconds=self.startat)
        self.profile = profile
        self.totaltime = profile.get_duration()
        self.state = "RUNNING"
        log.info("Running schedule %s starting at %d minutes" % (profile.name,startat))
        log.info("Starting")

    def abort_run(self):
        self.reset()
        self.save_automatic_restart_state()

    def kiln_must_catch_up(self):
        '''shift the whole schedule forward in time by one time_step
        to wait for the kiln to catch up'''
        if config.kiln_must_catch_up == True:
            temp = self.board.temp_sensor.temperature + \
                config.thermocouple_offset
            # kiln too cold, wait for it to heat up
            if self.target_temperature - temp > config.pid_control_window:
                log.info("kiln must catch up, too cold, shifting schedule")
                self.start_time = datetime.datetime.now() - datetime.timedelta(milliseconds = self.runtime * 1000)
            # kiln too hot, wait for it to cool down
            if temp - self.target_temperature > config.pid_control_window:
                log.info("kiln must catch up, too hot, shifting schedule")
                self.start_time = datetime.datetime.now() - datetime.timedelta(milliseconds = self.runtime * 1000)

    def update_runtime(self):

        runtime_delta = datetime.datetime.now() - self.start_time
        if runtime_delta.total_seconds() < 0:
            runtime_delta = datetime.timedelta(0)

        self.runtime = runtime_delta.total_seconds()

    def update_target_temp(self):
        self.target_temperature = self.profile.get_target_temperature(self.runtime)

    def update_target_vacuum(self):
        pressure = self.profile.get_target_values(self.runtime)[1]
        power = self.profile.get_target_values(self.runtime)[2]
        self.target_vacuum = pressure * pow(10, power)


    def reset_if_emergency(self):
        '''reset if the temperature is way TOO HOT, or other critical errors detected'''
        if (self.board.temp_sensor.temperature + config.thermocouple_offset >=
            config.emergency_shutoff_temp):
            log.info("emergency!!! temperature too high")
            if config.ignore_temp_too_high == False:
                self.abort_run()

        if self.board.temp_sensor.noConnection:
            log.info("emergency!!! lost connection to thermocouple")
            if config.ignore_lost_connection_tc == False:
                self.abort_run()

        if self.board.temp_sensor.unknownError:
            log.info("emergency!!! unknown thermocouple error")
            if config.ignore_unknown_tc_error == False:
                self.abort_run()

        if self.board.temp_sensor.bad_percent > 30:
            log.info("emergency!!! too many errors in a short period")
            if config.ignore_too_many_tc_errors == False:
                self.abort_run()

    def reset_if_schedule_ended(self):
        if self.runtime > self.totaltime:
            log.info("schedule ended, shutting down")
            log.info("total cost = %s%.2f" % (config.currency_type,self.cost))
            self.abort_run()

    def get_state(self):
        temp = 0
        try:
            temp = self.board.temp_sensor.temperature + config.thermocouple_offset
        except AttributeError as error:
            # this happens at start-up with a simulated oven
            temp = 0
            pass

        state = {
            'cost': self.cost,
            'runtime': self.runtime,
            'temperature': temp,
            'target_temperature': self.target_temperature,
            'target_pressure': self.target_pressure,
            'target_power': self.target_power,
            'target_vacuum': self.target_vacuum,
            'state': self.state,
            'heat': self.heat,
            'totaltime': self.totaltime,
            'kwh_rate': config.kwh_rate,
            'currency_type': config.currency_type,
            'profile': self.profile.name if self.profile else None,
            'pidstats': self.pid.pidstats,
        }
        return state

    def save_state(self):
        with open(config.automatic_restart_state_file, 'w', encoding='utf-8') as f:
            json.dump(self.get_state(), f, ensure_ascii=False, indent=4)

    def state_file_is_old(self):
        '''returns True is state files is older than 15 mins default
                   False if younger
                   True if state file cannot be opened or does not exist
        '''
        if os.path.isfile(config.automatic_restart_state_file):
            state_age = os.path.getmtime(config.automatic_restart_state_file)
            now = time.time()
            minutes = (now - state_age)/60
            if(minutes <= config.automatic_restart_window):
                return False
        return True

    def save_automatic_restart_state(self):
        # only save state if the feature is enabled
        if not config.automatic_restarts == True:
            return False
        self.save_state()

    def should_i_automatic_restart(self):
        # only automatic restart if the feature is enabled
        if not config.automatic_restarts == True:
            return False
        if self.state_file_is_old():
            duplog.info("automatic restart not possible. state file does not exist or is too old.")
            return False

        with open(config.automatic_restart_state_file) as infile:
            d = json.load(infile)
        if d["state"] != "RUNNING":
            duplog.info("automatic restart not possible. state = %s" % (d["state"]))
            return False
        return True

    def automatic_restart(self):
        with open(config.automatic_restart_state_file) as infile: d = json.load(infile)
        startat = d["runtime"]/60
        filename = "%s.json" % (d["profile"])
        profile_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'storage','profiles',filename))

        log.info("automatically restarting profile = %s at minute = %d" % (profile_path,startat))
        with open(profile_path) as infile:
            profile_json = json.dumps(json.load(infile))
        profile = Profile(profile_json)
        self.run_profile(profile,startat=startat)
        self.cost = d["cost"]
        time.sleep(1)
        self.ovenwatcher.record(profile)

    def set_ovenwatcher(self,watcher):
        log.info("ovenwatcher set in oven class")
        self.ovenwatcher = watcher

    def run(self):
        while True:
            if self.state == "IDLE":
                if self.should_i_automatic_restart() == True:
                    self.automatic_restart()
                time.sleep(1)
                continue
            if self.state == "RUNNING":
                self.update_cost()
                self.save_automatic_restart_state()
                self.kiln_must_catch_up()
                self.update_runtime()
                self.update_target_temp()
                self.heat_then_cool()
                self.reset_if_emergency()
                self.reset_if_schedule_ended()

class SimulatedOven(Oven):

    def __init__(self):
        self.board = BoardSimulated()
        self.t_env = config.sim_t_env
        self.c_heat = config.sim_c_heat
        self.c_oven = config.sim_c_oven
        self.p_heat = config.sim_p_heat
        self.R_o_nocool = config.sim_R_o_nocool
        self.R_ho_noair = config.sim_R_ho_noair
        self.R_ho = self.R_ho_noair

        # set temps to the temp of the surrounding environment
        self.t = self.t_env # deg C temp of oven
        self.t_h = self.t_env #deg C temp of heating element

        super().__init__()

        # start thread
        self.start()
        log.info("SimulatedOven started")

    def heating_energy(self,pid):
        # using pid here simulates the element being on for
        # only part of the time_step
        self.Q_h = self.p_heat * self.time_step * pid

    def temp_changes(self):
        #temperature change of heat element by heating
        self.t_h += self.Q_h / self.c_heat

        #energy flux heat_el -> oven
        self.p_ho = (self.t_h - self.t) / self.R_ho

        #temperature change of oven and heating element
        self.t += self.p_ho * self.time_step / self.c_oven
        self.t_h -= self.p_ho * self.time_step / self.c_heat

        #temperature change of oven by cooling to environment
        self.p_env = (self.t - self.t_env) / self.R_o_nocool
        self.t -= self.p_env * self.time_step / self.c_oven
        self.temperature = self.t
        self.board.temp_sensor.temperature = self.t

    def heat_then_cool(self):
        pid = self.pid.compute(self.target_temperature, self.board.temp_sensor.temperature + config.thermocouple_offset)
        pwm_value = int(100 * pid)  # Convertir la valeur du régulateur PID en valeur PWM entre 0 et 100
        # Utilisez la valeur de PWM calculée pour contrôler la chauffe, par exemple :
        self.heating_energy(pid)
        self.temp_changes()

        # self.heat is for the front end to display if the heat is on
        self.heat = 0.0
        if pwm_value > 0:
            self.heat = 1

        log.info("simulation: -> %dW heater: %.0f -> %dW oven: %.0f -> %dW env"            % (int(self.p_heat * pid),
            self.t_h,
            int(self.p_ho),
            self.t,
            int(self.p_env)))

        time_left = self.totaltime - self.runtime

        try:
            log.info("temp=%.2f, target_temperature=%.2f, error=%.2f, pid=%.2f, p=%.2f, i=%.2f, d=%.2f, pwm_value=%.2f, pwm_value=%.2f, run_time=%d, total_time=%d, time_left=%d" %
                (self.pid.pidstats['ispoint'],
                self.pid.pidstats['setpoint'],
                self.pid.pidstats['err'],
                self.pid.pidstats['pid'],
                self.pid.pidstats['p'],
                self.pid.pidstats['i'],
                self.pid.pidstats['d'],
                pwm_value,
                self.runtime,
                self.totaltime,
                time_left))
        except KeyError:
            pass

        # we don't actually spend time heating & cooling during
        # a simulation, so sleep.
        time.sleep(self.time_step)


class RealOven(Oven):

    def __init__(self):
        self.board = Board()
        self.output = Output()
        self.reset()

        # call parent init
        Oven.__init__(self)

        # start thread
        self.start()

    def reset(self):
        super().reset()
        self.output.cool(0)

    def heat_then_cool(self):
        pid = self.pid.compute(self.target_temperature, self.board.temp_sensor.temperature + config.thermocouple_offset)
        pwm_value = int(100 * pid)  # Convertir la valeur du régulateur PID en valeur PWM entre 0 et 100
        # Utilisez la valeur de PWM calculée pour contrôler la chauffe, par exemple :
        self.output.heat(pwm_value)

        # self.heat is for the front end to display if the heat is on
        self.heat = 0.0
        if pwm_value > 0:
            self.heat = 1.0

        if pwm_value > 1:
            self.output.heat(1)
        time_left = self.totaltime - self.runtime
        try:
            log.info("temp=%.2f, target_temperature=%.2f, error=%.2f, pid=%.2f, p=%.2f, i=%.2f, d=%.2f, pwm_value=%.2f, run_time=%d, total_time=%d, time_left=%d" %
                (self.pid.pidstats['ispoint'],
                self.pid.pidstats['setpoint'],
                self.pid.pidstats['err'],
                self.pid.pidstats['pid'],
                self.pid.pidstats['p'],
                self.pid.pidstats['i'],
                self.pid.pidstats['d'],
                pwm_value,
                self.runtime,
                self.totaltime,
                time_left))
        except KeyError:
            pass

class Profile():
    def __init__(self, json_data):
        obj = json.loads(json_data)
        self.name = obj["name"]
        self.data = sorted([(t, (temperature, pressure, power)) for t, temperature, pressure, power in obj["data"]])

    def get_duration(self):
        return max([t for (t, (temperature, pressure, power)) in self.data])

    def get_surrounding_points(self, time):
        if time > self.get_duration():
            return (None, None)

        prev_point = None
        next_point = None

        for i in range(len(self.data)):
            if time < self.data[i][0]:
                prev_point = self.data[i-1]
                next_point = self.data[i]
                break

        return (prev_point, next_point)

    def get_target_temperature(self, time):
        if time > self.get_duration():
            return 0

        (prev_point, next_point) = self.get_surrounding_points(time)
        incl_temperature = (next_point[1][0] - prev_point[1][0]) / (next_point[0] - prev_point[0])
        temperature = prev_point[1][0] + (time - prev_point[0]) * incl_temperature
        return temperature

    def get_target_pressure(self, time):
        if time > self.get_duration():
            return 0

        (prev_point, next_point) = self.get_surrounding_points(time)
        incl_pressure = (next_point[1][1] - prev_point[1][1]) / (next_point[0] - prev_point[0])
        pressure = prev_point[1][1] + (time - prev_point[0]) * incl_pressure
        return pressure

    def get_target_power(self, time):
        if time > self.get_duration():
            return 0

        (prev_point, next_point) = self.get_surrounding_points(time)
        incl_power = (next_point[1][2] - prev_point[1][2]) / (next_point[0] - prev_point[0])
        power = prev_point[1][2] + (time - prev_point[0]) * incl_power
        return power





class PID():

    def __init__(self, ki=1, kp=1, kd=1):
        self.ki = ki
        self.kp = kp
        self.kd = kd
        self.lastNow = datetime.datetime.now()
        self.iterm = 0
        self.lastErr = 0
        self.pidstats = {}

    # FIX - this was using a really small window where the PID control
    # takes effect from -1 to 1. I changed this to various numbers and
    # settled on -50 to 50 and then divide by 50 at the end. This results
    # in a larger PID control window and much more accurate control...
    # instead of what used to be binary on/off control.
    def compute(self, setpoint, ispoint):
        now = datetime.datetime.now()
        timeDelta = (now - self.lastNow).total_seconds()

        window_size = 100

        error = float(setpoint - ispoint)

        # this removes the need for config.stop_integral_windup
        # it turns the controller into a binary on/off switch
        # any time it's outside the window defined by
        # config.pid_control_window
        icomp = 0
        output = 0
        out4logs = 0
        dErr = 0
        if error < (-1 * config.pid_control_window):
            log.info("kiln outside pid control window, max cooling")
            output = 0
            # it is possible to set self.iterm=0 here and also below
            # but I dont think its needed
        elif error > (1 * config.pid_control_window):
            log.info("kiln outside pid control window, max heating")
            output = 1
        else:
            icomp = (error * timeDelta * (1/self.ki))
            self.iterm += (error * timeDelta * (1/self.ki))
            dErr = (error - self.lastErr) / timeDelta
            output = self.kp * error + self.iterm + self.kd * dErr
            output = sorted([-1 * window_size, output, window_size])[1]
            out4logs = output
            output = float(output / window_size)
            
        self.lastErr = error
        self.lastNow = now

        # no active cooling
        if output < 0:
            output = 0

        self.pidstats = {
            'time': time.mktime(now.timetuple()),
            'timeDelta': timeDelta,
            'setpoint': setpoint,
            'ispoint': ispoint,
            'err': error,
            'errDelta': dErr,
            'p': self.kp * error,
            'i': self.iterm,
            'd': self.kd * dErr,
            'kp': self.kp,
            'ki': self.ki,
            'kd': self.kd,
            'pid': out4logs,
            'out': output,
        }

        return output
