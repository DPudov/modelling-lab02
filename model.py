from math import pi

import scipy.interpolate as interpolate
import numpy as np

Mode = 1
LengthBetweenElectrodesDefault = 12.0
LampRadiusDefault = 0.35
InductanceDefault = 187 * 1e-6
CapacityDefault = 268 * 1e-6
ResistanceDefault = 0.25
VoltageDefault = 1400.0
CurrentDefault = 0.5
EdgeTemperature = 2000.0
TimeDefault = 0.0
TimeEndDefault = 600. * 1e-6
TimeStep = 1.e-6
UsingGas = True


def log_interp1d(xx, yy):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(logx, logy, fill_value='extrapolate')
    log_interp = lambda dd: np.power(10.0, lin_interp(np.log10(dd)))
    return log_interp


z_axis = np.linspace(0., 1., 51)


def simps(f, a, b, N, current):
    dx = (b - a) / N
    x = z_axis
    y = f(x, current)
    S = dx / 3 * np.sum(y[0:-1:2] + 4 * y[1::2] + y[2::2])
    return S


currents = [0.5, 1., 5., 10., 50., 200., 400., 800., 1200.]
temperatures_for_currents = [6400, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010]
m_powers = [0.4, 0.55, 1.7, 3, 11, 32, 40, 41, 39]
current_temperature_interp = log_interp1d(currents, temperatures_for_currents)
current_m_interp = log_interp1d(currents, m_powers)

temperatures_for_resistance = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
resistances = [0.031, 0.27, 2.05, 6.06, 12.0, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]
temperature_resistance_interp = log_interp1d(temperatures_for_resistance, resistances)


def current_temperature_dataset():
    return list(zip(currents, temperatures_for_currents, m_powers))


def temperature_resistance_dataset():
    return list(zip(temperatures_for_resistance, resistances))


def length_between_electrodes():
    return LengthBetweenElectrodesDefault


def lamp_radius():
    return LampRadiusDefault


def resistance_static():
    return ResistanceDefault


def inductance():
    return InductanceDefault


def capacity():
    return CapacityDefault


def time_start():
    return TimeDefault


def time_end():
    return TimeEndDefault


def temperature(z, current):
    m = current_m_interp(abs(current))
    T0 = current_temperature_interp(abs(current))
    return T0 + (EdgeTemperature - T0) * (z ** m)


def sigma_temperature(temperature):
    return temperature_resistance_interp(temperature)


def sigma(z, current):
    return sigma_temperature(temperature(z, current))


def resistance(z, current):
    res = []
    for item in z:
        res.append(sigma(item, current) * item)
    return res


def lamp_resistance(current):
    if UsingGas:
        length = length_between_electrodes()
        radius = lamp_radius()
        return length / (2. * pi * radius * radius * simps(resistance, 0., 1., 50, current))
    else:
        return 0.


def initial_vector():
    return [CurrentDefault,
            VoltageDefault]


def vector_function(current, voltage):
    lamp = lamp_resistance(current)
    return [(voltage - current * (resistance_static() + lamp)) / inductance(),
            -current / capacity(), lamp]


def runge_kutta_2(alpha=1.):
    res_timeline = []
    res_currents = []
    res_voltages = []
    res_resistances = []
    res_temperatures = []
    t = time_start()
    prev = initial_vector()

    res_timeline.append(t)
    res_currents.append(prev[0])
    res_voltages.append(prev[1])
    res_resistances.append(lamp_resistance(prev[0]))
    res_temperatures.append(temperature(0, prev[0]))
    while t < time_end():
        v0 = vector_function(prev[0], prev[1])
        v1 = vector_function(prev[0] + TimeStep / (2.0 * alpha),
                             prev[1] + TimeStep / (2.0 * alpha * v0[1]))
        prev[0] += TimeStep * ((1. - alpha) * v0[0] + alpha * v1[0])
        prev[1] += TimeStep * ((1. - alpha) * v0[1] + alpha * v1[1])
        t += TimeStep

        res_timeline.append(t)
        res_currents.append(prev[0])
        res_voltages.append(prev[1])
        res_resistances.append(lamp_resistance(prev[0]))
        res_temperatures.append(temperature(0, prev[0]))
        print(t)
    return res_timeline, res_currents, res_voltages, res_resistances, res_temperatures


def runge_kutta_4():
    res_timeline = []
    res_currents = []
    res_voltages = []
    res_resistances = []
    res_temperatures = []
    t = time_start()
    prev = initial_vector()

    res_timeline.append(t)
    res_currents.append(prev[0])
    res_voltages.append(prev[1])
    res_resistances.append(lamp_resistance(prev[0]))
    res_temperatures.append(temperature(0, prev[0]))
    while t < time_end():
        v1 = vector_function(prev[0], prev[1])
        k1 = v1[0]
        q1 = v1[1]

        v2 = vector_function(prev[0] + k1 * TimeStep / 2.0, prev[1] + q1 * TimeStep / 2.0)
        k2 = v2[0]
        q2 = v2[1]

        v3 = vector_function(prev[0] + k2 * TimeStep / 2.0, prev[1] + q2 * TimeStep / 2.0)
        k3 = v3[0]
        q3 = v3[1]

        v4 = vector_function(prev[0] + k3 * TimeStep, prev[1] + q3 * TimeStep)
        k4 = v4[0]
        q4 = v4[1]

        prev[0] += TimeStep * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        prev[1] += TimeStep * (q1 + 2 * q2 + 2 * q3 + q4) / 6.0
        t += TimeStep

        res_timeline.append(t)
        res_currents.append(prev[0])
        res_voltages.append(prev[1])
        res_resistances.append(lamp_resistance(prev[0]))
        res_temperatures.append(temperature(0, prev[0]))
        print(t)
    return res_timeline, res_currents, res_voltages, res_resistances, res_temperatures
