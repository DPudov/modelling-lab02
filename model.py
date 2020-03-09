from math import pi

import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import numpy as np

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


def log_interp1d(xx, yy):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(logx, logy, fill_value='extrapolate')
    log_interp = lambda dd: np.power(10.0, lin_interp(np.log10(dd)))
    return log_interp


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
    return sigma(z, current) * z


def lamp_resistance(current):
    length = length_between_electrodes()
    radius = lamp_radius()

    return length / (2. * pi * radius * radius * integrate.quad(resistance, 0., 1., args=(current))[0])


def initial_vector():
    return [CurrentDefault,
            VoltageDefault]


def vector_function(time, current, voltage):
    # print("Current: ", current)
    # print("Voltage ", voltage)
    return [(voltage - current * (resistance_static() + lamp_resistance(current))) / inductance(),
            -current / capacity()]


def runge_kutta_4():
    res_timeline = []
    res_currents = []
    res_voltages = []
    t = time_start()
    prev = initial_vector()

    res_timeline.append(t)
    res_currents.append(prev[0])
    res_voltages.append(prev[1])
    while t < time_end():
        v1 = vector_function(t, prev[0], prev[1])
        k1 = v1[0]
        q1 = v1[1]

        v2 = vector_function(t + TimeStep / 2.0, prev[0] + k1 * TimeStep / 2.0, prev[1] + q1 * TimeStep / 2.0)
        k2 = v2[0]
        q2 = v2[1]

        v3 = vector_function(t + TimeStep / 2.0, prev[0] + k2 * TimeStep / 2.0, prev[1] + q2 * TimeStep / 2.0)
        k3 = v3[0]
        q3 = v3[1]

        v4 = vector_function(t + TimeStep, prev[0] + k3 * TimeStep, prev[1] + q3 * TimeStep)
        k4 = v4[0]
        q4 = v4[1]

        prev[0] += TimeStep * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        prev[1] += TimeStep * (q1 + 2 * q2 + 2 * q3 + q4) / 6.0
        t += TimeStep

        res_timeline.append(t)
        res_currents.append(prev[0])
        res_voltages.append(prev[1])
        print(t)
    return res_timeline, res_currents, res_voltages


r = runge_kutta_4()
plt.subplot(2, 1, 1)
plt.plot(r[0], r[1])
plt.ylabel("Current")

plt.subplot(2, 1, 2)
plt.plot(r[0], r[2])
plt.ylabel("Voltage")
plt.xlabel("time")
plt.show()
