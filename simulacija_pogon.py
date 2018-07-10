import numpy as np
import math
import podaci
import motor
# import time as tm


def modulo(vector):
    return math.sqrt(np.sum(vector ** 2))


def jed_vec(ugao):
    return np.array((math.cos(ugao), math.sin(ugao)))


def izracunaj(r_, v_, t, brod, ugao, motor_uklj, k=np.zeros((3, 2)), step=0.0):
    r_, v_, brod = np.array((r_, v_, brod)) + k
    r = modulo(r_)
    a_ = - podaci.grav_par[0] * r_ / (r ** 3)
    delta_mass = 0
    if motor_uklj:
        a_ = a_ + motor.thrust(r, t + step) / np.sum(brod) * jed_vec(ugao)
        delta_mass = -motor.flow_rate(r, t + step)
    return np.array((v_, a_, np.array((0, delta_mass))))


def runge_kuta4(r_, v_, t, brod, ugao, step, motor_uklj):
    k1 = step * izracunaj(r_, v_, t, brod, ugao, motor_uklj)
    k2 = step * izracunaj(r_, v_, t, brod, ugao, motor_uklj, k1/2, step/2)
    k3 = step * izracunaj(r_, v_, t, brod, ugao, motor_uklj, k2/2, step/2)
    k4 = step * izracunaj(r_, v_, t, brod, ugao, motor_uklj, k3, step)
    k = k1 + 2 * k2 + 2 * k3 + k4
    return np.array((r_, v_, brod)) + k/6


def simulacija(r_, v_, brod, uglovi, y_max):  # simulira kretanje tela samo pocetnom brzinom u gravitacionom polju sunca
    # m = m_ukupna
    motor_uklj = True
    n = len(uglovi)
    t_max = y_max * 365.25 * 24 * 3600
    # print(t_max)
    len_nizova = int(t_max / (6*3600))
    _r = np.zeros((len_nizova, 2))
    _v = np.zeros((len_nizova, 2))
    _time = np.zeros(len_nizova)
    _step = np.zeros(len_nizova)
    r = modulo(r_)
    v = modulo(v_)
    time = 0.0
    limit = 0.007
    i = 0
    # prev_sol = 0
    while time < t_max:
        # start = tm.process_time()
        ind = math.floor(time / t_max * n)
        dry_mass, fuel_mass = brod
        if fuel_mass <= 0:
            motor_uklj = False
        a_ = -podaci.grav_par[0] * r_ / (r ** 3)
        if motor_uklj:
            a_ = a_ + motor.thrust(r, time) / (dry_mass + fuel_mass) * jed_vec(uglovi[ind])
        a = modulo(a_)
        if a == 0:
            step = 3600*12
        else:
            step = math.ceil((v/a)*limit)

        if motor.flow_rate(r, time) * step > fuel_mass and motor_uklj:
            step = math.ceil(fuel_mass / motor.flow_rate(r, time))
        
        if math.floor((time + step) / t_max * n) != ind:
            step = math.ceil((ind+1) * t_max/n - time)
        
        (r_, v_, brod) = runge_kuta4(r_, v_, time, brod, uglovi[ind], step, motor_uklj)

        r = modulo(r_)
        v = modulo(v_)
        if i == len_nizova:
            _r = np.append(_r, np.empty(100))
            _v = np.append(_v, np.empty(100))
            _time = np.append(_time, np.empty(100))
            _step = np.append(_step, np.empty(100))
        _r[i] = r_
        _v[i] = v_
        _time[i] = time+step
        _step[i] = step
        time = _time[i]
        # print(time, i)
        i = i + 1
        # print(tm.process_time() - start)
    return _r[:i-1], _v[:i-1], _step[:i-1]
