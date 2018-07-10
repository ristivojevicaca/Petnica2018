import numpy as np
import math
import podaci
import solar


def modulo(vector):
    return math.sqrt(np.sum(vector ** 2))


def jed_vec(ugao):
    return np.array((math.cos(ugao), math.sin(ugao)))


def izracunaj(r_, v_, t, brod, ugao, k=np.zeros((3, 2)), step=0.0):
    dry_mass, fuel_mass = brod
    mass = (dry_mass, fuel_mass)
    r_, v_, mass = np.array((r_, v_, mass)) + step * k
    a_ = solar.thrust(modulo(r_), t + step)/np.sum(mass) * jed_vec(ugao) - podaci.grav_par[0] * r_ / (modulo(r_) ** 3)
    delta_mass = -solar.flow_rate(modulo(r_), t + step)
    return np.array((v_, a_, np.array((0, delta_mass))))


def runge_kuta4(r_, v_, t, brod, ugao, step):
    dry_mass, fuel_mass = brod
    k1 = step * izracunaj(r_, v_, t, brod, ugao)
    k2 = step * izracunaj(r_, v_, t, brod, ugao, k1/2, step/2)
    k3 = step * izracunaj(r_, v_, t, brod, ugao, k2/2, step/2)
    k4 = step * izracunaj(r_, v_, t, brod, ugao, k3, step)
    k = k1 + 2 * k2 + 2 * k3 + k4
    return np.array((r_, v_, np.array((dry_mass, fuel_mass)))) + k/6


def simulacija(r_, v_, brod, uglovi, n):  # simulira kretanje tela samo pocetnom brzinom u gravitacionom polju sunca
    # m = m_ukupna
    _r = np.empty((n, 2))
    _v = np.empty((n, 2))
    _time = np.empty(n)
    r = modulo(r_)
    v = modulo(v_)
    time = 0.0
    limit = 0.007
    # prev_sol = 0
    for i in range(n):
        dry_mass, fuel_mass = brod
        if fuel_mass <= 0:
            break
        a_ = -podaci.grav_par[0] * r_ / (r ** 3)
        a_ = a_ + solar.thrust(r, time) / (dry_mass + fuel_mass) * jed_vec(uglovi[i])
        a = modulo(a_)
        if a == 0:
            step = 3600*12
        else:
            step = math.ceil((v/a)*limit)

        if solar.flow_rate(r, time) * step > fuel_mass:
            step = fuel_mass / solar.flow_rate(r, time)

        (r_, v_, brod) = runge_kuta4(r_, v_, time, brod, uglovi[i], step)
        r = modulo(r_)
        v = modulo(v_)
        _r[i] = r_
        _v[i] = v_
        _time[i] = time+step
        time = _time[i]
    return _r, _v
