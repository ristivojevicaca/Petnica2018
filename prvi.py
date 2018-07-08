import math
import numpy as np
from numpy import array as u_niz
import matplotlib.pyplot as plt


def simulacija(x, y, vx, vy, gm, n, step):  # simulira kretanje tela samo pocetnom brzinom u gravitacionom polju sunca
    _x = np.empty(n)
    _y = np.empty(n)
    _t = np.empty(n)
    data = (x, y, vx, vy)
    for i in range(n):
        data = runge_kuta4(data, gm, step)
        _x[i] = data[0]
        _y[i] = data[1]
        _t[i] = i * step
        # Teta[i] = ugao_teta(x_sim,y_sim);
    # plt.plot (u_niz(T)/(3600*24),u_niz(X)/1e9,label='x');
    # plt.plot(u_niz(T)/(3600*24),u_niz(Y)/1e9,label='y');
    # pokusaj da se ugao pod kojim se nalazi brod aproksimira polinomom - NE RADI
    stepen = 5
    koef_x = np.polyfit(_t, _x, stepen)
    koef_y = np.polyfit(_t, _y, stepen)
    # f = splrep(T[0:n-1:128],Teta[0:n-1:128],k=3);
    polinom_x = np.poly1d(koef_x)
    polinom_y = np.poly1d(koef_y)
    vrednosti_x = list(map(polinom_x, _t))
    vrednosti_y = list(map(polinom_y, _t))
    plt.figure()
    plt.plot(u_niz(_t), u_niz(_x), 'r')
    plt.plot(_t, vrednosti_x, 'r--')
    plt.plot(u_niz(_t), u_niz(_y), 'b')
    plt.plot(_t, vrednosti_y, 'b--')
    plt.figure()
    plt.plot(_x, _y)
    plt.axis('equal')

    # vrednosti = splev(T,f)
    # plt.figure();
    # plt.plot(Teta,'r--',ms=2,linewidth=1);
    # plt.plot(vrednosti,'k.',linewidth=0.1,ms=2);
    # plt.show()
    # print(Teta)
    # plt.legend()
    print(100 * ((-gm / math.sqrt(data[0] ** 2 + data[1] ** 2) + 0.5 * (data[2] ** 2 + data[3] ** 2)) / (
            -gm / math.sqrt(x ** 2 + y ** 2)
            + 0.5 * (vx * vx + vy * vy)) - 1), '%')
    print(n * step)


def izracunaj(gm, r, v0, k, step):  # k je niz dva koef.- jedan za izvod polozaja drugi za izvod brzine
    r = r + k[1] * step  # r je niz sa dve elem - x i y koord
    r_ = np.sqrt(r[0] * r[0] + r[1] * r[1])
    a = -1 * (r / r_) * gm / (r_ * r_)
    v = v0 + k[0] * step
    # print(r,R,a,v,'End\n')
    # print('a ', a,'shape', np.shape (a),'b ',v,'shape ', np.shape(v));
    return u_niz([a, v])


def runge_kuta4(data, gm, step):  # runge kuta 4. reda. mislim da i njegovo objasnjeneje imas na driveu.
    x, y, vx, vy = data
    r = u_niz([x, y])
    v = u_niz([vx, vy])
    k1 = izracunaj(gm, r, v, u_niz([[0, 0], [0, 0]]), 0)
    k2 = izracunaj(gm, r, v, k1, step / 2)
    k3 = izracunaj(gm, r, v, k2, step / 2)
    k4 = izracunaj(gm, r, v, k3, step)
    k = k1 + 2 * k2 + 2 * k3 + k4
    # print(k[0][0],k[0][0]);
    # print('rk: ',[r+(step/6)*k[1],v+(step/6)*k[0]]);
    return [r + (step / 6) * k[1], v + (step / 6) * k[0]]
