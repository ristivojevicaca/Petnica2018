from numpy.polynomial.chebyshev import chebfit, chebval
import matplotlib.pyplot as plt
from math import pi
import numpy as np
import math
import time
from podaci import r0_, v0_, n
import simulacija_pogon
import prvi
# import podaci
# import genetski_algoritam
# import drugi


def to_binary(niz_num):
    duz = len(niz_num)
    kod = []
    for j in range(duz):
            bitan_deo = [int(digit) for digit in bin(niz_num[j])[2:]]
            bitan_deo = [0 for _ in range(16 - len(bitan_deo))] + bitan_deo
            kod = kod + bitan_deo
    return kod

# from bitstring import BitArray
# >>> bitlist=[1,0,0,0,0,0,0,0]
# >>> b = BitArray(bitlist)
# >>> b.uint


def to_uint16(nizbin):  # ova funkcija prevodi niz 1 i 0 u nizniz dekadnih brojeva. inverzna funkcija prosloj
    duz = len(nizbin)//16
    # print(duz)
    niz_num = []
    for j in range(duz):
        niz = nizbin[16*j:16*j+16]
        niz.reverse()
        suma = 0
        stepen = 1
        for i in range(16):
            suma = suma + stepen * niz[i]
            stepen = stepen << 1
        niz_num.append(suma)
    return niz_num

# print(toBinary(0.15625))
# gornje dve funkcije su potrebne jer se 'genetski kod' putanja zapisuje u obliku dugackog binarnog broja,
# koje je u stvari spojeni niz zapisa uglova pod kojima radi motor

# podaci potrebni za racunanje pozicije planete po formulama. imas formule u jednom pdf - u na drajvu.
# redosled podataka je:
# masa - u jedinicama
# a - velika poluosa i njena promena po veku
# e - ekscentricitet i promena po veku
# L - srednja longituda i promena po veku
# malo teta - longtituda perhiela i promena po veku
# Teta - longtitude of ascending node i promena po veku
# Sunce
# Merkur
# Venera
# Zemlja
# Mars
# Jupiter
# Saturn
# Uran
# Neptun
# podaciPlaneta = np.array
# konst.

# testiranje - rastojanje 150e9, brzina 29740
# def polozaj_planete (podaciPlaneta,i):

# simulacija kretanja bez motora


def ugao_teta(x, y):  # odredjuje ugao u odnosu na sunce pod kojim je svemirski brod. nije bitan za samu simulaciju.
    if y >= 0:
        return math.acos(x/math.sqrt(x*x+y*y))
    else:
        return 2*math.pi - math.acos(x/math.sqrt(x*x+y*y))


# RungeKutta svemirski brod u pokretu sa silom F
# testiranje - rastojanje 150e9, brzina 29740
# def polozaj_planete (podaciPlaneta,i):


def fitness(min_d_mars):  # ideja za jednostavnu fitnes funkciju - to je minimalno rastojnje od marsa
    return min_d_mars - 200


def crtanje_planeta():
    broj_tacaka = 365 * 5
    broj_planeta = 4
    t = np.linspace(0, 365 * 5, broj_tacaka)
    x = [[] for _ in range(broj_planeta)]
    y = [[] for _ in range(broj_planeta)]
    e_prev = [2.0 for _ in range(broj_planeta)]
    indeksi_planeta = [0, 1, 2, 3]  # indeksi planeta koje crtamo
    for i in range(broj_tacaka):
        for j in range(broj_planeta):  # koliko planeta crtamo
            rez = prvi.polozaj_planeta(indeksi_planeta[j], t[i], e_prev[j])
            x[j].append(rez[0])
            y[j].append(rez[1])
            e_prev[j] = rez[2]

    plt.plot(x[0], y[0], 'g', linewidth=0.8)
    plt.plot(x[1], y[1], 'y', linewidth=0.8)
    plt.plot(x[2], y[2], 'b', linewidth=0.8)
    plt.plot(x[3], y[3], 'r', linewidth=0.8)

    # plt.setp(putanje,markersize=2.5)
    plt.plot(0.0, 0.0, 'k*', markersize=7)
    plt.axis('scaled')
    plt.title('Putanje Merkura, Venere, Zemlje i Marsa oko Sunca')
    plt.xlabel('x koord [astronomska jedinica]')
    plt.ylabel('y koord [astronomska jedinica]')
    # plt.show()


def x_osa(a):
    return np.arange(a)


def inicijalizacija(n, br_segm, chebdeg):
    uglovi_matrica = np.multiply(np.random.random_sample((n, br_segm)), 2 * pi)
    print(uglovi_matrica)
    koeficijenti = np.array([chebfit(x_osa(br_segm), _, chebdeg) for _ in uglovi_matrica])
    temp = np.array([chebval(x_osa(br_segm), cheb) for cheb in koeficijenti])
    print(np.divide(temp, uglovi_matrica))
    snaga = np.random.random_integers(0, 1, (n, br_segm))
    return koeficijenti, snaga


def otpakivanje(koeficijenti, snaga):
    koef_bitovi = np.unpackbits(koeficijenti.view(np.uint8), axis=1)
    matrica = np.concatenate((koef_bitovi, snaga), axis=1)
    return matrica


def pakovanje(matrica, chebdeg):
    koef_bitovi = matrica[:, :(chebdeg + 1) * 64]
    snaga = matrica[:, (chebdeg + 1) * 64:]
    koeficijenti = np.packbits(koef_bitovi, axis=1).view(np.float_)
    return koeficijenti, snaga

# def kodiranje(koeficijenti, snaga, n, br_segm):


# uglovi = np.zeros(1000)
uglovi = np.ones(n) * pi

start = time.process_time()
_r, _v, _step = simulacija_pogon.simulacija(r0_, v0_, (4e2, 0), uglovi, n)
# _b = simulacija_pogon.simulacija(r0_, v0_, (4e2, 20), np.zeros(1000), 1000)
# _c = prvi.simulacija(150e9, 0, 0, 29780, podaci.grav_par[0], 1000)
print(time.process_time() - start)
plt.plot(_r[:, 0], _r[:, 1])
# plt.plot(_b[:, 0], _b[:, 1])
# plt.plot(_c[0], _c[1])
plt.axis('equal')
plt.show()
