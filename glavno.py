import math
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import chebfit, chebval
# import genetski_algoritam
import prvi
# import drugi


# Ova funkcija prevodi niz dekadnih brojeva u njihov binarni zapis

# def binary(num):
#    # Struct can provide us with the float packed into bytes. The '!' ensures that
#    # it's in network byte order (big-endian) and the 'f' says that it should be
#    # packed as a float. Alternatively, for double-precision, you could use 'd'.
#    packed = struct.pack('!f', num)
#    #print 'Packed: %s' % repr(packed)
#
#    # For each character in the returned string, we'll turn it into its corresponding
#    # integer code point
#    #
#    # [62, 163, 215, 10] = [ord(c) for c in '>\xa3\xd7\n']
#    integers = [ord(c) for c in packed]
#    #print 'Integers: %s' % integers
#
#    # For each integer, we'll convert it to its binary representation.
#    binaries = [bin(i) for i in integers]
#    #print 'Binaries: %s' % binaries
#
#    # Now strip off the '0b' from each of these
#    stripped_binaries = [s.replace('0b', '') for s in binaries]
#    #print 'Stripped: %s' % stripped_binaries
#
#    # Pad each byte's binary representation's with 0's to make sure it has all 8 bits:
#    #
#    # ['00111110', '10100011', '11010111', '00001010']
#    padded = [s.rjust(8, '0') for s in stripped_binaries]
#    #print 'Padded: %s' % padded
#
#    # At this point, we have each of the bytes for the network byte ordered float
#    # in an array as binary strings. Now we just concatenate them to get the total
#    # representation of the float:
#    return ''.join(padded)


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


G = 1.32712440018e20
info = [(0.38709927,   0.00000037, 0.20563593,  0.00001906, 252.25032350, 149472.67411175,  77.45779628),
        (0.72333566,   0.00000390, 0.00677672, -0.00004107, 181.97909950,  58517.81538729, 131.60246718),
        (1.00000261,   0.00000562, 0.01671123, -0.00004392, 100.46457166,  35999.37244981, 102.93768193),
        (1.52371034,   0.00001847, 0.09339410,  0.00007882,  -4.55343205,  19140.30268499, -23.94362959),
        (5.20288700,  -0.00011607, 0.04838624, -0.00013253,  34.39644051,   3034.74612775,  14.72847983),
        (9.53667594,  -0.00125060, 0.05386179, -0.00050991,  49.95424423,   1222.49362201,  92.59887831),
        (19.18916464, -0.00196176, 0.04725744, -0.00004397, 313.23810451,    428.48202785, 170.95427630),
        (30.06992276,  0.00026291, 0.00859048,  0.00005105, -55.12002969,    218.45945325,  44.96476227),
        (39.48211675, -0.00031596, 0.24882730,  0.00005170, 238.92903833,    145.20780515, 224.06891629)]

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
    koeficijenti = np.array([chebfit(x_osa(br_segm), uglovi, chebdeg) for uglovi in uglovi_matrica])
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


coef, uklj = inicijalizacija(5, 10, 8)
matr = otpakivanje(coef, uklj)
coef, uklj = pakovanje(matr, 8)
print(np.array([chebval(x_osa(10), cheb) for cheb in coef]))
