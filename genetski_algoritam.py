from numpy import random


# populacija je lista koja cuva podatke za sve jedinke u populaciji.
# svaka jedinka ima binarni kod (predstavljen nizom 0 i 1) i
# vrednost fitnes funkcije
# elitizam bira procenat najboljih u populaciji i njima menja isti procenat najgorih
def elitizam(populacija, n, p_elit):
    m = round(p_elit * n)
    populacija = sorted(populacija, key=lambda x: x[1], reverse=True)
    populacija = populacija[:n - m]
    populacija.extend(populacija[:m])
    return populacija


def uparivanje(lista, n):  # deli elemente liste na parove uzastopnih elem
    return [lista[i:i+2] for i in range(0, n-1, 2)]


# turnir je vrsenje prirodne selekcije. u svakom paru koji se takmici pobedjuje
# najbolji. imas objasnjen geneticki algoritam u predlogu projekta.
def tournament(lista_parova, n):
    pobednici = []
    for i in range(n):
        if lista_parova[i][0][1] > lista_parova[i][1][1]:
            pobednici.append(lista_parova[i][0])
        else:
            pobednici.append(lista_parova[i][1])
    return pobednici


# vrsi selekciju. promesa listu nasumicno, daje da se na osnovu nje biraju parovi i onda
#  radi turnir da odredi najbolje
def selekcija(populacija, n):

    random.shuffle(populacija)

    parovi = uparivanje(populacija, n)

    tournament1 = tournament(parovi, n // 2)
    # print(tournament1)
    random.shuffle(populacija)

    parovi = uparivanje(populacija, n)

    tournament2 = tournament(parovi, n // 2)  # turnir se radi dvaput i onda se rezultati turnira spajaju
    # print(tournament2)
    return tournament1 + tournament2


def kros_over(kod1, kod2, n):  # svaki bit novog koda generise izborom bita jednog od roditelja.
    novi_kod = []
    for i in range(n):
        if random.random() < 0.5:
            novi_kod.append(kod1[i])
        else:
            novi_kod.append(kod2[i])
    return novi_kod


# dvaput ukrsta svaki par da bi populacija dece
# imala isti broj clanova kao populacija roditelja
def ukrstanje_parova(parovi, n):
    potomstvo = []
    m = len(parovi[0][0][0])
    for i in range(n):
        potomstvo.append([kros_over(parovi[i][0][0], parovi[i][1][0], m)])
        potomstvo.append([kros_over(parovi[i][0][0], parovi[i][1][0], m)])
    return potomstvo


def ukrstanje(populacija, n):  # ukrsta parove jedinki - svojevrsno razmnozavanje
    random.shuffle(populacija)
    parovi = uparivanje(populacija, n)
    return ukrstanje_parova(parovi, n//2)


def izvrni_bit(string, indeks):  # mutira bitove
    if string[indeks] == 0:
        string[indeks] = 1
    else:
        string[indeks] = 0
    return string


def mutacija(populacija, n, p_mut):  # mutira odredjen procenat bitova u populaciji
    m = len(populacija[0][0])
    for i in range(n):
        br_prom = random.binomial(m, p_mut)
        if br_prom != 0:
            indeksi = random.choice(m, br_prom)
            for j in indeksi:
                izvrni_bit(populacija[i][0], j)


def genetski_algoritam(populacija, n, p_elit, p_mut):
    populacija = elitizam(populacija, n, p_elit)
    populacija = selekcija(populacija, n)
    populacija = ukrstanje(populacija, n)
    mutacija(populacija, n, p_mut)
    return populacija
