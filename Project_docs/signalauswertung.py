#signalauswertung: aus einem sondensignal u(t) den periodischen endzustand bestimmen.
#
#steht bewusst neben Animation.py und nicht in Validierung/: beide auswertungen - das
#bild aus Animation.py und die benchmark-serien - muessen dieselbe strouhal-zahl
#liefern, sonst vergleicht man zwei verschiedene messvorschriften. Validierung/
#benchmark.py importiert die funktionen hier und prueft sie in seinem selbsttest.

import math

import numpy as np

#zeitfenster am laufende, ueber das gemittelt wird, wenn es keine perioden gibt
#(stationaerer nachlauf bei kleinem Re). auch fuer die restschwankung
FENSTER_ENDE = 10.0


def nulldurchgaenge(t, s):
    """indizes k mit s[k] < 0 <= s[k+1] und die linear interpolierten zeitpunkte"""
    k = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    t_null = t[k] - s[k] * (t[k + 1] - t[k]) / (s[k + 1] - s[k])
    return k, t_null


def leere_auswertung(status):
    return dict(status=status, St=math.nan, periode=math.nan, periode_streuung=math.nan,
                amplitude=math.nan, mittelwert=math.nan, t_einsatz=math.nan,
                t_fenster_start=math.nan, t_fenster_ende=math.nan, n_perioden=0)


def werte_signal_aus(t, u, D, U_inf, tol_A=0.02, tol_T=0.01, min_perioden=4, min_amplitude=1e-3):
    """
    bestimmt aus dem sondensignal den periodischen endzustand:
    - perioden T_k zwischen aufwaerts-nulldurchgaengen, amplitude A_k = (max - min)/2 je periode
    - auswertefenster: die letzten perioden, deren A_k und T_k um hoechstens tol_A bzw. tol_T
      vom median der letzten drei perioden abweichen (der anlauf faellt so heraus)
    - St = D / (U_inf * mittlere periode im fenster)
    - t_einsatz: erster zeitpunkt, an dem |u - mittelwert| die halbe endamplitude erreicht
    zweimal durchlaufen: im zweiten durchgang werden die nulldurchgaenge um den mittelwert
    des fensters bestimmt.
    """
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    letztes_viertel = t >= t[0] + 0.75 * (t[-1] - t[0])
    if 0.5 * np.ptp(u[letztes_viertel]) < min_amplitude * U_inf:
        return leere_auswertung("keine_abloesung")

    mittel = 0.0
    for durchgang in range(2):
        k, t_null = nulldurchgaenge(t, u - mittel)
        if len(t_null) < 4:
            return leere_auswertung("nicht_periodisch")

        T = np.diff(t_null)
        A = np.array([0.5 * np.ptp(u[k[c] + 1:k[c + 1] + 1]) for c in range(len(T))])
        A_ref = np.median(A[-3:])
        T_ref = np.median(T[-3:])

        c_start = len(T)
        while (c_start > 0 and abs(A[c_start - 1] / A_ref - 1.0) <= tol_A
               and abs(T[c_start - 1] / T_ref - 1.0) <= tol_T):
            c_start -= 1
        if c_start == len(T):
            return leere_auswertung("nicht_periodisch")

        t_a, t_b = t_null[c_start], t_null[-1]
        im_fenster = (t >= t_a) & (t <= t_b)
        mittel = np.trapezoid(u[im_fenster], t[im_fenster]) / (t[im_fenster][-1] - t[im_fenster][0])

    T_fenster = T[c_start:]
    A_fenster = A[c_start:]
    amplitude = float(np.mean(A_fenster))
    einsatz = np.nonzero(np.abs(u - mittel) >= 0.5 * amplitude)[0]

    ergebnis = leere_auswertung("periodisch" if len(T_fenster) >= min_perioden else "nicht_periodisch")
    ergebnis.update(
        St=float(D / (U_inf * np.mean(T_fenster))),
        periode=float(np.mean(T_fenster)),
        periode_streuung=float(np.std(T_fenster) / np.mean(T_fenster)),
        amplitude=amplitude,
        mittelwert=float(mittel),
        t_einsatz=float(t[einsatz[0]]),
        t_fenster_start=float(t_a),
        t_fenster_ende=float(t_b),
        n_perioden=int(len(T_fenster)),
    )
    return ergebnis


def restschwankung(t, u, fenster=FENSTER_ENDE):
    """
    groesste abweichung des sondensignals vom endwert innerhalb der letzten <fenster>
    zeiteinheiten. klein heisst: der lauf ist in einen stationaeren zustand gelaufen und
    der status "keine_abloesung" ist physikalisch (kein anlauf, der noch nicht fertig ist)
    """
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    spaet = t >= t[-1] - fenster
    return float(np.max(np.abs(u[spaet] - u[-1])))


def sonden_lage(domain, cfg, abstand_in_D=2.0):
    """
    lage der nachlaufsonde (r = abstand_in_D * D, theta = 0) als zeilenindex i und gewicht w
    zwischen den xi-zeilen i und i+1. so sitzt die sonde auf jedem gitter an exakt derselben
    stelle; theta = 0 ist immer ein gitterpunkt.
    i + 2 muss noch im gebiet liegen, weil u_theta in zeile i + 1 zentral aus psi[i+2] folgt.
    """
    lage = math.log(abstand_in_D * cfg.D / cfg.R) / domain.dxi
    i = int(math.floor(lage))
    if i < 1 or i + 2 > domain.i_far:
        raise ValueError(f"sonde bei r = {abstand_in_D} D liegt nicht im inneren des gebiets")
    return i, lage - i
