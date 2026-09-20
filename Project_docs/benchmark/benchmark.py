"""
benchmark: serien von laeufen mit veraenderter orts- und zeitaufloesung sowie
gebietsgroesse. die serien und stufen stehen in presets.py, die beschreibung in README_Presets.md.

  python benchmark.py                          uebersicht ueber serien und stufen
  python benchmark.py gitter                   serie "gitter", stufe "schnell"
  python benchmark.py gitter --stufe voll      alle laeufe der serie
  python benchmark.py alle --stufe mittel      alle serien
  python benchmark.py zeit --parallel 4        laeufe parallel (zeitmessung dann nicht vergleichbar)
  python benchmark.py gitter --nur-auswerten   csv und diagramm aus vorhandenen ergebnissen
  python benchmark.py selbsttest               prueft sonde und auswertung an bekannten faellen

fertige laeufe werden gespeichert und beim naechsten aufruf uebersprungen (--neu erzwingt
eine neue rechnung). laeufe, die in mehreren serien vorkommen, werden nur einmal gerechnet.
"""
import argparse
import csv
import dataclasses
import json
import math
import os
import platform
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import scipy

#der loeser liegt eine ebene hoeher (Project_docs). der eintrag macht ihn importierbar,
#egal aus welchem ordner heraus benchmark/benchmark.py aufgerufen wird
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gitter import Config
from main import eine_Schleife
from presets import SERIEN, STUFEN, BENCHMARK_BASIS, BENCHMARK_T_END, laeufe_der_serie

SKRIPT_ORDNER = os.path.dirname(os.path.abspath(__file__))
STANDARD_AUSGABE = os.path.join(os.path.dirname(SKRIPT_ORDNER), "ergebnisse", "benchmark")
FORMAT_VERSION = 1

#zeitfenster am laufende, ueber das gemittelt wird, wenn es keine perioden gibt
#(stationaerer nachlauf bei kleinem Re). auch fuer die restschwankung
FENSTER_ENDE = 10.0


# ---------------------------------------------------------------------------
# messung waehrend des laufs
# ---------------------------------------------------------------------------

class Instabil(Exception):
    def __init__(self, t):
        super().__init__(f"rechnung bei t = {t:.3f} instabil")
        self.t = t


class Messung:
    """
    wird an eine_Schleife uebergeben und zeichnet waehrend des laufs auf:
    - quergeschwindigkeit u_theta an der sonde (r = 2D, theta = 0) nach jedem schritt
    - zeitintegral der wandwirbelstaerke (fuer den zeitlich gemittelten abloesewinkel)
    - zeitintegral der geschwindigkeit u_r auf der nachlaufachse (fuer die rueckstroemlaenge)
    psi und omega werden nur gelesen, der lauf selbst bleibt unveraendert.
    """

    def __init__(self, name, t_end, r_sonde_in_D=2.0, dt_kontrolle=0.05, dt_meldung=10.0, ausgabe=True):
        self.name = name
        self.t_end = t_end
        self.r_sonde_in_D = r_sonde_in_D
        self.dt_kontrolle = dt_kontrolle
        self.dt_meldung = dt_meldung
        self.ausgabe = ausgabe

    def start(self, domain, löser, psi, omega):
        self.uhr_start = time.perf_counter()
        self.domain = domain
        self.löser = löser
        cfg = domain.cfg

        #sonde zwischen zwei xi-zeilen linear interpolieren, damit sie auf jedem gitter
        #an exakt derselben stelle sitzt (theta = 0 ist immer ein gitterpunkt)
        lage = math.log(self.r_sonde_in_D * cfg.D / cfg.R) / domain.dxi
        self.i_sonde = int(math.floor(lage))
        self.w_sonde = lage - self.i_sonde
        if self.i_sonde < 1 or self.i_sonde + 2 > domain.i_far:
            raise ValueError(f"sonde bei r = {self.r_sonde_in_D} D liegt nicht im inneren des gebiets")

        self.t = [0.0]
        self.u = [self.sonde(psi)]

        #wandwirbelstaerke: laufendes integral (trapezregel) und kontrollpunkte, aus
        #denen spaeter der mittelwert ueber ein beliebiges zeitfenster folgt
        self.wand_alt = omega[domain.i_wall].copy()
        self.wand_integral = np.zeros(domain.n_theta)
        self.kontroll_t = [0.0]
        self.kontroll_integral = [self.wand_integral.copy()]

        #geschwindigkeit auf der nachlaufachse (theta = 0), genauso als laufendes integral.
        #daraus folgt spaeter die laenge des rueckstroemgebiets hinter dem zylinder
        self.achse_alt = self.achsengeschwindigkeit(psi)
        self.achse_integral = np.zeros(domain.n_xi)
        self.kontroll_achse = [self.achse_integral.copy()]

        self.naechste_kontrolle = self.dt_kontrolle
        self.naechste_meldung = self.dt_meldung

    def sonde(self, psi):
        #u_theta = -1/r * dpsi/dxi bei theta = 0, zentrale differenzen wie ops["D_xi"]
        d, i, w = self.domain, self.i_sonde, self.w_sonde
        u_i = -(psi[i + 1, 0] - psi[i - 1, 0]) / (2.0 * d.dxi * d.r[i])
        u_j = -(psi[i + 2, 0] - psi[i, 0]) / (2.0 * d.dxi * d.r[i + 1])
        return float((1.0 - w) * u_i + w * u_j)

    def achsengeschwindigkeit(self, psi):
        #u_r = (1/r) * dpsi/dtheta auf der achse theta = 0, zentrale differenzen wie ops["D_theta"].
        #dort zeigt die radiale richtung in x-richtung, u_r ist also die laengsgeschwindigkeit
        #im nachlauf: u_r < 0 heisst rueckstroemung
        d = self.domain
        return (psi[:, 1] - psi[:, -1]) / (2.0 * d.dtheta * d.r)

    def schritt(self, t, dt, psi, omega):
        u = self.sonde(psi)
        self.t.append(t)
        self.u.append(u)

        wand = omega[self.domain.i_wall]
        self.wand_integral += 0.5 * dt * (self.wand_alt + wand)
        self.wand_alt = wand.copy()

        achse = self.achsengeschwindigkeit(psi)
        self.achse_integral += 0.5 * dt * (self.achse_alt + achse)
        self.achse_alt = achse

        if not math.isfinite(u) or abs(u) > 1e3:
            raise Instabil(t)

        if t >= self.naechste_kontrolle:
            #gelegentlich das ganze feld pruefen, eine instabilitaet erreicht die sonde evtl. spaet
            if not np.all(np.isfinite(omega)) or np.max(np.abs(omega)) > 1e8:
                raise Instabil(t)
            self.kontroll_t.append(t)
            self.kontroll_integral.append(self.wand_integral.copy())
            self.kontroll_achse.append(self.achse_integral.copy())
            while self.naechste_kontrolle <= t:
                self.naechste_kontrolle += self.dt_kontrolle

        if self.ausgabe and t >= self.naechste_meldung:
            print(f"  [{self.name}] t = {t:6.1f} / {self.t_end:g}  "
                  f"({time.perf_counter() - self.uhr_start:5.0f} s)", flush=True)
            while self.naechste_meldung <= t:
                self.naechste_meldung += self.dt_meldung

    def abschluss(self):
        #letzten stand als kontrollpunkt sichern, damit das fenster bis zum ende reicht
        if self.kontroll_t[-1] < self.t[-1]:
            self.kontroll_t.append(self.t[-1])
            self.kontroll_integral.append(self.wand_integral.copy())
            self.kontroll_achse.append(self.achse_integral.copy())


# ---------------------------------------------------------------------------
# auswertung (reine funktionen, im selbsttest geprueft)
# ---------------------------------------------------------------------------

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


def mittel_aus_integral(kontroll_t, kontroll_integral, t_a, t_b):
    """zeitmittel ueber [t_a, t_b] aus dem laufenden integral, linear zwischen den kontrollpunkten"""
    kt = np.asarray(kontroll_t)
    KI = np.asarray(kontroll_integral)
    if t_a < kt[0] or t_b > kt[-1] or t_b <= t_a:
        return None

    def integral_bei(tz):
        j = min(max(int(np.searchsorted(kt, tz)) - 1, 0), len(kt) - 2)
        w = (tz - kt[j]) / (kt[j + 1] - kt[j])
        return (1.0 - w) * KI[j] + w * KI[j + 1]

    return (integral_bei(t_b) - integral_bei(t_a)) / (t_b - t_a)


def abloesewinkel(theta, omega_wand):
    """
    zeitlich gemittelter abloesewinkel in grad, gemessen vom vorderen staupunkt (theta = pi).
    oberseite: anliegende stroemung hat omega_wand < 0, hinter der abloesung (rueckstroemung)
    ist omega_wand > 0. unterseite spiegelbildlich. gibt (oben, unten) zurueck, nan ohne abloesung.
    """
    theta = np.asarray(theta)
    omega_wand = np.asarray(omega_wand)

    def suche(th, w):
        #erster wechsel von anliegend (w < 0) zu rueckstroemung (w >= 0). der letzte punkt ist
        #der hintere staupunkt (omega = 0 aus symmetrie) und zaehlt nicht, die suche endet davor
        for a in range(len(th) - 2):
            if w[a] < 0:
                if w[a + 1] >= 0:
                    th_null = th[a] + (th[a + 1] - th[a]) * w[a] / (w[a] - w[a + 1])
                    return abs(math.degrees(th_null) - 180.0)
        return math.nan

    oben = theta <= math.pi + 1e-12
    unten = theta >= math.pi - 1e-12
    winkel_oben = suche(theta[oben][::-1], omega_wand[oben][::-1])
    #unterseite von pi bis 2pi, der punkt theta = 2pi ist wieder theta = 0; vorzeichen gespiegelt
    th_u = np.append(theta[unten], 2.0 * math.pi)
    w_u = -np.append(omega_wand[unten], omega_wand[0])
    winkel_unten = suche(th_u, w_u)
    return winkel_oben, winkel_unten


def rueckstroemlaenge(xi, u_achse, R, D):
    """
    laenge des rueckstroemgebiets hinter dem zylinder, ab der zylinderrueckseite und in D.
    u_achse ist die laengsgeschwindigkeit auf der nachlaufachse theta = 0 ueber xi = ln(r/R).
    direkt an der wand ist u = 0 (haftbedingung), das gebiet beginnt also bei i = 1; gesucht
    ist der erste vorzeichenwechsel von rueckstroemung (u < 0) zu abstroemung (u >= 0).
    nan, wenn die stroemung schon am ersten punkt nach aussen zeigt (kein wirbelpaar).
    """
    u = np.asarray(u_achse)
    if len(u) < 3 or not (u[1] < 0.0):
        return math.nan
    for i in range(1, len(u) - 1):
        if u[i] < 0.0 <= u[i + 1]:
            #linear in xi interpolieren (dort liegen die punkte aequidistant), dann zurueck auf r
            w = u[i] / (u[i] - u[i + 1])
            xi_null = xi[i] + w * (xi[i + 1] - xi[i])
            return float((R * math.exp(xi_null) - R) / D)
    return math.nan


def restschwankung(t, u, fenster=10.0):
    """
    groesste abweichung des sondensignals vom endwert innerhalb der letzten <fenster>
    zeiteinheiten. klein heisst: der lauf ist in einen stationaeren zustand gelaufen und
    der status "keine_abloesung" ist physikalisch (kein anlauf, der noch nicht fertig ist)
    """
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    spaet = t >= t[-1] - fenster
    return float(np.max(np.abs(u[spaet] - u[-1])))


def richardson(h, f):
    """
    beobachtete konvergenzordnung p und extrapolierter wert fuer h -> 0 aus drei gittern mit
    konstantem verfeinerungsverhaeltnis (h grob -> fein). nan, wenn nicht monoton konvergent.
    """
    r = h[0] / h[1]
    if abs(h[1] / h[2] - r) > 1e-6 * r:
        return math.nan, math.nan
    e_grob, e_fein = f[0] - f[1], f[1] - f[2]
    if e_fein == 0 or e_grob * e_fein <= 0:
        return math.nan, math.nan
    p = math.log(e_grob / e_fein) / math.log(r)
    return p, f[2] + (f[2] - f[1]) / (r**p - 1.0)


# ---------------------------------------------------------------------------
# ein lauf
# ---------------------------------------------------------------------------

def lauf_name(cfg, t_end):
    return (f"nxi{cfg.n_xi:03d}_nth{cfg.n_theta:03d}_cfl{cfg.cfl_target:g}"
            f"_rmax{cfg.r_max:.4g}_Re{cfg.Re:g}_t{t_end:g}")


def messe_loesezeit(löser, domain, mindestdauer=0.3):
    omega = np.random.default_rng(0).standard_normal((domain.n_xi, domain.n_theta))
    löser.löse(omega)
    anzahl, uhr = 0, time.perf_counter()
    while anzahl < 10 or time.perf_counter() - uhr < mindestdauer:
        löser.löse(omega)
        anzahl += 1
    return (time.perf_counter() - uhr) / anzahl


def zerlegung_groesse(löser):
    #eintraege in den LR-faktoren. der FFT-loeser haelt statt einer grossen zerlegung
    #viele kleine (eine je fourier-mode)
    zerlegungen = [löser.lu] if hasattr(löser, "lu") else löser.löser_pro_mode
    return int(sum(z.L.nnz + z.U.nnz for z in zerlegungen))


def _json_wert(x):
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def fuehre_lauf_aus(name, cfg_werte, t_end, ordner, parallel, ausgabe=True):
    """rechnet einen lauf, wertet ihn aus und speichert <name>.json und <name>.npz"""
    cfg = Config(**cfg_werte)
    messung = Messung(name, t_end, ausgabe=ausgabe)
    instabil_bei = math.nan

    uhr0 = time.perf_counter()
    try:
        with np.errstate(over="ignore", invalid="ignore"):
            eine_Schleife(cfg, t_end=t_end, verbose=False, t_speicher=math.inf, messung=messung)
    except Instabil as fehler:
        instabil_bei = fehler.t
    uhr1 = time.perf_counter()
    messung.abschluss()

    t = np.array(messung.t)
    u = np.array(messung.u)
    schritte = len(t) - 1
    domain = messung.domain

    if math.isnan(instabil_bei):
        auswertung = werte_signal_aus(t, u, cfg.D, cfg.U_inf)
    else:
        auswertung = leere_auswertung("instabil")

    #zeitmittel von wandwirbelstaerke und achsengeschwindigkeit. bei periodischer abloesung
    #ueber das periodenfenster, sonst (stationaerer nachlauf) ueber die letzten FENSTER_ENDE
    #zeiteinheiten - so sind abloesewinkel und rueckstroemlaenge auch unterhalb der
    #kritischen reynolds-zahl definiert
    wand_mittel = achse_mittel = None
    winkel_oben = winkel_unten = laenge = math.nan
    fenster_art = "keins"
    if math.isnan(instabil_bei):
        if auswertung["status"] == "periodisch":
            t_a, t_b = auswertung["t_fenster_start"], auswertung["t_fenster_ende"]
            fenster_art = "perioden"
        else:
            t_b = float(t[-1])
            t_a = max(t_b - FENSTER_ENDE, 0.5 * t_b)
            fenster_art = "ende"
        wand_mittel = mittel_aus_integral(messung.kontroll_t, messung.kontroll_integral, t_a, t_b)
        achse_mittel = mittel_aus_integral(messung.kontroll_t, messung.kontroll_achse, t_a, t_b)
        if wand_mittel is not None:
            winkel_oben, winkel_unten = abloesewinkel(domain.theta, wand_mittel)
        if achse_mittel is not None:
            laenge = rueckstroemlaenge(domain.xi, achse_mittel, cfg.R, cfg.D)

    t_aufbau = messung.uhr_start - uhr0
    t_schleife = uhr1 - messung.uhr_start
    ms_pro_schritt = 1e3 * t_schleife / max(schritte, 1)
    t_loese_ms = 1e3 * messe_loesezeit(messung.löser, domain)
    dt_mittel = t[-1] / max(schritte, 1)
    nnz_lu = zerlegung_groesse(messung.löser)

    ergebnis = dict(
        format_version=FORMAT_VERSION,
        name=name,
        config=dataclasses.asdict(cfg),
        t_end=t_end,
        **auswertung,
        instabil_bei_t=instabil_bei,
        abloesewinkel_oben=winkel_oben,
        abloesewinkel_unten=winkel_unten,
        abloesewinkel=0.5 * (winkel_oben + winkel_unten),
        rueckstroemlaenge=laenge,
        fenster_art=fenster_art,
        #restschwankung relativ zu U_inf: < 1e-4 heisst stationaer eingelaufen
        restschwankung=restschwankung(t, u, FENSTER_ENDE) / cfg.U_inf if schritte > 0 else math.nan,
        schritte=schritte,
        t_erreicht=float(t[-1]),
        dt_mittel=dt_mittel,
        dt_min=float(np.min(np.diff(t))) if schritte > 0 else math.nan,
        t_aufbau_s=t_aufbau,
        t_schleife_s=t_schleife,
        ms_pro_schritt=ms_pro_schritt,
        t_loese_ms=t_loese_ms,
        poisson_anteil=4.0 * t_loese_ms / ms_pro_schritt,
        rechenzeit_pro_periode_s=ms_pro_schritt * 1e-3 * auswertung["periode"] / dt_mittel,
        unbekannte=cfg.n_xi * cfg.n_theta,
        poisson_loeser=type(messung.löser).__name__,
        nnz_LU=nnz_lu,
        speicher_LU_MB=nnz_lu * 12 / 1e6,     #8 byte wert + 4 byte index je eintrag
        parallel=parallel,
        zeitmessung_vergleichbar=(parallel == 1),
        rechner=dict(cpu=platform.processor(), python=platform.python_version(),
                     numpy=np.__version__, scipy=scipy.__version__),
        datum=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    ergebnis = {k: _json_wert(v) for k, v in ergebnis.items()}

    os.makedirs(ordner, exist_ok=True)
    np.savez(os.path.join(ordner, name + ".npz"), t=t, u=u, theta=domain.theta, xi=domain.xi,
             wand_mittel=wand_mittel if wand_mittel is not None else np.array([]),
             achse_mittel=achse_mittel if achse_mittel is not None else np.array([]))
    #erst in eine temporaere datei, damit ein abbruch keine halbe json hinterlaesst
    pfad = os.path.join(ordner, name + ".json")
    with open(pfad + ".tmp", "w", encoding="utf-8") as datei:
        json.dump(ergebnis, datei, indent=2, ensure_ascii=False)
    os.replace(pfad + ".tmp", pfad)
    return ergebnis


def ist_fertig(ordner, name, cfg, t_end):
    pfad = os.path.join(ordner, name + ".json")
    if not os.path.exists(pfad):
        return False
    with open(pfad, encoding="utf-8") as datei:
        alt = json.load(datei)
    return (alt.get("format_version") == FORMAT_VERSION and alt.get("t_end") == t_end
            and alt.get("config") == dataclasses.asdict(cfg))


def kurzbericht(e):
    St = f"{e['St']:.4f}" if e["St"] is not None else "  -   "
    winkel = f"{e['abloesewinkel']:.1f}°" if e["abloesewinkel"] is not None else "  -  "
    laenge = f"{e['rueckstroemlaenge']:.3f}" if e.get("rueckstroemlaenge") is not None else "  -  "
    return (f"{e['name']}: {e['status']}, St = {St}, abloesewinkel = {winkel}, L/D = {laenge}, "
            f"{e['schritte']} schritte, {e['ms_pro_schritt']:.2f} ms/schritt "
            f"({poisson_kurz(e.get('poisson_loeser'))}), gesamt {e['t_aufbau_s'] + e['t_schleife_s']:.0f} s")


def poisson_kurz(name):
    #kurzname des poisson-loesers fuer tabelle und meldungen
    return {"PoissonSolver": "LR", "PoissonSolverFFT": "FFT"}.get(name, name or "?")


# ---------------------------------------------------------------------------
# serien: ausfuehren, tabelle, diagramm
# ---------------------------------------------------------------------------

CSV_SPALTEN = ["name", "x", "stufe", "n_xi", "n_theta", "cfl_target", "r_max", "Re", "status", "St",
               "periode_streuung", "amplitude", "rueckstroemlaenge", "abloesewinkel",
               "abloesewinkel_oben", "abloesewinkel_unten", "restschwankung", "fenster_art",
               "t_einsatz", "n_perioden", "instabil_bei_t", "schritte",
               "dt_mittel", "t_aufbau_s", "t_schleife_s", "ms_pro_schritt", "t_loese_ms",
               "poisson_anteil", "rechenzeit_pro_periode_s", "unbekannte", "poisson_loeser",
               "nnz_LU", "zeitmessung_vergleichbar"]


def lade_serie(serie, t_end, laufordner):
    zeilen = []
    for stufe, cfg in laeufe_der_serie(serie, "voll"):
        name = lauf_name(cfg, t_end)
        pfad = os.path.join(laufordner, name + ".json")
        if not ist_fertig(laufordner, name, cfg, t_end):
            continue
        with open(pfad, encoding="utf-8") as datei:
            e = json.load(datei)
        e.update(x=SERIEN[serie]["x"](cfg), stufe=stufe, **dataclasses.asdict(cfg))
        zeilen.append(e)
    #nach dem serienparameter sortieren, nicht in der reihenfolge aus presets.py:
    #sonst laufen die linien im diagramm hin und her
    zeilen.sort(key=lambda e: e["x"])
    return zeilen


def werte_serie_aus(serie, t_end, ordner):
    laufordner = os.path.join(ordner, "laeufe")
    zeilen = lade_serie(serie, t_end, laufordner)
    if not zeilen:
        print(f"serie {serie}: noch keine ergebnisse")
        return

    pfad_csv = os.path.join(ordner, f"{serie}.csv")
    with open(pfad_csv, "w", newline="", encoding="utf-8") as datei:
        schreiber = csv.DictWriter(datei, fieldnames=CSV_SPALTEN, extrasaction="ignore")
        schreiber.writeheader()
        for z in zeilen:
            schreiber.writerow({k: ("" if z.get(k) is None else z.get(k)) for k in CSV_SPALTEN})

    print(f"\nserie {serie} ({SERIEN[serie]['beschreibung']})")
    print(f"  {SERIEN[serie]['x_name']:>10s}  {'status':16s} {'St':>7s} {'amplitude':>9s} "
          f"{'L/D':>6s} {'winkel':>7s} {'einsatz':>7s} {'perioden':>8s} {'dt_mittel':>9s} "
          f"{'ms/schr.':>8s} {'gesamt':>8s} {'poisson':>7s}")
    for z in zeilen:
        def f(wert, format_):
            return format_.format(wert) if wert is not None else "-"
        print(f"  {z['x']:>10.4g}  {z['status']:16s} {f(z['St'], '{:.4f}'):>7s} "
              f"{f(z['amplitude'], '{:.4f}'):>9s} {f(z.get('rueckstroemlaenge'), '{:.3f}'):>6s} "
              f"{f(z['abloesewinkel'], '{:.1f}'):>7s} "
              f"{f(z['t_einsatz'], '{:.1f}'):>7s} {z['n_perioden']:>8d} {z['dt_mittel']:>9.2e} "
              f"{z['ms_pro_schritt']:>8.2f} {z['t_aufbau_s'] + z['t_schleife_s']:>7.0f}s "
              f"{poisson_kurz(z.get('poisson_loeser')):>7s}")

    zusatz = ""
    if serie == "gitter":
        fam = {z["n_theta"]: z for z in zeilen if z["status"] == "periodisch"}
        if all(n in fam for n in (80, 160, 320)):
            h = [2 * math.pi / n for n in (80, 160, 320)]
            teile = []
            for groesse in ("St", "amplitude", "abloesewinkel"):
                werte = [fam[n][groesse] for n in (80, 160, 320)]
                if None in werte:
                    continue
                p, extra = richardson(h, werte)
                if math.isnan(p):
                    text = f"{groesse}: nicht monoton konvergent"
                else:
                    text = f"{groesse}: p = {p:.2f}, h→0: {extra:.4g}"
                print(f"  richardson {text}")
                teile.append(text)
            zusatz = "Richardson (80/160/320):  " + ";  ".join(teile)

    pfad_png = os.path.join(ordner, f"{serie}.png")
    zeichne_serie(pfad_png, serie, zeilen, zusatz)
    print(f"  gespeichert: {pfad_csv}, {pfad_png}")


#die ersten drei felder des diagramms. eine serie kann in presets.py mit "panels"
#andere groessen waehlen (die serie "reynolds" z.b. die rueckstroemlaenge)
STANDARD_PANELS = ("St", "amplitude", "abloesewinkel")
PANEL_TITEL = {
    "St": "Strouhal-Zahl",
    "amplitude": r"Amplitude $u_	heta$ bei $r = 2D$",
    "abloesewinkel": "Ablösewinkel [°] (vom Staupunkt)",
    "rueckstroemlaenge": "Rückströmlänge $L/D$ (ab Zylinderrückseite)",
}


def zeichne_referenz(ax, referenz, x):
    """literaturvergleich einzeichnen: kurve (aufrufbar) oder einzelne punkte (dict)"""
    if callable(referenz):
        xr = np.linspace(x.min(), x.max(), 200)
        yr = np.array([referenz(v) for v in xr])
        if np.any(np.isfinite(yr)):
            ax.plot(xr, yr, "k--", lw=1.2, zorder=0, label="Literatur")
    else:
        ax.plot(list(referenz), list(referenz.values()), "k*", ms=11, zorder=0, label="Literatur")


def zeichne_serie(pfad, serie, zeilen, zusatz):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    x = np.array([z["x"] for z in zeilen], dtype=float)
    status = [z["status"] for z in zeilen]

    def werte(schluessel):
        #.get, weil aeltere ergebnisdateien neuere groessen noch nicht enthalten
        return np.array([np.nan if z.get(schluessel) is None else z[schluessel] for z in zeilen],
                        dtype=float)

    referenzen = SERIEN[serie].get("referenz", {})
    fig, achsen = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, schluessel in zip(achsen.flat[:3], SERIEN[serie].get("panels", STANDARD_PANELS)):
        y = werte(schluessel)
        periodisch = np.array([s == "periodisch" for s in status])
        ax.plot(x[periodisch], y[periodisch], "o-", color="C0", label="periodisch")
        andere = ~periodisch & np.isfinite(y)
        if andere.any():
            ax.plot(x[andere], y[andere], "o", mfc="none", color="C1", label="ohne Ablösung")
        if schluessel in referenzen:
            zeichne_referenz(ax, referenzen[schluessel], x)
            ax.legend(fontsize=8)
        ax.set_title(PANEL_TITEL.get(schluessel, schluessel))

    def y_achse(achse, y):
        #logarithmisch nur, wenn die werte ueber mehr als eine groessenordnung reichen
        y = y[np.isfinite(y)]
        if len(y) and y.min() > 0 and y.max() / y.min() > 10:
            achse.set_yscale("log")
        achse.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        achse.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))

    ax = achsen.flat[3]
    gesamt = werte("t_aufbau_s") + werte("t_schleife_s")
    ax.plot(x, gesamt, "s-", color="C2", label="Gesamt [s]")
    y_achse(ax, gesamt)
    ax.set_ylabel("Rechenzeit gesamt [s]")
    ax2 = ax.twinx()
    ax2.plot(x, werte("ms_pro_schritt"), "^--", color="C3", label="ms / Schritt")
    y_achse(ax2, werte("ms_pro_schritt"))
    ax2.set_ylabel("ms pro Zeitschritt")
    titel = "Rechenaufwand"
    if not all(z["zeitmessung_vergleichbar"] for z in zeilen):
        titel += " (teils parallel gemessen, nicht vergleichbar)"
    ax.set_title(titel)
    linien = ax.get_legend_handles_labels()
    linien2 = ax2.get_legend_handles_labels()
    ax.legend(linien[0] + linien2[0], linien[1] + linien2[1], fontsize=8)

    for ax in achsen.flat:
        ax.set_xscale("log", base=2)
        ax.set_xticks(np.unique(x))
        ax.set_xticklabels([f"{v:.4g}" for v in np.unique(x)])
        ax.minorticks_off()
        ax.set_xlim(x.min() / 1.15, x.max() * 1.15)    #rand, damit markierungen am ende sichtbar bleiben
        ax.set_xlabel(SERIEN[serie]["x_name"])
        ax.grid(alpha=0.3)
    for z in zeilen:
        if z["status"] in ("instabil", "keine_abloesung", "nicht_periodisch"):
            #"keine_abloesung" ist bei kleinem Re das physikalisch richtige ergebnis und kein
            #fehlschlag - dann grau und als "stationaer" beschriften, sofern eingelaufen
            stationaer = (z["status"] != "instabil"
                          and (z.get("restschwankung") or 1.0) < 1e-3)
            farbe = "C7" if stationaer else "C3"
            text = "stationär" if stationaer else z["status"].replace("_", " ")
            achsen.flat[0].axvline(z["x"], color=farbe, ls=":", lw=1)
            achsen.flat[0].annotate(text, (z["x"], 0.02), xycoords=("data", "axes fraction"),
                                    rotation=90, fontsize=8, color=farbe, ha="right")
    if "St" not in SERIEN[serie].get("referenz", {}):
        achsen.flat[0].legend(fontsize=8)

    fig.suptitle(f"Benchmark „{serie}“: {SERIEN[serie]['beschreibung']}" + (f"\n{zusatz}" if zusatz else ""),
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(pfad, dpi=120)
    plt.close(fig)


def uebersicht():
    print(__doc__)
    for serie, s in SERIEN.items():
        print(f"serie '{serie}': {s['beschreibung']}")
        for stufe, cfg in laeufe_der_serie(serie, "voll"):
            print(f"  {stufe:8s} n_xi = {cfg.n_xi:3d}, n_theta = {cfg.n_theta:3d}, "
                  f"cfl_target = {cfg.cfl_target:<5g}, r_max = {cfg.r_max:.4g} ({cfg.r_max / cfg.D:.1f} D)")
        print()


def main():
    parser = argparse.ArgumentParser(description="benchmark-serien, siehe README_Presets.md")
    parser.add_argument("serie", nargs="?", choices=list(SERIEN) + ["alle", "selbsttest"])
    parser.add_argument("--stufe", choices=STUFEN, default="schnell")
    parser.add_argument("--parallel", type=int, default=1, help="anzahl gleichzeitiger laeufe")
    parser.add_argument("--neu", action="store_true", help="vorhandene ergebnisse neu rechnen")
    parser.add_argument("--nur-auswerten", action="store_true", help="nichts rechnen, nur csv/diagramm")
    parser.add_argument("--t-end", type=float, default=BENCHMARK_T_END)
    parser.add_argument("--ausgabe", default=STANDARD_AUSGABE, help="ergebnisordner")
    args = parser.parse_args()
    #zeilenweise ausgeben, sonst mischen sich im parallelbetrieb die ausgaben der prozesse falsch.
    #utf-8, damit grad- und pfeilzeichen auch beim umleiten in eine datei funktionieren
    #(sonst faellt windows auf cp1252 zurueck und die ausgabe bricht mit UnicodeEncodeError ab)
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

    if args.serie is None:
        uebersicht()
        return
    if args.serie == "selbsttest":
        sys.exit(0 if selbsttest() else 1)

    serien = list(SERIEN) if args.serie == "alle" else [args.serie]
    laufordner = os.path.join(args.ausgabe, "laeufe")

    if not args.nur_auswerten:
        #doppelte laeufe (gleiche config in mehreren serien) nur einmal
        aufgaben = {}
        for serie in serien:
            for stufe, cfg in laeufe_der_serie(serie, args.stufe):
                aufgaben.setdefault(lauf_name(cfg, args.t_end), cfg)

        offen = {}
        for name, cfg in aufgaben.items():
            if not args.neu and ist_fertig(laufordner, name, cfg, args.t_end):
                print(f"vorhanden, uebersprungen: {name}")
            else:
                offen[name] = cfg
        print(f"{len(offen)} von {len(aufgaben)} laeufen zu rechnen (stufe {args.stufe}, t_end = {args.t_end:g})")

        if args.parallel > 1 and offen:
            print(f"achtung: {args.parallel} laeufe gleichzeitig - die zeitmessungen sind dann nicht "
                  f"vergleichbar und werden entsprechend markiert")
            with ProcessPoolExecutor(max_workers=args.parallel) as pool:
                zukuenfte = {pool.submit(fuehre_lauf_aus, name, dataclasses.asdict(cfg), args.t_end,
                                         laufordner, args.parallel): name for name, cfg in offen.items()}
                for zukunft in as_completed(zukuenfte):
                    try:
                        print("fertig: " + kurzbericht(zukunft.result()), flush=True)
                    except Exception:
                        print(f"FEHLER in {zukuenfte[zukunft]}:\n{traceback.format_exc()}", flush=True)
        else:
            for name, cfg in offen.items():
                print(f"starte {name}", flush=True)
                try:
                    print("fertig: " + kurzbericht(
                        fuehre_lauf_aus(name, dataclasses.asdict(cfg), args.t_end, laufordner, 1)), flush=True)
                except Exception:
                    print(f"FEHLER in {name}:\n{traceback.format_exc()}", flush=True)

    for serie in serien:
        werte_serie_aus(serie, args.t_end, args.ausgabe)


# ---------------------------------------------------------------------------
# selbsttest
# ---------------------------------------------------------------------------

def selbsttest():
    from gitter import Domain
    from loeser import build_alle_operatoren, geschwindigkeit

    alles_ok = True

    def pruefe(bedingung, text):
        nonlocal alles_ok
        alles_ok &= bool(bedingung)
        print(f"  [{'ok' if bedingung else 'FEHLER'}] {text}")

    print("sonde:")
    #1) auf einer gitterzeile (w = 0) muss die sonde exakt geschwindigkeit() entsprechen
    cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=81, n_theta=160, dt=0.05)
    d = Domain(cfg)
    psi = np.random.default_rng(1).standard_normal((d.n_xi, d.n_theta))
    m = Messung("test", 1.0, ausgabe=False)
    m.start(d, None, psi, np.zeros_like(psi))
    _, u_theta = geschwindigkeit(psi, d, build_alle_operatoren(d))
    m.i_sonde, m.w_sonde = 30, 0.0
    abw = abs(m.sonde(psi) - u_theta[30, 0])
    pruefe(abw < 1e-12 * max(1.0, abs(u_theta[30, 0])), f"gitterzeile gegen geschwindigkeit(): abweichung {abw:.1e}")

    #2) psi = r^3 -> u_theta = -dpsi/dr = -3 r^2 = -12 bei r = 2, fehler 2. ordnung
    fehler = []
    for n in (41, 81, 161):
        d = Domain(Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=n, n_theta=2 * (n - 1), dt=0.05))
        psi = np.tile(d.r[:, None] ** 3, (1, d.n_theta))
        m = Messung("test", 1.0, ausgabe=False)
        m.start(d, None, psi, np.zeros_like(psi))
        fehler.append(abs(m.sonde(psi) + 12.0))
    ordnung = [math.log2(fehler[0] / fehler[1]), math.log2(fehler[1] / fehler[2])]
    pruefe(fehler[1] < 0.01 * 12 and 1.8 < min(ordnung) and max(ordnung) < 2.2,
           f"psi = r^3 bei r = 2D: fehler {fehler[0]:.1e} / {fehler[1]:.1e} / {fehler[2]:.1e}, "
           f"ordnung {ordnung[0]:.2f} / {ordnung[1]:.2f} (soll 2)")

    print("signalauswertung:")
    t = np.arange(0.0, 150.0, 0.005)
    St_soll, A_soll, versatz = 0.1644, 0.55, 0.01
    huelle = A_soll / (1.0 + np.exp(-(t - 40.0) / 3.0))
    u = versatz + huelle * np.sin(2 * math.pi * St_soll * t)
    e = werte_signal_aus(t, u, 1.0, 1.0)
    pruefe(e["status"] == "periodisch" and abs(e["St"] - St_soll) < 1e-5,
           f"anlaufendes sinussignal: St = {e['St']:.6f} (soll {St_soll})")
    pruefe(abs(e["amplitude"] / A_soll - 1) < 0.01, f"amplitude = {e['amplitude']:.4f} (soll {A_soll})")
    pruefe(abs(e["mittelwert"] - versatz) < 1e-3, f"mittelwert = {e['mittelwert']:.5f} (soll {versatz})")
    pruefe(40.0 <= e["t_einsatz"] <= 40.0 + 1.0 / St_soll,
           f"einsatz bei t = {e['t_einsatz']:.2f} (halbe amplitude ab t = 40)")
    pruefe(e["t_fenster_start"] > 40.0 + 3.0 * math.log(49.0) - 1.0 / St_soll,
           f"fenster beginnt bei t = {e['t_fenster_start']:.1f} (huelle erst ab t = 51.7 innerhalb 2 %)")

    e = werte_signal_aus(t, 0.3 * np.exp(-t / 10.0) * np.sin(2 * math.pi * St_soll * t), 1.0, 1.0)
    pruefe(e["status"] == "keine_abloesung", f"abklingendes signal: status {e['status']}")
    kurz = t < 15.0    #volle amplitude, aber nur 2.5 perioden
    e = werte_signal_aus(t[kurz], A_soll * np.sin(2 * math.pi * St_soll * t[kurz]), 1.0, 1.0)
    pruefe(e["status"] == "nicht_periodisch", f"zu kurzes signal (2.5 perioden): status {e['status']}")

    print("wandwirbelstaerke und abloesewinkel:")
    #laufendes integral von f(t) = a + b cos(wt), mittel ueber ganze perioden = a
    a, b, w = 2.0, 5.0, 2 * math.pi * St_soll
    kt = np.arange(0.0, 100.0, 0.05)
    KI = (a * kt + b / w * np.sin(w * kt))[:, None]
    t_a = 20.0
    mittel = mittel_aus_integral(kt, KI, t_a, t_a + 10.0 / St_soll)[0]
    pruefe(abs(mittel - a) < 1e-3 * b, f"zeitmittel ueber 10 perioden: {mittel:.6f} (soll {a})")

    d = Domain(BENCHMARK_BASIS)
    th_s = math.radians(180.0 - 117.0)
    omega_wand = -np.sin(d.theta) * (math.cos(th_s) - np.cos(d.theta))
    oben, unten = abloesewinkel(d.theta, omega_wand)
    pruefe(abs(oben - 117.0) < 0.1 and abs(unten - 117.0) < 0.1,
           f"vorgegebene abloesung bei 117°: oben {oben:.3f}°, unten {unten:.3f}°")
    oben, unten = abloesewinkel(d.theta, -np.sin(d.theta))
    pruefe(math.isnan(oben) and math.isnan(unten), "anliegende stroemung ohne vorzeichenwechsel: keine abloesung")

    print("rueckstroemlaenge:")
    d = Domain(BENCHMARK_BASIS)
    #u linear in xi mit nulldurchgang bei r = R + 1.0 D: die lineare interpolation in xi
    #muss diesen punkt exakt treffen
    xi_null = math.log((BENCHMARK_BASIS.R + BENCHMARK_BASIS.D) / BENCHMARK_BASIS.R)
    u_achse = d.xi - xi_null
    u_achse[0] = 0.0                      #an der wand gilt die haftbedingung
    L = rueckstroemlaenge(d.xi, u_achse, BENCHMARK_BASIS.R, BENCHMARK_BASIS.D)
    pruefe(abs(L - 1.0) < 1e-12, f"vorgegebener nulldurchgang bei L/D = 1: gemessen {L:.12f}")

    ohne = np.ones(d.n_xi)
    ohne[0] = 0.0
    pruefe(math.isnan(rueckstroemlaenge(d.xi, ohne, BENCHMARK_BASIS.R, BENCHMARK_BASIS.D)),
           "stroemung ueberall nach aussen: keine rueckstroemlaenge")

    #eine negative zone weiter aussen, die den zylinder nicht beruehrt, zaehlt nicht
    getrennt = np.ones(d.n_xi)
    getrennt[0] = 0.0
    getrennt[20:30] = -1.0
    pruefe(math.isnan(rueckstroemlaenge(d.xi, getrennt, BENCHMARK_BASIS.R, BENCHMARK_BASIS.D)),
           "abgeloeste negative zone ohne wandkontakt: keine rueckstroemlaenge")

    print("restschwankung:")
    t_s = np.arange(0.0, 150.0, 0.005)
    rest_ab = restschwankung(t_s, 0.3 * np.exp(-t_s / 10.0) * np.sin(2 * math.pi * 0.1644 * t_s))
    rest_per = restschwankung(t_s, 0.55 * np.sin(2 * math.pi * 0.1644 * t_s))
    pruefe(rest_ab < 1e-4, f"abgeklungenes signal: restschwankung {rest_ab:.2e} (soll klein)")
    pruefe(rest_per > 0.9 * 0.55, f"periodisches signal: restschwankung {rest_per:.4f} (soll ~ amplitude)")

    print("richardson:")
    p, extra = richardson([1.0, 0.5, 0.25], [0.17 - 0.02 * hh**2 for hh in (1.0, 0.5, 0.25)])
    pruefe(abs(p - 2.0) < 1e-9 and abs(extra - 0.17) < 1e-12, f"f = 0.17 - 0.02 h^2: p = {p:.6f}, h->0: {extra:.6f}")

    print("selbsttest " + ("bestanden" if alles_ok else "FEHLGESCHLAGEN"))
    return alles_ok


if __name__ == "__main__":
    main()
