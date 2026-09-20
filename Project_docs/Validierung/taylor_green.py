"""
validierung am taylor-green-wirbel - dem pflichtbeispiel der aufgabenstellung und dem
einzigen fall, in dem advektion, poisson-loesung und rk4 GEMEINSAM gegen eine exakte
loesung der nichtlinearen gleichungen gestellt werden koennen.

  python Validierung/taylor_green.py

idee: die log-polaren gleichungen des loesers enthalten den kartesischen fall als exakten
spezialfall. mit r = 1 (metrischer vorfaktor 1) und xi = x, theta = y wird

    laplace       (1/r^2)(d2/dxi2 + d2/dtheta2)           ->  d2/dx2 + d2/dy2
    geschwind.    u_r = (1/r) dpsi/dtheta                 ->  u =  dpsi/dy
                  u_theta = -(1/r) dpsi/dxi               ->  v = -dpsi/dx
    advektion     (1/r)(u_r dom/dxi + u_theta dom/dtheta) ->  u dom/dx + v dom/dy
    zeitschritt   aus r*dxi und r*dtheta                  ->  aus dx und dy

zeile fuer zeile die kartesische formulierung. deshalb laufen geschwindigkeit(),
berechne_rhs(), cfl_zeitschritt() und rk4() hier UNVERAENDERT - es ist derselbe code, der
auch die zylinderumstroemung rechnet, und kein nachbau. ersetzt werden nur die beiden
teile, die es periodisch nicht gibt: der poisson-loeser (keine dirichlet-raender) und die
randbedingungen (rk4 bekommt dafuer eine leere funktion uebergeben).

exakte loesung auf [0, 2pi)^2, in beiden richtungen periodisch:

    psi = sin(x) sin(y) exp(-2 nu t),    omega = 2 psi

sie loest die nichtlinearen gleichungen, weil omega ein festes vielfaches von psi ist
(beltrami-stroemung): dann steht der gradient von psi senkrecht auf der geschwindigkeit
und u.grad(omega) verschwindet identisch.

WAS DIESER TEST NICHT PRUEFT: mit r = 1 und periodischen raendern sind die metrik 1/r^2,
die thom-formel an der wand, die fernfeld-randbedingung und die einseitigen randzeilen von
build_D_xi/build_D2_xi abgeschaltet. und weil die nichtlinearitaet im taylor-green-fall
identisch verschwindet, sagt er ueber die GENAUIGKEIT des advektionsterms nichts aus -
dafuer steht unten ein eigener test mit einer frei gewaehlten loesung, bei der omega kein
vielfaches von psi ist. die zylinderspezifischen teile decken selbsttest_loeser.py
(abnahmetest gegen die potentialstroemung) und der benchmark (serie reynolds,
literaturvergleich) ab.
"""

import numpy as np
import scipy.sparse as sp

from pruefung import Pruefung      #setzt den suchpfad auf Project_docs

from gitter import Config
from operatoren import build_D_theta, build_D2_theta
from loeser import berechne_rhs, cfl_zeitschritt, geschwindigkeit, rk4


# ---------------------------------------------------------------------------
# das periodische quadrat in der schnittstelle eines Domain
# ---------------------------------------------------------------------------

class KartesischesGebiet:
    """
    quadratisches, in beide richtungen periodisches gebiet [0, 2pi)^2. bietet genau die
    attribute, die der loeser von einem Domain benutzt - mit r = 1 und vorfaktor = 1.

    i_wall und i_far fehlen bewusst: apply_bc, build_rhs und build_poisson_matrix duerfen
    hier nicht verwendet werden (es gibt keine wand und keinen fernrand). wer es doch
    versucht, bekommt einen AttributeError statt still falscher ergebnisse.
    """

    def __init__(self, n, cfg):
        self.cfg = cfg
        self.n_xi = self.n_theta = n            #xi = x, theta = y
        self.dxi = self.dtheta = 2.0 * np.pi / n
        self.xi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        self.theta = self.xi.copy()
        #r = 1: damit fallen alle metrischen faktoren des loesers auf 1
        self.r = np.ones(n)
        self.vorfaktor = np.ones(n)

    def flatten(self, feld2d):
        return feld2d.reshape(-1)

    def unflatten(self, feld1d):
        return feld1d.reshape(self.n_xi, self.n_theta)

    def gitter(self):
        return np.meshgrid(self.xi, self.theta, indexing="ij")     #X, Y


def operatoren_periodisch(gebiet):
    """
    operator-dict wie build_alle_operatoren, aber in BEIDEN richtungen periodisch.
    build_D_theta und build_D2_theta lesen nur (n, h) aus dem gebiet und sind periodisch -
    auf dem quadratischen gitter ist das in x- und y-richtung dieselbe matrix.
    """
    D = build_D_theta(gebiet)
    D2 = build_D2_theta(gebiet)
    I = sp.identity(gebiet.n_xi, format="csr")
    return {
        "D_xi": sp.kron(D, I, format="csr"),        #d/dx
        "D_theta": sp.kron(I, D, format="csr"),     #d/dy
        #der metrische vorfaktor ist 1, der laplace also die reine summe
        "L": sp.kron(D2, I, format="csr") + sp.kron(I, D2, format="csr"),
    }


class PoissonPeriodisch:
    """
    loest nabla^2 psi = -omega auf dem periodischen quadrat, als ersatz fuer PoissonSolver.
    build_poisson_matrix passt hier nicht: sie ersetzt erste und letzte xi-zeile durch
    dirichlet-zeilen (wand, fernfeld), die es periodisch nicht gibt. ohne sie ist die
    matrix singulaer - psi ist nur bis auf eine konstante bestimmt, was die
    geschwindigkeiten nicht beruehrt.

    geloest wird spektral: D2 ist in beiden richtungen zyklisch und wird von der DFT
    diagonalisiert. invertiert werden die eigenwerte der DISKRETEN matrix (aus ihrer ersten
    zeile, genau wie in PoissonSolverFFT) und nicht die des analytischen operators -k^2 -
    sonst waere die poisson-loesung spektral genau und der rest der rechnung 2. ordnung.
    """

    def __init__(self, gebiet):
        self.gebiet = gebiet
        erste_zeile = np.asarray(build_D2_theta(gebiet)[0, :].todense()).ravel()
        lam = np.fft.fft(erste_zeile).real          #eigenwerte der zyklischen matrix
        self.eigenwerte = lam[:, None] + lam[None, :]
        #die nullmode (konstanter anteil) ist nicht bestimmt; die 1.0 verhindert die
        #division durch null, ihr koeffizient wird in löse() auf 0 gesetzt
        self.eigenwerte[0, 0] = 1.0

    def löse(self, omega):
        psi_dach = -np.fft.fft2(omega) / self.eigenwerte
        psi_dach[0, 0] = 0.0
        return np.real(np.fft.ifft2(psi_dach))


def ohne_randbedingungen(psi, omega, gebiet):
    """ein periodisches gebiet hat keinen rand - psi und omega bleiben unveraendert"""
    return psi, omega


# ---------------------------------------------------------------------------
# exakte loesung und ein lauf
# ---------------------------------------------------------------------------

def exakte_loesung(gebiet, t, nu):
    X, Y = gebiet.gitter()
    psi = np.sin(X) * np.sin(Y) * np.exp(-2.0 * nu * t)
    return psi, 2.0 * psi


def taylor_green_lauf(n, Re, t_end, dt):
    """rechnet den taylor-green-wirbel mit rk4 aus dem loeser bis t_end"""
    #R = 0.5 macht D = 1, damit ist nu = U_inf * D / Re = 1/Re. r_max wird nicht benutzt
    cfg = Config(R=0.5, r_max=2.0, U_inf=1.0, Re=Re, n_xi=n, n_theta=n, dt=dt)
    gebiet = KartesischesGebiet(n, cfg)
    ops = operatoren_periodisch(gebiet)
    löser = PoissonPeriodisch(gebiet)

    #startzustand: omega exakt, psi passend dazu aus der poisson-gleichung - das setzt rk4
    #voraus, genauso wie main.Anfangsbedingungen es fuer den zylinder macht
    _, omega = exakte_loesung(gebiet, 0.0, cfg.nu)
    psi = löser.löse(omega)

    t = 0.0
    while t < t_end - 1e-12:
        schritt = min(dt, t_end - t)
        psi, omega = rk4(psi, omega, cfg, löser, ops, schritt, gebiet,
                         randbedingungen=ohne_randbedingungen)
        t += schritt
    return gebiet, psi, omega, t


def l2_fehler(feld, exakt):
    return float(np.linalg.norm(feld - exakt) / np.linalg.norm(exakt))


def ordnungen(fehler):
    #die gitterweite halbiert sich je stufe
    return [np.log2(fehler[k] / fehler[k + 1]) for k in range(len(fehler) - 1)]


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

GITTER = (16, 32, 64)
NU = 0.1            #ueber Re = 1/nu = 10
T_ENDE = 0.5
DT = 0.005          #fest fuer alle gitter, sonst steckt im ortsfehler auch ein zeitfehler


def main():
    p = Pruefung()

    # --- 1) der ersatz-poisson-loeser loest wirklich das diskrete system -----
    p.abschnitt("poisson-loeser (periodisch)")
    gebiet = KartesischesGebiet(32, Config(R=0.5, r_max=2.0, U_inf=1.0, Re=10.0,
                                           n_xi=32, n_theta=32, dt=DT))
    ops = operatoren_periodisch(gebiet)
    omega = np.random.default_rng(0).standard_normal((32, 32))
    omega -= omega.mean()       #periodisch nur loesbar, wenn der mittelwert verschwindet
    psi = PoissonPeriodisch(gebiet).löse(omega)
    residuum = np.abs(gebiet.unflatten(ops["L"] @ gebiet.flatten(psi)) + omega).max()
    p.pruefe(residuum < 1e-10, f"residuum |L psi + omega| = {residuum:.1e} (soll ~ 0)")

    # --- 2) advektionsterm an einer loesung, die nicht beltrami ist ----------
    #der taylor-green-fall prueft den advektionsterm nicht, dort ist er identisch null.
    #deshalb hier psi und omega unabhaengig voneinander vorgeben und berechne_rhs mit
    #nu = 0 gegen die analytische ableitung stellen:
    #  psi = sin x sin y   ->  u =  dpsi/dy =  sin x cos y,  v = -dpsi/dx = -cos x sin y
    #  omega = sin 2x cos 3y
    #  u.grad(omega) = 2 sin x cos y cos 2x cos 3y + 3 cos x sin y sin 2x sin 3y
    p.abschnitt("advektionsterm (berechne_rhs, nu = 0)")
    fehler = []
    for n in GITTER:
        cfg = Config(R=0.5, r_max=2.0, U_inf=1.0, Re=np.inf, n_xi=n, n_theta=n, dt=DT)
        g = KartesischesGebiet(n, cfg)
        X, Y = g.gitter()
        psi = np.sin(X) * np.sin(Y)
        omega = np.sin(2 * X) * np.cos(3 * Y)
        advektion_exakt = (2 * np.sin(X) * np.cos(Y) * np.cos(2 * X) * np.cos(3 * Y)
                           + 3 * np.cos(X) * np.sin(Y) * np.sin(2 * X) * np.sin(3 * Y))
        #mit nu = 0 ist die rechte seite genau die negative advektion
        advektion = -berechne_rhs(omega, psi, g, cfg, operatoren_periodisch(g))
        fehler.append(l2_fehler(advektion, advektion_exakt))
    ord_adv = ordnungen(fehler)
    print("       " + ",  ".join(f"n = {n}: {e:.3e}" for n, e in zip(GITTER, fehler)))
    p.pruefe(1.9 < min(ord_adv) and max(ord_adv) < 2.1,
             f"konvergenzordnung {ord_adv[0]:.2f} / {ord_adv[1]:.2f} (soll 2)")

    #die geschwindigkeit, aus der die advektion gebildet wird, muss selbst stimmen
    g = KartesischesGebiet(64, Config(R=0.5, r_max=2.0, U_inf=1.0, Re=10.0,
                                      n_xi=64, n_theta=64, dt=DT))
    X, Y = g.gitter()
    u, v = geschwindigkeit(np.sin(X) * np.sin(Y), g, operatoren_periodisch(g))
    fehler_u = max(l2_fehler(u, np.sin(X) * np.cos(Y)), l2_fehler(v, -np.cos(X) * np.sin(Y)))
    p.pruefe(fehler_u < g.dxi**2, f"geschwindigkeit aus psi: fehler {fehler_u:.2e} "
                                  f"(schranke h^2 = {g.dxi**2:.2e})")

    # --- 3) taylor-green mit reibung: ordnung im ort -------------------------
    p.abschnitt(f"taylor-green, nu = {NU}, t = {T_ENDE}, dt = {DT} fest")
    fehler = []
    for n in GITTER:
        gebiet, psi, omega, t = taylor_green_lauf(n, 1.0 / NU, T_ENDE, DT)
        _, omega_exakt = exakte_loesung(gebiet, t, NU)
        fehler.append(l2_fehler(omega, omega_exakt))
        #der feste zeitschritt muss auf jedem gitter innerhalb der stabilitaetsgrenze liegen
        u, v = geschwindigkeit(psi, gebiet, operatoren_periodisch(gebiet))
        dt_grenze = cfl_zeitschritt(u, v, gebiet, gebiet.cfg)
        print(f"       n = {n:3d}: rel. L2-fehler omega = {fehler[-1]:.3e}   "
              f"(cfl-grenze dt = {dt_grenze:.4f})")
        p.pruefe(DT < dt_grenze, f"n = {n:3d}: dt = {DT} unter der cfl-grenze {dt_grenze:.4f}")

    ord_tg = ordnungen(fehler)
    p.pruefe(1.9 < min(ord_tg) and max(ord_tg) < 2.1,
             f"konvergenzordnung {ord_tg[0]:.2f} / {ord_tg[1]:.2f} (soll 2)")

    #die abklingrate der exakten loesung ist 2 nu. diskret weicht sie um O(h^2) ab - das
    #ist derselbe fehler wie oben, hier nur als physikalische groesse gelesen
    gebiet, psi, omega, t = taylor_green_lauf(64, 1.0 / NU, T_ENDE, DT)
    _, omega_start = exakte_loesung(gebiet, 0.0, NU)
    rate = -np.log(np.linalg.norm(omega) / np.linalg.norm(omega_start)) / t
    p.pruefe(abs(rate / (2.0 * NU) - 1.0) < 0.01,
             f"abklingrate n = 64: {rate:.6f} gegen 2 nu = {2 * NU} "
             f"({100 * (rate / (2 * NU) - 1):+.2f} %)")

    # --- 4) ohne reibung: stationaer, enstrophie erhalten --------------------
    #bei nu = 0 ist die loesung stationaer. jede aenderung ist dann rein numerisch - so
    #weist man kuenstliche diffusion nach (upwind 1. ordnung wuerde hier sofort auffallen)
    p.abschnitt("taylor-green ohne reibung (nu = 0, stationaere loesung)")
    gebiet, psi, omega, t = taylor_green_lauf(64, np.inf, T_ENDE, DT)
    _, omega_exakt = exakte_loesung(gebiet, t, 0.0)
    fehler_stat = l2_fehler(omega, omega_exakt)
    p.pruefe(fehler_stat < 1e-12, f"rel. L2-fehler nach {int(T_ENDE / DT)} schritten: "
                                  f"{fehler_stat:.2e} (soll maschinengenau)")

    enstrophie = 0.5 * np.sum(omega**2) * gebiet.dxi * gebiet.dtheta
    enstrophie_start = 0.5 * np.sum(omega_exakt**2) * gebiet.dxi * gebiet.dtheta
    drift = abs(enstrophie / enstrophie_start - 1.0)
    p.pruefe(drift < 1e-12, f"enstrophie-drift: {drift:.2e} (soll 0, sonst numerische diffusion)")

    return p.fazit("taylor-green")


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
