"""
method of manufactured solutions auf dem ECHTEN log-polaren gitter.

  python Validierung/mms_logpolar.py

warum es diesen test zusaetzlich zu taylor_green.py braucht: dort ist r = 1 und beide
richtungen sind periodisch. genau dadurch sind die metrik 1/r^2, die einseitigen randzeilen
von build_D_xi/build_D2_xi und die thom-formel abgeschaltet - der taylor-green-test kann
ueber sie also nichts sagen. hier laeuft dieselbe rechte seite (berechne_rhs) auf dem
gitter der zylinderrechnung.

das prinzip von MMS: statt eine loesung zu suchen, gibt man sich eine beliebige glatte
funktion VOR, setzt sie in die diskreten operatoren ein und vergleicht mit der analytisch
differenzierten. psi und omega muessen dabei keine physikalische loesung sein und auch die
randbedingungen nicht erfuellen - geprueft wird allein, ob der code das ableitet, was er
ableiten soll.
"""

import numpy as np

from pruefung import Pruefung      #setzt den suchpfad auf Project_docs

from gitter import Config, Domain
from loeser import apply_wall_bc, berechne_rhs, build_alle_operatoren

GITTER = (41, 81, 161)     #dxi halbiert sich exakt


def _domain(n_xi):
    cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0,
                 n_xi=n_xi, n_theta=2 * (n_xi - 1), dt=1e-3)
    return cfg, Domain(cfg)


def herstellte_loesung(dom, nu):
    """
    psi = sin(2 theta) exp(-xi),  omega = cos(3 theta) exp(-xi/2)

    omega ist bewusst KEIN vielfaches von psi: sonst verschwaende der advektionsterm
    identisch (beltrami) und bliebe ungeprueft - genau die luecke des taylor-green-falls.
    zurueckgegeben wird auch die analytisch gebildete rechte seite
        domega/dt = -(1/r)(u_r domega/dxi + u_theta domega/dtheta) + nu (1/r^2) lap(omega)
    """
    XI = dom.xi[:, None]
    TH = dom.theta[None, :]
    r = dom.r[:, None]

    psi = np.sin(2 * TH) * np.exp(-XI)
    omega = np.cos(3 * TH) * np.exp(-0.5 * XI)

    psi_xi = -psi
    psi_theta = 2 * np.cos(2 * TH) * np.exp(-XI)
    omega_xi = -0.5 * omega
    omega_theta = -3 * np.sin(3 * TH) * np.exp(-0.5 * XI)

    u_r = psi_theta / r
    u_theta = -psi_xi / r

    advektion = (u_r * omega_xi + u_theta * omega_theta) / r
    diffusion = nu * (0.25 * omega - 9.0 * omega) / r**2     #lap = (om_xixi + om_thth)/r^2

    return psi, omega, -advektion + diffusion


def ordnungen(fehler):
    return [np.log2(fehler[k] / fehler[k + 1]) for k in range(len(fehler) - 1)]


def main():
    p = Pruefung()

    # --- 1) berechne_rhs im inneren, mit metrik ------------------------------
    p.abschnitt("berechne_rhs auf dem log-polaren gitter (metrik, advektion, diffusion)")
    fehler = []
    for n_xi in GITTER:
        cfg, dom = _domain(n_xi)
        psi, omega, rhs_exakt = herstellte_loesung(dom, cfg.nu)
        rhs = berechne_rhs(omega, psi, dom, cfg, build_alle_operatoren(dom))
        #die beiden randzeilen sind bedeutungslos, dort setzt apply_bc omega
        innen = slice(1, -1)
        fehler.append(float(np.max(np.abs(rhs[innen] - rhs_exakt[innen]))
                            / np.max(np.abs(rhs_exakt[innen]))))
        print(f"       n_xi = {n_xi:3d}: rel. max-fehler = {fehler[-1]:.3e}")
    ord_rhs = ordnungen(fehler)
    p.pruefe(1.9 < min(ord_rhs) and max(ord_rhs) < 2.1,
             f"konvergenzordnung {ord_rhs[0]:.2f} / {ord_rhs[1]:.2f} (soll 2)")

    # --- 2) advektion und diffusion einzeln ---------------------------------
    #wenn 1) schieflaeuft, zeigt das hier, an welchem der beiden terme es liegt
    p.abschnitt("dieselbe loesung, terme einzeln")
    for name, Re_wert in (("nur advektion (nu = 0)", np.inf), ("nur diffusion (psi = 0)", 100.0)):
        fehler = []
        for n_xi in GITTER:
            cfg, dom = _domain(n_xi)
            cfg.Re = Re_wert
            psi, omega, _ = herstellte_loesung(dom, cfg.nu)
            if name.startswith("nur diffusion"):
                psi = np.zeros_like(psi)
            _, _, rhs_exakt = herstellte_loesung(dom, cfg.nu)
            if name.startswith("nur diffusion"):
                #ohne psi bleibt von der exakten rechten seite nur die diffusion
                r = dom.r[:, None]
                rhs_exakt = cfg.nu * (0.25 - 9.0) * omega / r**2
            rhs = berechne_rhs(omega, psi, dom, cfg, build_alle_operatoren(dom))
            innen = slice(1, -1)
            fehler.append(float(np.max(np.abs(rhs[innen] - rhs_exakt[innen]))
                                / np.max(np.abs(rhs_exakt[innen]))))
        o = ordnungen(fehler)
        p.pruefe(1.9 < min(o) and max(o) < 2.1,
                 f"{name:24s}: {fehler[0]:.2e} / {fehler[1]:.2e} / {fehler[2]:.2e}, "
                 f"ordnung {o[0]:.2f} / {o[1]:.2f}")

    # --- 3) thom-formel: ist sie wirklich 1. ordnung? -----------------------
    #psi = xi^2 exp(-xi) sin(3 theta) erfuellt an der wand psi = 0 und dpsi/dxi = 0,
    #also genau das, was die thom-formel voraussetzt. exakt ist dort
    #   omega_w = -(1/R^2) d2psi/dxi2 = -(2/R^2) sin(3 theta)
    #die gemessene ordnung belegt die aussage "thom ist 1. ordnung in omega_w" -
    #sie ist damit die einzige stelle des loesers, die nicht 2. ordnung ist
    p.abschnitt("thom-formel an der wand (erwartet: 1. ordnung)")
    fehler = []
    for n_xi in GITTER:
        cfg, dom = _domain(n_xi)
        XI = dom.xi[:, None]
        TH = dom.theta[None, :]
        psi = (XI**2) * np.exp(-XI) * np.sin(3 * TH)
        omega = np.zeros_like(psi)
        apply_wall_bc(psi, omega, dom)
        exakt = -(2.0 / cfg.R**2) * np.sin(3 * dom.theta)
        fehler.append(float(np.max(np.abs(omega[dom.i_wall] - exakt)) / np.max(np.abs(exakt))))
        print(f"       n_xi = {n_xi:3d}: rel. fehler omega_wand = {fehler[-1]:.3e}")
    ord_thom = ordnungen(fehler)
    p.pruefe(0.9 < min(ord_thom) and max(ord_thom) < 1.2,
             f"konvergenzordnung {ord_thom[0]:.2f} / {ord_thom[1]:.2f} (soll 1 - "
             f"das ist die bekannte schwaeche der thom-formel, kein fehler)")

    return p.fazit("MMS log-polar")


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
