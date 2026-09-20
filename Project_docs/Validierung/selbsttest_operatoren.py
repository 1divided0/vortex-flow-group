"""
selbsttest der differenzenmatrizen aus operatoren.py: jede matrix wird auf eine funktion
losgelassen, deren ableitung analytisch bekannt ist.

verglichen wird immer das GANZE feld einschliesslich der randzeilen - wer die raender
ausklammert, uebersieht genau die fehler, die dort sitzen (in diesem projekt waren zwei
von sechs gefundenen operator-fehlern randfehler).

  python Validierung/selbsttest_operatoren.py
"""

import numpy as np
import scipy.sparse as sp

from pruefung import Pruefung      #setzt den suchpfad auf Project_docs

from gitter import Config, Domain
from operatoren import (build_D_xi, build_D2_xi, build_D_theta, build_D2_theta,
                        build_upwind_xi, build_upwind_theta, build_laplacian)


def main():
    p = Pruefung()
    cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=80, n_theta=160, dt=1e-3)
    dom = Domain(cfg)
    I_xi = sp.identity(dom.n_xi, format="csr")
    I_theta = sp.identity(dom.n_theta, format="csr")

    # --- ableitungen in xi: f = xi^3 -> f' = 3 xi^2, f'' = 6 xi ---------------
    #die matrizen wirken nur auf die xi-richtung (groesse n_xi). fuer ein volles,
    #geflachtes 2d-feld muessen sie per kronecker-produkt mit der identitaet in
    #theta-richtung auf die volle groesse gebracht werden - dasselbe prinzip wie
    #in build_laplacian und build_alle_operatoren
    p.abschnitt("ableitungen in xi (f = xi^3, einseitige randformeln)")
    f = np.tile(dom.xi[:, None]**3, (1, dom.n_theta))

    df = dom.unflatten(sp.kron(build_D_xi(dom), I_theta, format="csr") @ dom.flatten(f))
    fehler = np.max(np.abs(df - 3.0 * dom.xi[:, None]**2))
    #zentrale differenz: fehler h^2/6 * f''' = h^2, einseitig am rand das doppelte
    p.pruefe(fehler < 3.0 * dom.dxi**2, f"D_xi:  max fehler {fehler:.2e} (schranke {3.0 * dom.dxi**2:.2e})")

    d2f = dom.unflatten(sp.kron(build_D2_xi(dom), I_theta, format="csr") @ dom.flatten(f))
    fehler = np.max(np.abs(d2f - 6.0 * dom.xi[:, None]))
    #beide formeln (zentral und einseitig) sind fuer kubische funktionen exakt,
    #uebrig bleibt nur rundungsfehler, durch 1/h^2 verstaerkt
    p.pruefe(fehler < 1e-6, f"D2_xi: max fehler {fehler:.2e} (fuer kubische funktionen exakt)")

    # --- ableitungen in theta: g = sin(theta), periodisch --------------------
    p.abschnitt("ableitungen in theta (g = sin(theta), periodisch)")
    g = np.tile(np.sin(dom.theta)[None, :], (dom.n_xi, 1))

    dg = dom.unflatten(sp.kron(I_xi, build_D_theta(dom), format="csr") @ dom.flatten(g))
    fehler = np.max(np.abs(dg - np.cos(dom.theta)[None, :]))
    p.pruefe(fehler < dom.dtheta**2, f"D_theta:  max fehler {fehler:.2e} (schranke {dom.dtheta**2:.2e})")

    d2g = dom.unflatten(sp.kron(I_xi, build_D2_theta(dom), format="csr") @ dom.flatten(g))
    fehler = np.max(np.abs(d2g + np.sin(dom.theta)[None, :]))
    p.pruefe(fehler < dom.dtheta**2, f"D2_theta: max fehler {fehler:.2e} (schranke {dom.dtheta**2:.2e})")

    # --- laplace mit metrik: h = xi^2 sin(theta) -----------------------------
    #nabla^2 h = (1/r^2) (h_xixi + h_thetatheta) = (1/r^2) (2 - xi^2) sin(theta)
    p.abschnitt("laplace-operator inkl. metrik 1/r^2 (h = xi^2 sin(theta))")
    h_feld = (dom.xi[:, None]**2) * np.sin(dom.theta)[None, :]
    lap = dom.unflatten(build_laplacian(dom) @ dom.flatten(h_feld))
    lap_exakt = dom.vorfaktor[:, None] * (2.0 - dom.xi[:, None]**2) * np.sin(dom.theta)[None, :]
    fehler = np.max(np.abs(lap - lap_exakt))
    #der fehler steckt allein im theta-anteil (in xi ist die funktion quadratisch und
    #damit exakt), verstaerkt um den groessten metrischen vorfaktor 1/R^2 an der wand
    schranke = dom.dtheta**2 * np.max(dom.xi**2) * np.max(dom.vorfaktor)
    p.pruefe(fehler < schranke, f"laplace: max fehler {fehler:.2e} (schranke {schranke:.2e})")

    # --- upwind-matrizen: zeilensumme muss exakt null sein -------------------
    #eine ableitungsmatrix muss ein konstantes feld auf null abbilden. genau das war
    #in den randzeilen der xi-upwind-matrizen einmal verletzt
    p.abschnitt("upwind-matrizen (zeilensummen)")
    for name, (fw, bw) in (("xi", build_upwind_xi(dom)), ("theta", build_upwind_theta(dom))):
        summe = max(np.abs(np.asarray(fw.sum(axis=1))).max(),
                    np.abs(np.asarray(bw.sum(axis=1))).max())
        p.pruefe(summe < 1e-12, f"upwind {name:5s}: max |zeilensumme| {summe:.1e} (soll 0)")

    return p.fazit("selbsttest operatoren")


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
