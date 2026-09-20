"""
selbsttest des poisson-loesers und der randbedingungen aus loeser.py.

abnahmetest: bei omega = 0 ist die potentialstroemung um den zylinder die exakte loesung.
sie ist damit der einzige fall, in dem das ZUSAMMENSPIEL aus laplace-operator, metrik,
dirichlet-randzeilen und wandbedingung an einer geschlossenen formel geprueft werden kann.
zusaetzlich wird die konvergenzordnung gemessen statt geglaubt, und der FFT-loeser gegen
die LR-zerlegung gestellt.

  python Validierung/selbsttest_loeser.py
"""

import numpy as np

from pruefung import Pruefung      #setzt den suchpfad auf Project_docs

from gitter import Config, Domain
from loeser import (PoissonSolver, PoissonSolverFFT, wähle_poisson_löser, build_poisson_matrix,
                    build_rhs, build_alle_operatoren, geschwindigkeit, potential_flow_psi)


def potentialstroemung(dom):
    """exakte loesung fuer omega = 0 auf dem ganzen gebiet: psi = U sin(theta) (r - R^2/r)"""
    return potential_flow_psi(dom, dom.r[:, None])


def main():
    p = Pruefung()

    # --- abnahmetest und konvergenzordnung ----------------------------------
    #n_xi - 1 verdoppelt sich, damit sich dxi exakt halbiert
    p.abschnitt("poisson-loeser gegen die potentialstroemung (omega = 0)")
    fehler = []
    for n_xi in (41, 81, 161):
        cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0,
                     n_xi=n_xi, n_theta=2 * (n_xi - 1), dt=1e-3)
        dom = Domain(cfg)
        psi = PoissonSolver(dom).löse(np.zeros((dom.n_xi, dom.n_theta)))
        fehler.append(np.max(np.abs(psi - potentialstroemung(dom))))
        print(f"       n_xi = {n_xi:3d}: max |psi - psi_exakt| = {fehler[-1]:.3e}")

    ordnung = [np.log2(fehler[0] / fehler[1]), np.log2(fehler[1] / fehler[2])]
    p.pruefe(1.8 < min(ordnung) and max(ordnung) < 2.2,
             f"konvergenzordnung {ordnung[0]:.2f} / {ordnung[1]:.2f} (soll 2)")

    # --- geschwindigkeit an der wand ----------------------------------------
    #reibungsfrei gilt an der zylinderwand analytisch u_theta = -2 U sin(theta).
    #das prueft die ableitung dpsi/dxi genau in der zeile, in der die einseitige
    #randformel steht - dort, wo spaeter auch die thom-formel ansetzt
    p.abschnitt("wandgeschwindigkeit der potentialstroemung")
    cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=81, n_theta=160, dt=1e-3)
    dom = Domain(cfg)
    psi = PoissonSolver(dom).löse(np.zeros((dom.n_xi, dom.n_theta)))
    _, u_theta = geschwindigkeit(psi, dom, build_alle_operatoren(dom))
    fehler_wand = np.max(np.abs(u_theta[dom.i_wall] + 2.0 * cfg.U_inf * np.sin(dom.theta)))
    p.pruefe(fehler_wand < 2e-3,
             f"81x160: max |u_theta + 2 U sin(theta)| = {fehler_wand:.2e} (schranke 2.0e-03)")

    # --- die beiden poisson-loeser muessen dasselbe system loesen ------------
    p.abschnitt("FFT-loeser gegen LR-zerlegung")
    for n_xi, n_theta in ((41, 80), (81, 160), (121, 240)):
        cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=n_xi, n_theta=n_theta, dt=1e-3)
        dom = Domain(cfg)
        omega = np.random.default_rng(0).standard_normal((n_xi, n_theta))

        psi_lr = PoissonSolver(dom).löse(omega)
        psi_fft = PoissonSolverFFT(dom).löse(omega)

        #residuum im urspruenglichen gleichungssystem, nicht nur vergleich der beiden:
        #zwei loeser koennen sich auch einig sein und trotzdem beide falsch liegen
        residuum = np.abs(build_poisson_matrix(dom) @ dom.flatten(psi_fft) - build_rhs(dom, omega)).max()
        abweichung = np.abs(psi_fft - psi_lr).max()
        groesse = np.abs(psi_lr).max()

        gewaehlt = type(wähle_poisson_löser(dom)).__name__
        p.pruefe(abweichung < 1e-8 * groesse and residuum < 1e-8,
                 f"{n_xi:3d}x{n_theta:<3d}: |psi_FFT - psi_LR| = {abweichung:.1e} "
                 f"(|psi| bis {groesse:.1f}), residuum = {residuum:.1e}, gewaehlt: {gewaehlt}")

    return p.fazit("selbsttest loeser")


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
