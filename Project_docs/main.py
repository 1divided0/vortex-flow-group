import numpy as np
import time
import os
import sys

from config import Config
from domain import Domain
from Poisson import PoissonSolver
from Randbedingungen import apply_bc
from zeitintegration import build_alle_operatoren,  geschwindigkeit, cfl_zeitschritt, rk4

def Anfangsbedingungen(domain, cfg, stoerung=0.5):

    R_grid = domain.r[:, None]
    Theta_grid = domain.theta[None, :]
    psi = cfg.U_inf * np.sin(Theta_grid) * (R_grid - cfg.R**2 / R_grid)
    omega = np.zeros_like(psi)

    # gewollte kleine stoerung: gitter und randbedingungen sind spiegelsymmetrisch,
    # ohne stoerung bleibt die loesung symmetrisch und es gibt keine wirbelstrasse.
    # ein kleiner wirbel leicht oberhalb der achse hinter dem zylinder bricht die
    # symmetrie kontrolliert (stoerung=0 schaltet das ab)
    X, Y = domain.kartesisch()
    x0, y0, breite = 1.5 * cfg.D, 0.3 * cfg.D, 0.1 * cfg.D**2
    omega += stoerung * np.exp(-((X - x0)**2 + (Y - y0)**2) / breite)

    return apply_bc(psi, omega, domain)

def  eine_Schleife(cfg, t_end, max_steps = None, Snapshotrange=50, verbose=True, t_speicher=0.0):
    # t_speicher: snapshots erst ab dieser zeit sichern. der anlauf bis zur
    # periodischen abloesung (bei Re=100 etwa t=50) wird fuer die auswertung
    # nicht gebraucht und wuerde die snapshot-datei nur unnoetig aufblaehen
    domain = Domain(cfg)
    löser = PoissonSolver(domain)
    ops = build_alle_operatoren(domain)

    psi, omega = Anfangsbedingungen(domain, cfg)

    snapshots = {"t": [], "psi": [], "omega": []}
    t = 0.0
    step = 0    # int, sonst bricht die ausgabe mit {step:6d} beim 100. schritt ab
    t_start_uhr = time.time()

    while t < t_end:
        if max_steps is not None and step >= max_steps:
            break

        u_r, u_theta = geschwindigkeit(psi, domain, ops)
        dt = cfl_zeitschritt(u_r, u_theta, domain, cfg)
        dt = min(dt, cfg.dt, t_end - t)

        psi, omega = rk4(psi, omega, cfg, löser, ops, dt, domain)

        t += dt
        step += 1

        if step % Snapshotrange == 0 and t >= t_speicher:
            snapshots["t"].append(t)
            snapshots["psi"].append(psi.copy())
            snapshots["omega"].append(omega.copy())

        if verbose and step % 1000 == 0:
            elapsed = time.time() - t_start_uhr
            print(f"Schritt {step:6d} t={t: .5f} dt={dt: .2e} "
                  f"max|omega|={np.max(np.abs(omega)):.3f}  ({elapsed:.1f}s Rechenzeit)")

    if verbose:
        elapsed = time.time() - t_start_uhr
        print(f"Fertig: {step} Schritte, t={t:.5f}, Rechenzeit={elapsed:.1f}s")
    return domain, snapshots, psi, omega 



# TEstt_bereich


if __name__ == "__main__":
    # Testkonfiguration mit bewusst KLEINEM Gitter und kurzer Simulations-
    # zeit -- reiner Funktionstest, dass die komplette Kette (Domain,
    # Solver, Zeitintegration) durchlaeuft. Fuer eine physikalisch
    # aussagekraeftige Simulation (z.B. um die Karmansche Wirbelstrasse
    # bei Re~100 zu sehen) braucht es ein feineres Gitter und eine viel
    # laengere Simulationszeit (mehrere Ablösezyklen) -- das dauert
    # entsprechend deutlich laenger und sollte lokal, nicht als
    # Schnelltest, laufen.
    # dt ist nur die obergrenze, der cfl-schritt bestimmt den tatsaechlichen schritt.
    #
    #   python main.py        schnelltest (kleines gitter, ~1000 schritte, wenige sekunden)
    #   python main.py lang   produktionslauf fuer visualisation.py: wirbelstrasse bei
    #         cd                Re=100, gitter 80x160, t=0..100, snapshots ab t=60
    #                         (dauert ca. 2-3 minuten, datei ca. 40 MB)
    if len(sys.argv) > 1 and sys.argv[1] == "lang":
        cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=80, n_theta=160, dt=0.05, cfl_target=0.5)
        domain, snapshots, psi_final, omega_final = eine_Schleife(
            cfg, t_end=100.0, Snapshotrange=20, verbose=True, t_speicher=60.0
        )
    else:
        cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=50, n_theta=100, dt=0.05, cfl_target=0.5)
        domain, snapshots, psi_final, omega_final = eine_Schleife(
            cfg, t_end=12.0, Snapshotrange=10, verbose=True
        )
 
    print()
    print(f"Anzahl gespeicherter Snapshots: {len(snapshots['t'])}")
    print(f"psi an der Wand am Ende exakt 0: {np.allclose(psi_final[domain.i_wall, :], 0.0, atol=1e-10)}")
    print(f"alle Werte endlich: {np.all(np.isfinite(psi_final)) and np.all(np.isfinite(omega_final))}")
 
    # Snapshots auf die Festplatte sichern, damit postprocessing.py sie
    # spaeter laden kann, OHNE die Simulation nochmal laufen lassen zu
    # muessen (siehe Diskussion: Simulieren und Visualisieren als zwei
    # getrennte Phasen)
    skript_ordner = os.path.dirname(os.path.abspath(__file__))
    pfad = os.path.join(skript_ordner, "simulation_snapshots.npz")

    # die config wird mitgespeichert, damit visualisation.py das gitter exakt
    # rekonstruieren kann. float32 halbiert die dateigroesse und reicht fuer
    # die auswertung voellig aus
    np.savez(
    pfad,
    t=np.array(snapshots["t"]),
    psi=np.array(snapshots["psi"], dtype=np.float32),
    omega=np.array(snapshots["omega"], dtype=np.float32),
    R=cfg.R, r_max=cfg.r_max, U_inf=cfg.U_inf, Re=cfg.Re,
    n_xi=cfg.n_xi, n_theta=cfg.n_theta, dt=cfg.dt, cfl_target=cfg.cfl_target,
    )
    
    print("Snapshots gespeichert in simulation_snapshots.npz")
                
        




    