import numpy as np
import time 
import os

from config import Config
from domain import Domain
from Poisson import PoissonSolver
from Randbedingungen import apply_bc
from zeitintegration import build_alle_operatoren,  geschwindigkeit, cfl_zeitschritt, rk4

def Anfangsbedingungen(domain, cfg):

    R_grid = domain.r[:, None]
    Theta_grid = domain.theta[None, :]
    psi = cfg.U_inf * np.sin(Theta_grid) * (R_grid - cfg.R**2 / R_grid)
    omega = np.zeros_like(psi)
    return apply_bc(psi, omega, domain)

def  eine_Schleife(cfg, t_end, max_steps = None, Snapshotrange=50, verbose=True):
    domain = Domain(cfg)
    löser = PoissonSolver(domain)
    ops = build_alle_operatoren(domain)

    psi, omega = Anfangsbedingungen(domain, cfg)

    snapshots = {"t": [], "psi": [], "omega": []}
    t = 0.0
    step = 0.0
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

        if step % Snapshotrange == 0:
            snapshots["t"].append(t)
            snapshots["psi"].append(psi.copy())
            snapshots["omega"].append(omega.copy())

            if verbose and step % 100 == 0:
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
    cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=50, n_theta=100, dt=1e-3, cfl_target=0.5)
 
    domain, snapshots, psi_final, omega_final = eine_Schleife(
        cfg, t_end=0.05, Snapshotrange=10, verbose=True
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

    np.savez(
    pfad,
    t=np.array(snapshots["t"]),
    psi=np.array(snapshots["psi"]),
    omega=np.array(snapshots["omega"]),
    )
    
    print("Snapshots gespeichert in simulation_snapshots.npz")
                
        




    