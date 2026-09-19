import numpy as np
import time
import os
import sys

from gitter import Config, Domain
from loeser import wähle_poisson_löser, apply_bc, build_alle_operatoren, geschwindigkeit, cfl_zeitschritt, rk4


#voreinstellungen fuer einzellaeufe, beschrieben in README_Presets.md.
#die serien des benchmarks stehen in benchmark/presets.py
EINZELLAEUFE = {
    "schnell": dict(
        beschreibung="funktionstest: kleines gitter, t = 0..12, wenige sekunden",
        cfg=Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=50, n_theta=100, dt=0.05, cfl_target=0.5),
        t_end=12.0, Snapshotrange=10, t_speicher=0.0,
    ),
    "lang": dict(
        beschreibung=("produktionslauf fuer Animation.py: Re = 100, 160x320, t = 0..100, "
                      "snapshots ab t = 60. ca. 15 min, ~83000 schritte, ~1660 snapshots "
                      "-> datei ca. 680 MB, Animation.py braucht dann mehrere GB RAM"),
        cfg=Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=160, n_theta=320, dt=0.05, cfl_target=0.5),
        t_end=100.0, Snapshotrange=20, t_speicher=60.0,
    ),
}


def Anfangsbedingungen(domain, cfg, löser, stoerung=0.5):
    # startzustand: potentialstroemung plus kleine stoerung in omega.
    # psi wird ueber poisson aus omega berechnet, damit psi und omega von anfang
    # an zueinander passen (das setzt rk4 voraus). bei omega = 0 liefert poisson
    # wegen des dirichlet-fernrands genau die diskrete potentialstroemung
    omega = np.zeros((domain.n_xi, domain.n_theta))

    # gewollte kleine stoerung: gitter und randbedingungen sind spiegelsymmetrisch,
    # ohne stoerung bleibt die loesung symmetrisch und es gibt keine wirbelstrasse.
    # ein kleiner wirbel leicht oberhalb der achse hinter dem zylinder bricht die
    # symmetrie kontrolliert (stoerung=0 schaltet das ab)
    X, Y = domain.kartesisch()
    x0, y0, breite = 1.5 * cfg.D, 0.3 * cfg.D, 0.1 * cfg.D**2
    omega += stoerung * np.exp(-((X - x0)**2 + (Y - y0)**2) / breite)

    psi = löser.löse(omega)
    return apply_bc(psi, omega, domain)


def eine_Schleife(cfg, t_end, max_steps=None, Snapshotrange=50, verbose=True, t_speicher=0.0, messung=None):
    # t_speicher: snapshots erst ab dieser zeit sichern. der anlauf bis zur
    # periodischen abloesung (bei Re=100 etwa t=50) wird fuer die auswertung
    # nicht gebraucht und wuerde die snapshot-datei nur unnoetig aufblaehen
    # messung: optionales objekt mit start(domain, löser, psi, omega) und
    # schritt(t, dt, psi, omega), z.b. aus benchmark.py. es darf psi und omega nur
    # lesen. ohne messung (None) laeuft die schleife exakt wie vorher
    domain = Domain(cfg)
    löser = wähle_poisson_löser(domain)
    ops = build_alle_operatoren(domain)

    psi, omega = Anfangsbedingungen(domain, cfg, löser)
    if messung is not None:
        messung.start(domain, löser, psi, omega)

    snapshots = {"t": [], "psi": [], "omega": []}
    t = 0.0
    step = 0    # int, sonst bricht die ausgabe mit {step:6d} beim 1000. schritt ab
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

        if messung is not None:
            messung.schritt(t, dt, psi, omega)

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


if __name__ == "__main__":
    # die voreinstellungen stehen oben in EINZELLAEUFE und in README_Presets.md.
    # dt ist nur die obergrenze, der cfl-schritt bestimmt den tatsaechlichen schritt.
    #
    #   python main.py          = python main.py schnell (funktionstest, wenige sekunden)
    #   python main.py lang     produktionslauf fuer Animation.py (ca. 15 minuten, datei ca. 680 MB)
    #
    # achtung: alle varianten schreiben nach simulation_snapshots.npz, der
    # schnelltest ueberschreibt also einen vorhandenen produktionslauf
    name = sys.argv[1] if len(sys.argv) > 1 else "schnell"
    if name not in EINZELLAEUFE:
        print(f"unbekannte voreinstellung '{name}'. vorhanden:")
        for n, p in EINZELLAEUFE.items():
            print(f"  {n:10s} {p['beschreibung']}")
        sys.exit(1)

    preset = EINZELLAEUFE[name]
    cfg = preset["cfg"]
    domain, snapshots, psi_final, omega_final = eine_Schleife(
        cfg, t_end=preset["t_end"], Snapshotrange=preset["Snapshotrange"],
        verbose=True, t_speicher=preset["t_speicher"]
    )

    print()
    print(f"Anzahl gespeicherter Snapshots: {len(snapshots['t'])}")
    print(f"psi an der Wand am Ende exakt 0: {np.allclose(psi_final[domain.i_wall, :], 0.0, atol=1e-10)}")
    print(f"alle Werte endlich: {np.all(np.isfinite(psi_final)) and np.all(np.isfinite(omega_final))}")

    # snapshots auf die festplatte sichern, damit Animation.py sie spaeter laden
    # kann, ohne die simulation nochmal laufen zu lassen (simulieren und
    # auswerten sind zwei getrennte phasen)
    skript_ordner = os.path.dirname(os.path.abspath(__file__))
    pfad = os.path.join(skript_ordner, "simulation_snapshots.npz")

    # die config wird mitgespeichert, damit Animation.py das gitter exakt
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
