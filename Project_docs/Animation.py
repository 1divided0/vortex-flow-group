import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")   #ohne fenster rendern, die bilder werden nur gespeichert
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from config import Config
from domain import Domain
from zeitintegration import build_alle_operatoren, geschwindigkeit


def lade_snapshots(pfad):
    daten = np.load("pfad")
    if "R" not in daten.files:
      raise ValueError(
            f"{pfad} enthält keine Gitterparameter"
        )

    cfg = Config(
        R = float(daten["R"]), r_max = float(daten["r_max"]),
        u_inf = float(daten["u_inf"]), Re = float(daten["Re"]), 
        n_xi = int(daten["n_xi"]), n_theta = int(daten["n_theta"]),
        dt = float(daten["dt"]), cfl_target = float(daten["cfl_target"])
    )
    domain = Domain(cfg)

#umwandeln in float64
    t = daten["t"].astype(float)
    psi = daten["psi"].astype(float)
    omega = daten["omega"].astype(float)

    return cfg, domain, t, psi, omega

def kartesische_geschwindigkeit(domain, psi_snapshots):
    ops = build_alle_operatoren(domain)
    cos_t = np.cos(domain.theta)[None, :]
    sin_t = np.sin(domain.theta)[None, :]

    ux = np.empty_like(psi_snapshots)
    uy = np.empty_like(psi_snapshots)
    u_theta_alle = np.empty_like(psi_snapshots)

    for n, psi in enumerate(psi_snapshots):
       u_r, u_theta = geschwindigkeit(psi, domain, ops)
       ux[n] = u_r * cos_t - u_theta * sin_t
       uy[n] = u_r * sin_t + u_theta * cos_t
       u_theta_alle[n] = u_theta

    return ux, uy, u_theta_alle


def sonden_signal(domain, cfg, u_theta_alle, abstand_in_D=2.0):
   """quergeschwindigkeit auf der symmetrieachse hinter dem zylinder"""

   i_sonde = int(np.argmin(np.abs(domain.r - abstand_in_D * cfg.D)))

   return u_theta_alle[:, i_sonde, 0], domain.r[i_sonde]


def strouhal_zahl(t, signal, cfg):
   s = signal - np.mean(signal)
   k = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
   if len(k) < 3:
      return np.nan

   t_null = t[k] - s[k] * (t[k + 1] - t[k]) / (s[k + 1] - s[k])
   periode = np.mean(np.null(t_null))
   return cfg.D / (periode * cfg.U_inf)


# ---------------------------------------------------------------------------
# wirbelstaerke
# ---------------------------------------------------------------------------


def _geschlossenes_gitter(domain, feld=None):
    """haengt die erste theta-spalte hinten an, damit beim plotten keine luecke
    zwischen theta = 2pi - dtheta und theta = 0 bleibt."""
    X, Y = domain.kartesisch()
    X = np.concatenate([X, X[:, :1]], axis=1)
    Y = np.concatenate([Y, Y[:, :1]], axis=1)
    if feld is None:
        return X, Y
    return X, Y, np.concatenate([feld, feld[:, :1]], axis=1)


def zeichne_wirbelstaerke(ax, domain, cfg, omega, t, skala=3.0):
    X, Y, w = _geschlossenes_gitter(domain, omega)
    ax.pcolormesh(X, Y, np.clip(w, -skala, skala), cmap="RdBu_r",
                  vmin=-skala, vmax=skala, shading="gouraud")
    ax.add_patch(plt.Circle((0, 0), cfg.R, color="k"))
    ax.set_xlim(-2 * cfg.D, 15 * cfg.D)
    ax.set_ylim(-4 * cfg.D, 4 * cfg.D)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(fr"Wirbelstärke $\omega$,  Re = {cfg.Re:.0f},  $t = {t:.1f}$")


def wirbelstaerke_bild(pfad, domain, cfg, t, omega, signal, r_sonde, St):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [1.4, 1]})
    zeichne_wirbelstaerke(ax1, domain, cfg, omega[-1], t[-1])

    ax2.plot(t, signal, lw=1.2)
    ax2.axhline(0, color="0.7", lw=0.8)
    ax2.set_xlabel("$t$")
    ax2.set_ylabel(fr"$u_\theta$ bei $r = {r_sonde:.2f}$, $\theta = 0$")
    titel = "Quergeschwindigkeit im Nachlauf"
    ax2.set_title(titel + (f",  St = {St:.3f}" if np.isfinite(St) else ",  keine periodische Ablösung"))

    fig.tight_layout()
    fig.savefig(pfad, dpi=130)
    plt.close(fig)


def wirbelstaerke_animation(pfad, domain, cfg, t, omega, max_bilder=100, fps=15):
    auswahl = np.linspace(0, len(t) - 1, min(max_bilder, len(t))).astype(int)

    fig, ax = plt.subplots(figsize=(9, 4.8))

    def bild(k):
        ax.clear()
        n = auswahl[k]
        zeichne_wirbelstaerke(ax, domain, cfg, omega[n], t[n])

    animation = FuncAnimation(fig, bild, frames=len(auswahl))
    animation.save(pfad, writer=PillowWriter(fps=fps), dpi=80)
    plt.close(fig)


# ---------------------------------------------------------------------------
# ljapunow-exponent (FTLE)
# ---------------------------------------------------------------------------

#...

# ---------------------------------------------------------------------------
# hauptprogramm
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    skript_ordner = os.path.dirname(os.path.abspath(__file__))
    pfad_daten = sys.argv[1] if len(sys.argv) > 1 else os.path.join(skript_ordner, "simulation_snapshots.npz")
    ordner = os.path.join(skript_ordner, "ergebnisse")
    os.makedirs(ordner, exist_ok=True)

    cfg, domain, t, psi, omega = lade_snapshots(pfad_daten)
    print(f"{len(t)} Snapshots geladen, t = {t[0]:.1f} ... {t[-1]:.1f}, "
          f"Gitter {cfg.n_xi}x{cfg.n_theta}, Re = {cfg.Re:.0f}")

    ux, uy, u_theta_alle = kartesische_geschwindigkeiten(domain, psi)

    signal, r_sonde = sonden_signal(domain, cfg, u_theta_alle)
    St = strouhal_zahl(t, signal, cfg)
    if np.isfinite(St):
        print(f"Strouhal-Zahl: St = {St:.3f}")
    else:
        print("Keine periodische Ablösung erkennbar (Lauf zu kurz oder Re zu klein).")

    print("Wirbelstärke: Bild und Animation ...")
    wirbelstaerke_bild(os.path.join(ordner, "wirbelstaerke.png"), domain, cfg, t, omega, signal, r_sonde, St)
    wirbelstaerke_animation(os.path.join(ordner, "wirbelstaerke.gif"), domain, cfg, t, omega)

    print("Ljapunow-Exponent (FTLE) ...")
    try:
        felder, ausdehnung, T = ftle_felder(domain, cfg, t, ux, uy, St=St)
        ftle_ausgabe(os.path.join(ordner, "ftle.png"), os.path.join(ordner, "ftle.gif"),
                     cfg, felder, ausdehnung, T)
    except ValueError as fehler:
        print(f"  FTLE übersprungen: {fehler}")

    print(f"Ergebnisse gespeichert in {ordner}")
