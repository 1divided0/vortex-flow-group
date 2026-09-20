import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")   #ohne fenster rendern, die bilder werden nur gespeichert
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from gitter import Config, Domain
from loeser import build_alle_operatoren, geschwindigkeit
from signalauswertung import sonden_lage, werte_signal_aus


def lade_snapshots(pfad):
    #snapshot-datei aus main.py laden und daraus gitter und config rekonstruieren
    daten = np.load(pfad)
    if "R" not in daten.files:
        raise ValueError(f"{pfad} enthält keine Gitterparameter")

    cfg = Config(
        R = float(daten["R"]), r_max = float(daten["r_max"]),
        U_inf = float(daten["U_inf"]), Re = float(daten["Re"]),
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
    #u_x, u_y aus u_r, u_theta fuer alle snapshots
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
    """
    quergeschwindigkeit auf der symmetrieachse hinter dem zylinder.
    die lage wird zwischen zwei xi-zeilen interpoliert (sonden_lage) statt auf den
    naechsten gitterpunkt gerundet - damit misst diese sonde auf jedem gitter an
    derselben stelle, und an genau derselben wie die sonde des benchmarks
    """
    i, w = sonden_lage(domain, cfg, abstand_in_D)
    signal = (1.0 - w) * u_theta_alle[:, i, 0] + w * u_theta_alle[:, i + 1, 0]
    return signal, abstand_in_D * cfg.D


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


def wirbelstaerke_bild(pfad, domain, cfg, t, omega, signal, r_sonde, auswertung):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [1.4, 1]})
    zeichne_wirbelstaerke(ax1, domain, cfg, omega[-1], t[-1])

    ax2.plot(t, signal, lw=1.2)
    ax2.axhline(0, color="0.7", lw=0.8)
    #das ausgewertete periodenfenster sichtbar machen - die strouhal-zahl stammt
    #nur von dort und nicht vom ganzen signal
    St = auswertung["St"]
    if np.isfinite(St):
        ax2.axvspan(auswertung["t_fenster_start"], auswertung["t_fenster_ende"],
                    color="C1", alpha=0.12, lw=0, label="Auswertefenster")
        ax2.legend(fontsize=8, loc="lower right")
    ax2.set_xlabel("$t$")
    ax2.set_ylabel(fr"$u_\theta$ bei $r = {r_sonde:.2f}$, $\theta = 0$")
    titel = "Quergeschwindigkeit im Nachlauf"
    if np.isfinite(St):
        titel += (f",  St = {St:.4f},  Amplitude = {auswertung['amplitude']:.3f}"
                  f"  ({auswertung['n_perioden']} Perioden)")
    else:
        titel += f",  keine periodische Ablösung (Status {auswertung['status']})"
    ax2.set_title(titel, fontsize=10)

    fig.tight_layout()
    #dpi bestimmt allein die pixelzahl der abgabe-bilder und kostet keine rechenzeit:
    #200 dpi ergeben bei figsize (10, 8) genau 2000x1600 px
    fig.savefig(pfad, dpi=200)
    plt.close(fig)


def wirbelstaerke_animation(pfad, domain, cfg, t, omega, max_bilder=100, fps=15):
    auswahl = np.linspace(0, len(t) - 1, min(max_bilder, len(t))).astype(int)

    fig, ax = plt.subplots(figsize=(9, 4.8))

    def bild(k):
        ax.clear()
        n = auswahl[k]
        zeichne_wirbelstaerke(ax, domain, cfg, omega[n], t[n])

    animation = FuncAnimation(fig, bild, frames=len(auswahl))
    #dpi 80 ergab nur 720x384 px, 130 dpi sind 1170x624 px. gemessen an der
    #hauptloesung waechst das gif dadurch von etwa 7 auf 12.6 MB bei 100 bildern -
    #wird das zu gross, max_bilder senken statt dpi
    animation.save(pfad, writer=PillowWriter(fps=fps), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# ljapunow-exponent (FTLE)
# ---------------------------------------------------------------------------

def geschwindigkeit_bei(domain, cfg, ux, uy, x, y):
    """
    bilineare interpolation eines geschwindigkeitsfeldes vom (xi, theta)-gitter
    auf beliebige kartesische punkte.
    """
    r = np.hypot(x, y)
    theta = np.mod(np.arctan2(y, x), 2.0 * np.pi)
    xi = np.log(np.maximum(r, 1e-12) / cfg.R)

    #theta ist periodisch: der rechte nachbar von n_theta-1 ist wieder 0
    jt = theta / domain.dtheta
    j0 = np.floor(jt).astype(int) % domain.n_theta
    j1 = (j0 + 1) % domain.n_theta
    wj = jt - np.floor(jt)

    it = np.clip(xi / domain.dxi, 0.0, domain.n_xi - 1 - 1e-9)
    i0 = np.floor(it).astype(int)
    i1 = i0 + 1
    wi = it - i0

    def bilinear(F):
        return ((1 - wi) * (1 - wj) * F[i0, j0] + (1 - wi) * wj * F[i0, j1]
                + wi * (1 - wj) * F[i1, j0] + wi * wj * F[i1, j1])

    vx, vy = bilinear(ux), bilinear(uy)

    #im zylinder haftbedingung, ausserhalb des rechengebiets ungestoerte anstroemung
    innen = r < cfg.R
    aussen = r > cfg.r_max
    vx[innen], vy[innen] = 0.0, 0.0
    vx[aussen], vy[aussen] = cfg.U_inf, 0.0
    return vx, vy


def flussabbildung(domain, cfg, t, ux, uy, n_start, schritte, X0, Y0, richtung):
    """
    transportiert ein partikelgitter ueber `schritte` snapshot-intervalle.
    """
    X, Y = X0.copy(), Y0.copy()
    for k in range(schritte):
        n = n_start + richtung * k
        m = n + richtung
        h = t[m] - t[n]    #bei rueckwaerts negativ

        vx1, vy1 = geschwindigkeit_bei(domain, cfg, ux[n], uy[n], X, Y)
        vx2, vy2 = geschwindigkeit_bei(domain, cfg, ux[m], uy[m], X + h * vx1, Y + h * vy1)

        X += 0.5 * h * (vx1 + vx2)
        Y += 0.5 * h * (vy1 + vy2)
    return X, Y


def ljapunow_exponent(X, Y, delta, T):
    """
    sigma = 1/|T| * ln( sqrt( lambda_max(J^T J) ) )
    """
    J11 = np.gradient(X, delta, axis=1)   #dX/dx0
    J12 = np.gradient(X, delta, axis=0)   #dX/dy0
    J21 = np.gradient(Y, delta, axis=1)   #dY/dx0
    J22 = np.gradient(Y, delta, axis=0)   #dY/dy0

    a = J11**2 + J21**2
    b = J11 * J12 + J21 * J22
    d = J12**2 + J22**2
    lambda_max = 0.5 * (a + d) + np.sqrt((0.5 * (a - d))**2 + b**2)

    return np.log(np.sqrt(np.maximum(lambda_max, 1e-30))) / abs(T)


def _normieren(sigma):
    """skalierung auf [0, 1]: der median wird zu 0, die obersten 0.5 %
    werden gekappt, damit einzelne extremwerte die farbskala nicht dominieren."""
    werte = sigma[np.isfinite(sigma)]
    unten, oben = np.percentile(werte, [50, 99.5])
    return np.nan_to_num(np.clip((sigma - unten) / (oben - unten), 0.0, 1.0))


def ftle_felder(domain, cfg, t, ux, uy, T=8.0, delta_in_D=0.01, anzahl=30, St=np.nan):
    """berechnet FTLE-felder fuer mehrere startzeiten t0 ueber eine abloeseperiode.

    delta_in_D ist die aufloesung des *partikelgitters* und damit die bildaufloesung
    des FTLE-bildes - sie haengt nicht am CFD-gitter. mit 0.04 waren es 338x151
    datenpunkte, die imshow auf ~1100 px hochinterpoliert hat (faktor 3.3, sichtbar
    unscharf). 0.01 ergibt 1350x601 punkte, also 16-mal so viele partikel.
    kosten: gut 45 s je feld, mit anzahl = 30 also etwa 25 min. wird das zu lang,
    zuerst anzahl senken (weniger gif-phasen), erst danach delta_in_D erhoehen
    """
    D = cfg.D
    x0 = np.arange(-1.5 * D, 12.0 * D, delta_in_D * D)
    y0 = np.arange(-3.0 * D, 3.0 * D + 1e-9, delta_in_D * D)
    X0, Y0 = np.meshgrid(x0, y0)
    im_zylinder = np.hypot(X0, Y0) < cfg.R

    abstand = np.mean(np.diff(t))
    schritte = int(round(T / abstand))
    if 2 * schritte >= len(t):
        raise ValueError(
            f"fuer T = {T} werden {2 * schritte + 1} snapshots gebraucht, vorhanden "
            f"sind {len(t)}. laengeren lauf speichern oder T verkleinern."
        )

    #startzeiten: nur dort, wo vorwaerts und rueckwaerts genug snapshots liegen,
    #und moeglichst genau eine abloeseperiode, damit die animation sich wiederholt
    n_min, n_max = schritte, len(t) - 1 - schritte
    periode = D / (St * cfg.U_inf) if np.isfinite(St) else t[n_max] - t[n_min]
    t_bis = min(t[n_min] + periode, t[n_max])
    startzeiten = np.linspace(t[n_min], t_bis, anzahl, endpoint=np.isnan(St))

    felder = []
    for t_start in startzeiten:
        n0 = int(np.argmin(np.abs(t - t_start)))
        n0 = min(max(n0, n_min), n_max)

        Xv, Yv = flussabbildung(domain, cfg, t, ux, uy, n0, schritte, X0, Y0, +1)
        Xr, Yr = flussabbildung(domain, cfg, t, ux, uy, n0, schritte, X0, Y0, -1)

        sigma_vor = ljapunow_exponent(Xv, Yv, delta_in_D * D, T)
        sigma_rueck = ljapunow_exponent(Xr, Yr, delta_in_D * D, T)
        sigma_vor[im_zylinder] = np.nan
        sigma_rueck[im_zylinder] = np.nan

        felder.append((t[n0], _normieren(sigma_vor), _normieren(sigma_rueck)))
        print(f"  FTLE bei t0 = {t[n0]:.2f} berechnet")

    ausdehnung = [x0[0], x0[-1], y0[0], y0[-1]]
    return felder, ausdehnung, T


def _ftle_rgb(sigma_vor, sigma_rueck):
    """
    ueberlagert beide felder: weiss als hintergrund, vorwaerts (instabil)
    nimmt blau und gruen weg -> rot, rueckwaerts (stabil) nimmt rot und gruen
    weg -> blau.
    """
    R = 1.0 - sigma_rueck
    G = 1.0 - 0.85 * sigma_vor - 0.85 * sigma_rueck
    B = 1.0 - sigma_vor
    return np.clip(np.dstack([R, G, B]), 0.0, 1.0)


def zeichne_ftle(ax, cfg, feld, ausdehnung, T):
    t0, sigma_vor, sigma_rueck = feld
    ax.imshow(_ftle_rgb(sigma_vor, sigma_rueck), origin="lower",
              extent=ausdehnung, interpolation="bilinear")
    ax.add_patch(plt.Circle((0, 0), cfg.R, color="k"))
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(fr"Ljapunow-Exponent,  Re = {cfg.Re:.0f},  $t_0 = {t0:.1f}$,  $T = \pm{T:.0f}$")


def ftle_ausgabe(pfad_png, pfad_gif, cfg, felder, ausdehnung, T, fps=12):
    #fps 6 bei 12 phasen war eine 2-sekunden-schleife und ruckelte sichtbar.
    #30 phasen bei 12 fps ergeben 2.5 s je abloeseperiode und laufen rund
    fig, ax = plt.subplots(figsize=(10, 4.8))
    zeichne_ftle(ax, cfg, felder[0], ausdehnung, T)
    fig.tight_layout()
    fig.savefig(pfad_png, dpi=200)

    def bild(k):
        ax.clear()
        zeichne_ftle(ax, cfg, felder[k], ausdehnung, T)

    animation = FuncAnimation(fig, bild, frames=len(felder))
    animation.save(pfad_gif, writer=PillowWriter(fps=fps), dpi=130)
    plt.close(fig)

# ---------------------------------------------------------------------------
# hauptprogramm
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    #umlaute sollen auch dann funktionieren, wenn die ausgabe in eine datei umgeleitet
    #wird - windows faellt sonst auf cp1252 zurueck und bricht mit UnicodeEncodeError ab
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, OSError):
        pass

    skript_ordner = os.path.dirname(os.path.abspath(__file__))
    pfad_daten = sys.argv[1] if len(sys.argv) > 1 else os.path.join(skript_ordner, "simulation_snapshots.npz")
    ordner = os.path.join(skript_ordner, "ergebnisse")
    os.makedirs(ordner, exist_ok=True)

    cfg, domain, t, psi, omega = lade_snapshots(pfad_daten)
    print(f"{len(t)} Snapshots geladen, t = {t[0]:.1f} ... {t[-1]:.1f}, "
          f"Gitter {cfg.n_xi}x{cfg.n_theta}, Re = {cfg.Re:.0f}")

    ux, uy, u_theta_alle = kartesische_geschwindigkeit(domain, psi)

    signal, r_sonde = sonden_signal(domain, cfg, u_theta_alle)
    #dieselbe auswertung wie im benchmark: sie sucht sich das periodenfenster am
    #ende des signals selbst. wird der anlauf mitgemittelt, liegt St bis zu 2 % zu tief
    auswertung = werte_signal_aus(t, signal, cfg.D, cfg.U_inf)
    St = auswertung["St"]
    if np.isfinite(St):
        print(f"Strouhal-Zahl: St = {St:.4f}  (Amplitude {auswertung['amplitude']:.3f}, "
              f"{auswertung['n_perioden']} Perioden aus t = {auswertung['t_fenster_start']:.1f} "
              f"... {auswertung['t_fenster_ende']:.1f})")
    else:
        print(f"Keine periodische Ablösung auswertbar (Status {auswertung['status']}): "
              f"Lauf zu kurz, Re zu klein oder zu wenig Snapshots nach dem Einschwingen.")

    print("Wirbelstärke: Bild und Animation ...")
    wirbelstaerke_bild(os.path.join(ordner, "wirbelstaerke.png"), domain, cfg, t, omega,
                       signal, r_sonde, auswertung)
    wirbelstaerke_animation(os.path.join(ordner, "wirbelstaerke.gif"), domain, cfg, t, omega)

    print("Ljapunow-Exponent (FTLE) ...")
    try:
        felder, ausdehnung, T = ftle_felder(domain, cfg, t, ux, uy, St=St)
        ftle_ausgabe(os.path.join(ordner, "ftle.png"), os.path.join(ordner, "ftle.gif"),
                     cfg, felder, ausdehnung, T)
    except ValueError as fehler:
        print(f"  FTLE übersprungen: {fehler}")

    print(f"Ergebnisse gespeichert in {ordner}")
