#gitter: konfiguration, log-polares koordinatensystem und das rechengebiet.
#fasst die frueheren dateien config.py, log_polar_grid.py und domain.py zusammen

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# konfiguration
# ---------------------------------------------------------------------------

@dataclass
class Config:
    R: float           #zylinderradius
    r_max: float       #äußerer rand

    U_inf: float       #anströmgeschwindigkeit
    Re: float          #reynoldszahl (mit D = 2R)

    #auflösung
    n_xi: int          #anzahl der gitterpunkte in radialer richtung (xi-richtung)
    n_theta: int       #anzahl der gitterpunkte in umfangsrichtung (theta-richtung)

    dt: float          #maximale zeitschrittweite - den tatsaechlichen schritt bestimmt die cfl-bedingung,
                       #ein zu kleines dt hier bremst die rechnung unnoetig aus

    cfl_target: float = 0.5   #cfl wert (numerischer sicherheitswert)

    @property
    def D(self) -> float:    #zylinderdurchmesser für Re
        return 2.0 * self.R

    @property
    def nu(self) -> float:   #kinematische viskosität aus Re = U_inf * D / nu
        return self.U_inf * self.D / self.Re


# ---------------------------------------------------------------------------
# log-polares gitter:  xi = ln(r/R), theta periodisch
# ---------------------------------------------------------------------------

def erzeuge_log_polar_feld(R, r_max, n_xi, n_theta):
    #xi = ln(r/R) aequidistant von der wand (xi = 0) bis r_max,
    #theta periodisch ohne den doppelten punkt 2pi
    xi_max = np.log(r_max / R)
    xi = np.linspace(0.0, xi_max, n_xi)
    dxi = xi[1] - xi[0]

    theta = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    dtheta = theta[1] - theta[0]

    r = R * np.exp(xi)

    return xi, theta, r, dxi, dtheta


def rücktrafo(r, theta):
    #(r, theta) -> kartesische koordinaten, form (n_xi, n_theta)
    R_grid, Theta_grid = np.meshgrid(r, theta, indexing='ij')
    X = R_grid * np.cos(Theta_grid)
    Y = R_grid * np.sin(Theta_grid)
    return X, Y


def vorfaktor(r):
    #metrischer vorfaktor des laplace-operators in log-polarkoordinaten
    return 1.0 / r**2


# ---------------------------------------------------------------------------
# rechengebiet
# ---------------------------------------------------------------------------

class Domain:

    def __init__(self, config):
        # config wird als Ganzes gespeichert (nicht nur einzelne Werte
        # rauskopiert), damit später jede Funktion, die ein Domain-Objekt
        # bekommt, auch automatisch Zugriff auf z.B. cfg.Re, cfg.dt hat
        self.cfg = config

        #trennung von gittererzeugung und gitter arbeit
        self.xi, self.theta, self.r, self.dxi, self.dtheta = erzeuge_log_polar_feld(
            config.R, config.r_max, config.n_xi, config.n_theta
        )
        self.n_xi = config.n_xi
        self.n_theta = config.n_theta
        #metrischer vorfaktor wird einmal berechnet und gespeichert statt immer wieder neu
        self.vorfaktor = vorfaktor(self.r)

        #rand indizes
        self.i_wall = 0
        self.i_far = self.n_xi - 1

        #in und outflow
        #die toleranz ist noetig: cos(90°) ergibt +6e-17, cos(270°) aber -1.8e-16.
        #ohne toleranz waere oben ausstrom und unten einstrom, also bekaemen die
        #beiden spiegelpunkte durch reinen rundungszufall verschiedene randbedingungen
        self.is_inflow = np.cos(self.theta) < -1e-12
        self.is_outflow = ~self.is_inflow

    def flatten(self, field2d):
        #2d feld (n_xi, n_theta) -> 1d vektor, zeilenweise: k = i*n_theta + j
        return field2d.reshape(-1)

    def unflatten(self, field1d):
        return field1d.reshape(self.n_xi, self.n_theta)

    def kartesisch(self):
        return rücktrafo(self.r, self.theta)
