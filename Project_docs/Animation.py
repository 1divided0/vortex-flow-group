import numpy as np
import matplotlib as plt
import os 
import sys
from domain import Domain
from config import Config
from matplotlib.animation import FuncAnimation, PillowWriter
from log_polar_grid import erzeuge_log_polar_feld, rücktrafo
from zeitintegration import geschwindigkeit, build_alle_operatoren


def lade_snapshots(pfad):
    daten = np.load("pfad")

    t = daten["t"]
    psi = daten["psi"]
    omega = daten["omega"]

     #if "R" not in daten.files:
        #  raise ValueError(
       #     f"{pfad} enthält keine Parameter"
       # )

    cfg = Config(R = float(daten["R"]), r_max = float(daten["r_max"]), u_inf = float(daten["u_inf"]), Re = float(daten["Re"]), 
                 n_xi = int(daten["n_xi"]), n_theta = int(daten["n_theta"]), dt = float(daten["dt"]), cfl_target = float(daten["cfl_target"]))

    domain = Domain(cfg)
    t = daten["t"].astype(float)
    psi = daten["psi"].astype(float)
    omega = daten["omega"].astype(float)

    return cfg, domain, t, psi, omega

cfg, domain, t, psi, omega = lade_snapshots(
    "simulation_snapshots.npz"
)

r = domain.r
theta = domain.theta

X, Y = rücktrafo(r, theta)
    









