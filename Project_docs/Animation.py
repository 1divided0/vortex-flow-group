import numpy as np
import matplotlib as plt
from matplotlib.animation import FuncAnimation, PillowWriter

daten = np.load("simulation_snapshots.npz")

t = daten["t"]
psi = daten["psi"]
omega = daten["omega"]



