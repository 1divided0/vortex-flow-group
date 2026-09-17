#grid

import numpy as np


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
