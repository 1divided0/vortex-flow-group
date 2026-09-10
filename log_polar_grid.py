#grid

import numpy as np 

#gridfunktion mit paramteren aus der config datei



def erzeuge_log_polar_feld(R, r_max, n_xi, n_theta):
    xi_max = np.log(r_max / R)
    xi = np.linspace(0.0, xi_max, n_xi)
    dxi = xi[1] - xi[0]

    theta = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    dtheta = theta[1] - theta[0]

    r = R * np.exp(xi)

    return xi, theta, r, dxi, dtheta

def rücktrafo(r, theta):
    R_grid, Theta_grid = np.meshgrid(r, theta, indexing='ij')
    X = R_grid * np.cos(Theta_grid)
    Y = R_grid * np.sin(Theta_grid)
    return X, Y

def vorfaktor(r):
    return 1.0 / r**2 

