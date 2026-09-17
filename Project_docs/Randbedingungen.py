import numpy as np


def potential_flow_psi(domain, r_value):
    # psi = U_inf * sin(theta) * (r - R^2/r), zum wiederverwenden
    cfg = domain.cfg
    return cfg.U_inf * np.sin(domain.theta) * (r_value - cfg.R**2 / r_value)


def apply_wall_bc(psi, omega, domain):
    # BC am Zylinder mit no slip und no penetration
    cfg = domain.cfg
    i_w = domain.i_wall

    #kein durchfluss
    psi[i_w] = 0.0

    #no slip ueber die thom-formel: omega_w = -2 (psi_1 - psi_w) / (R^2 dxi^2)
    omega[i_w, :] = -2.0 * (psi[i_w + 1] - psi[i_w, :]) / (cfg.R**2 * domain.dxi**2)
    return psi, omega


def apply_farfield_bc(psi, omega, domain):
    i_f = domain.i_far

    psi_potential = potential_flow_psi(domain, domain.r[i_f])

    inflow = domain.is_inflow
    outflow = domain.is_outflow

    #psi: auf dem ganzen fernrand potentialströmung (passend zu Poisson.py),
    #denn psi legt den durchfluss fest
    psi[i_f, :] = psi_potential

    #omega: nur hier ein-/ausstrom unterscheiden. einströmseite ungestört,
    #ausströmseite nullgradient, damit wirbelstärke das gebiet verlassen kann
    omega[i_f, inflow] = 0.0
    omega[i_f, outflow] = omega[i_f - 1, outflow]

    return psi, omega


def apply_bc(psi, omega, domain):
    #bündelt die bc's. achtung: psi und omega werden direkt (in place) veraendert

    psi, omega = apply_wall_bc(psi, omega, domain)
    psi, omega = apply_farfield_bc(psi, omega, domain)
    return psi, omega
