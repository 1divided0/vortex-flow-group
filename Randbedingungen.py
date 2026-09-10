import numpy as np

def potential_flow_psi(domain, r_value):
    # ist psi = U_inf * sintheta(r-R^2/r) zim wiederverwenden 
    cfg = domain.cfg
    return cfg.U_inf * np.sin(domain.theta) * (r_value - cfg.R**2 / r_value)

def apply_wall_bc(psi, omega, domain):
    # BC am Zylinder mit no slip und no penetration
    cfg = domain.cfg
    i_w = domain.i_wall

    #kein durchfluss
    psi[i_w] = 0.0 

    #no slip
    omega[i_w, :] = -2.0 * (psi[i_w + 1] - psi[i_w, :]) / (cfg.R**2 * domain.dxi**2)
    return psi, omega

def apply_farfield_bc(psi, omega, domain):
    i_f = domain.i_far

    psi_potential = potential_flow_psi(domain, domain.r[i_f])

    inflow = domain.is_inflow
    outflow = domain.is_outflow

    #einströmseite:
    #maske sodass nur theta werte überschrieben werden (inflow)

    psi[i_f, inflow] = psi_potential[inflow]
    omega[i_f, inflow] = 0.0

    #ausstromseite: mit neumann bc (randwerte einfach wie vorheriger wert)
    psi[i_f, outflow] = psi[i_f - 1, outflow]
    omega[i_f, outflow] = omega[i_f - 1, outflow]

    return psi, omega

def apply_bc(psi, omega, domain):
    #bündelt die bc's

    psi, omega = apply_wall_bc(psi, omega, domain)
    psi, omega = apply_farfield_bc(psi, omega, domain)
    return psi, omega







