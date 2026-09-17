import numpy as np
import scipy.sparse as sp

from operatoren import build_D_xi, build_D2_xi, build_D_theta, build_D2_theta
from operatoren import build_upwind_xi, build_upwind_theta, build_laplacian
from Randbedingungen import apply_bc

def build_alle_operatoren(domain):
    #sammelt alle operatoren in nem dictionary

    n_theta = domain.n_theta
    n_xi = domain.n_xi
    I_xi = sp.identity(n_xi, format='csr')
    I_theta = sp.identity(n_theta, format='csr')

    D_xi = build_D_xi(domain)
    D_theta = build_D_theta(domain)
    fwd_xi, bwd_xi = build_upwind_xi(domain)
    fwd_theta, bwd_theta = build_upwind_theta(domain)

    return {
        "D_xi": sp.kron(D_xi, I_theta, format='csr'),
        "D_theta": sp.kron(I_xi, D_theta, format='csr'),
        "fwd_xi": sp.kron(fwd_xi, I_theta, format='csr'),
        "bwd_xi": sp.kron(bwd_xi, I_theta, format='csr'),
        "fwd_theta": sp.kron(I_xi, fwd_theta, format='csr'),
        "bwd_theta": sp.kron(I_xi, bwd_theta, format='csr'),
        "L": build_laplacian(domain),
    }


def geschwindigkeit(psi, domain, ops): # ops muss dann noch definiert werden (sind die operatoren in der funktion zuvor)
    #berechnet aus psi die geschwindigkeiten 
    # u_r = 1/r * dpsi/dtheta
    # u_theta = -1/r * dpsi/dr
    #zentrale differenzen nicht upwind das kommt seperat weil andere methode und erst später benötigt für advektion und co

    dpsi_dtheta = domain.unflatten(ops["D_theta"] @ domain.flatten(psi))
    dpsi_dxi = domain.unflatten(ops["D_xi"] @ domain.flatten(psi))

    r = domain.r[:, None] #format (n_xi, 1) für broadcasting über alle theta spalten (broadcasting ist in numpy ne funktion welche arrays auf die benötigte größe anpasst)
    u_r = dpsi_dtheta / r
    u_theta = -dpsi_dxi / r

    return u_r, u_theta

def upwind_ableitung(field, velocity, forward_op, backward_op, domain): #auch hier müssen die paramter noch in main.py(später) definiert werden
    # hier wird punktweise ableitungen für upwind berechnet 
    # bei v >= 0 rückwerts differenz
    #bei v < 0 vorwärts differenz 

    d_forward = domain.unflatten(forward_op @ domain.flatten(field))
    d_backward = domain.unflatten(backward_op @ domain.flatten(field))
    return np.where(velocity >= 0.0, d_backward, d_forward)

def berechne_rhs(omega, psi, domain, cfg, ops):  #rhs sthet hier für die rechte seite der gleichung in dem fall für die rechte seite der poisson gleichung
    #berechnet domega/dt = -advektion + diffuision  an jedem gitterpunkt (innen)
    #nur von i bis i_far-1 notwendig rest ist BC 
    # advektion in xi theta ist d/dr = 1/r d/dxi mit (u.grad(omega)
    # u.grad(omega) = 1/r *(u_r * domega/dxi + u_theta * domega/dtheta)
    #Diffusion mit nu * Laplace(omega) 

    u_r, u_theta = geschwindigkeit(psi, domain, ops) #ops auch hier erst später definiert

    # zentrale differenzen statt upwind: upwind 1. ordnung wirkt wie eine
    # kuenstliche viskositaet ~ |u|*h/2. weil das log-gitter nach aussen groeber
    # wird, war die im nachlauf mehrfach so gross wie nu -> effektiv Re ~ 20,
    # und darunter loesen sich keine wirbel ab. rk4 ist mit zentralen
    # differenzen bei Re = 100 stabil. upwind_ableitung bleibt fuer ein
    # spaeteres hybridschema bei hoeheren Re erhalten.
    domega_dxi = domain.unflatten(ops["D_xi"] @ domain.flatten(omega))
    domega_dtheta = domain.unflatten(ops["D_theta"] @ domain.flatten(omega))

    r = domain.r[:, None]

    advektion = (u_r * domega_dxi + u_theta * domega_dtheta) / r

    diffusion = domain.unflatten(ops["L"] @ domain.flatten(omega)) * cfg.nu

    rhs = -advektion + diffusion

    return rhs

def cfl_zeitschritt(u_r, u_theta, domain, cfg):
    # advektiv:  dt <= cfl_target * dx_phys / |u|
    # diffusiv:  dt <= cfl_target * dx_phys^2 / (4*nu)
    # dx_phys ist physikalischer gitterabstand
    # in xi richtung: r*dxi (weil dr = r*dxi)
    # in theta richtung: r*dtheta 
    # eps verhindert divison durch 0, weil mathe sagt nein

    eps = 1e-10
    r = domain.r[:, None]

    dx_xi_phys = r * domain.dxi
    dx_dtheta_phys = r * domain.dtheta 

    dt_adv_xi = cfg.cfl_target * dx_xi_phys / (np.abs(u_r) + eps)
    dt_adv_theta = cfg.cfl_target * dx_dtheta_phys / (np.abs(u_theta) + eps)

    dx_min_phys = np.minimum(dx_xi_phys, dx_dtheta_phys)
    dt_diff = cfg.cfl_target * dx_min_phys**2 / (4.0 * cfg.nu)

    return float(np.min([dt_adv_xi.min(), dt_adv_theta.min(), dt_diff.min()]))

def rk4(psi, omega, cfg, poisson_solver, ops, dt, domain):



    def berechne(omega_stage):
        psi_stage = poisson_solver.löse(omega_stage)
        psi_stage, omega_stage = apply_bc(psi_stage, omega_stage, domain)
        k = berechne_rhs(omega_stage, psi_stage, domain, cfg, ops)
        return psi_stage, omega_stage, k

    _, omega0, k1 = berechne(omega)
    _, _, k2 = berechne(omega0 + 0.5 * dt * k1)
    _, _, k3 = berechne(omega0 + 0.5 * dt * k2)
    _, _, k4 = berechne(omega0 + dt * k3)

    omega_next = omega0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    psi_next = poisson_solver.löse(omega_next)
    psi_next, omega_next = apply_bc(psi_next, omega_next, domain)

    return psi_next, omega_next


