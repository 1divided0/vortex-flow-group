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
    d_backward = domain.unflatten(backward_op @ domain.unflatten(field))
    return np.where(velocity >= 0.0, d_backward, d_forward)

def berechne_rhs(omega, psi, domain, cfg, ops):  #rhs sthet hier für die rechte seite der gleichung in dem fall für die rechte seite der poisson gleichung
    #berechnet domega/dt = -advektion + diffuision  an jedem gitterpunkt (innen)
    #nur von i bis i_far-1 notwendig rest ist BC 
    # advektion in xi theta ist d/dr = 1/r d/dxi mit (u.grad(omega)
    # u.grad(omega) = 1/r *(u_r * domega/dxi + u_theta * domega/dtheta)
    #Diffusion mit nu * Laplace(omega) 

    u_r, u_theta = geschwindigkeit(psi, domain, ops) #ops auch hier erst später definiert 
    domega 