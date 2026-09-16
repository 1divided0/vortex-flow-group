import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from operatoren import build_laplacian
from Randbedingungen import potential_flow_psi


def build_poisson_matrix(domain):

    n_xi = domain.n_xi
    n_theta = domain.n_theta

    A = build_laplacian(domain).tolil()

    i_w = domain.i_wall
    i_f = domain.i_far 

    # wand mit dirichlet psi = 0 
    for j in range(n_theta):
        idx = i_w * n_theta + j
        A.rows[idx] = [idx]
        A.data[idx] = [1.0]

    # fernfeld
    for j in range(n_theta):
        idx = i_f * n_theta + j
        if domain.is_inflow[j]:
            # inflow: dirichlet wie bei der wand
            A.rows[idx] = [idx]
            A.data[idx] = [1.0]
        else:
            # Neumann nullgradient aus der Formel mit psi[i_far] - psi[i_far-1] = 0
            idx_next = (i_f - 1) * n_theta + j

            cols = sorted([idx, idx_next])
            vals = [1.0 if c == idx else -1.0 for c in cols]
            A.rows[idx] = cols
            A.data[idx] = vals

    return A.tocsr() #csr (compressed sparse row format ist wichtig für später)

def build_rhs(domain, omega):
    #rechte seite des GLS 

    n_xi = domain.n_xi
    n_theta = domain.n_theta
    i_w = domain.i_wall
    i_f = domain.i_far

    rhs = np.zeros((n_xi, n_theta))

    # innere Punkte: quellterm der poissongleichung 
    rhs[1:i_f, :] = -omega[1:i_f, :]

    # wand dirichlet zielwert = 0 -> rhs bleibt 0

    #fernfeld, einströmseite -> dirichlet zielwert = potentialströmung

    psi_potential = potential_flow_psi(domain, domain.r[i_f])
    rhs[i_f, domain.is_inflow] = psi_potential[domain.is_inflow]

    #fernfeld, ausstömseite -> neumann zielwert = 0 (rhs bleibt 0)

    return domain.flatten(rhs)

class PoissonSolver:
    # LR - zerlegung nach dem prinzip aus den Vorlesungen
    #als klasse weil wir es nur einmal brauchen und ne normale funktion zwischen verschiedenen aufrufen keine daten speichert/ beibehält

    def __init__(self, domain):
        self.domain = domain
        A = build_poisson_matrix(domain)
        self.lu = spla.splu(A)   #LR zerlegung für die matrix A

    def löse(self, omega):
        # löst die poissongleichung für ein gegebenes omega
        rhs = build_rhs(self.domain, omega)
        psi_flat = self.lu.solve(rhs)
        return self.domain.unflatten(psi_flat)

    

 


