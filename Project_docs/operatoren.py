import numpy as np
import scipy.sparse as sp


def build_D_xi(domain):
    #erste ableitung in xi richtung
    n = domain.n_xi
    h = domain.dxi

    D = sp.diags(
        diagonals=[-1.0, 1.0],
        offsets=[-1, 1],
        shape=(n, n),
        format='lil'
    ) / (2.0 * h)

    #einseitige differenzen 2. ordnung am rand, damit die fehlerordnung
    #des gesamten operators bei 2 bleibt
    D[0, 0] = -3.0 / (2.0 * h)
    D[0, 1] = 4.0 / (2.0 * h)
    D[0, 2] = -1.0 / (2.0 * h)

    D[-1, -1] = 3.0 / (2.0 * h)
    D[-1, -2] = -4.0 / (2.0 * h)
    D[-1, -3] = 1.0 / (2.0 * h)

    return D.tocsr()


def build_D2_xi(domain):
    #zweite ableitung in xi richtung, zentrale differenzen 2. ordnung
    n = domain.n_xi
    h = domain.dxi

    D2 = sp.diags(
        diagonals=[1.0, -2.0, 1.0],
        offsets=[-1, 0, 1],
        shape=(n, n),
        format='lil'
    ) / (h**2)

    #hier keine identitaetszeilen fuer dirichlet eintragen: in build_laplacian
    #wuerde kron(I_xi, D2_theta) sie wieder verfaelschen und die metrik sie mit
    #1/r^2 skalieren. die randzeilen ersetzt Poisson.build_poisson_matrix komplett.
    #einseitige zweite ableitung 2. ordnung: (2f0 - 5f1 + 4f2 - f3)/h^2
    D2[0, 0] = 2.0 / h**2
    D2[0, 1] = -5.0 / h**2
    D2[0, 2] = 4.0 / h**2
    D2[0, 3] = -1.0 / h**2

    D2[-1, -1] = 2.0 / h**2
    D2[-1, -2] = -5.0 / h**2
    D2[-1, -3] = 4.0 / h**2
    D2[-1, -4] = -1.0 / h**2

    return D2.tocsr()


def build_D_theta(domain):
    #erste ableitung in theta richtung
    n = domain.n_theta
    h = domain.dtheta

    D = sp.diags(
        diagonals=[-1.0, 1.0],
        offsets=[-1, 1],
        shape=(n, n),
        format='lil'
    ) / (2.0 * h)

    #periodische randbedingungen
    D[0, -1] = -1.0 / (2.0 * h)
    D[-1, 0] = 1.0 / (2.0 * h)

    return D.tocsr()


def build_D2_theta(domain):
    #zweite ableitung in theta richtung
    n = domain.n_theta
    h = domain.dtheta

    D2 = sp.diags(
        diagonals=[1.0, -2.0, 1.0],
        offsets=[-1, 0, 1],
        shape=(n, n),
        format='lil'
    ) / (h**2)

    #periodische randbedingungen
    D2[0, -1] = 1.0 / (h**2)
    D2[-1, 0] = 1.0 / (h**2)

    return D2.tocsr()


def build_upwind_xi(domain):
    #vorwaerts- und rueckwaertsdifferenzen 1. ordnung in xi richtung
    n = domain.n_xi
    h = domain.dxi

    #in xi richtung gibt es keine periodizitaet: der vorwaertsdifferenz fehlt in
    #der letzten zeile der nachbar, der rueckwaertsdifferenz in der ersten.
    #ohne korrektur waere die zeilensumme dort ungleich null
    forward = sp.diags([-1.0, 1.0], [0, 1], shape=(n, n), format='lil') / h
    forward[-1, -1] = 1.0 / h
    forward[-1, -2] = -1.0 / h

    backward = sp.diags([-1.0, 1.0], [-1, 0], shape=(n, n), format='lil') / h
    backward[0, 0] = -1.0 / h
    backward[0, 1] = 1.0 / h

    return forward.tocsr(), backward.tocsr()


def build_upwind_theta(domain):
    #vorwaerts- und rueckwaertsdifferenzen 1. ordnung in theta richtung (periodisch)
    n = domain.n_theta
    h = domain.dtheta

    forward = sp.diags([-1.0, 1.0], [0, 1], shape=(n, n), format='lil') / h
    forward[-1, 0] = 1.0 / h

    backward = sp.diags([-1.0, 1.0], [-1, 0], shape=(n, n), format='lil') / h
    backward[0, -1] = -1.0 / h

    return forward.tocsr(), backward.tocsr()


def build_laplacian(domain):
    #2d laplace operator via kronecker: (1/r^2) * (d2/dxi2 + d2/dtheta2)
    n_xi = domain.n_xi
    n_theta = domain.n_theta

    D2_xi = build_D2_xi(domain)
    D2_theta = build_D2_theta(domain)

    I_xi = sp.identity(n_xi, format='csr')
    I_theta = sp.identity(n_theta, format='csr')

    L_uniform = sp.kron(D2_xi, I_theta, format='csr') + sp.kron(I_xi, D2_theta, format='csr')

    metric = np.repeat(domain.vorfaktor, n_theta)
    metric_diag = sp.diags(metric, format='csr')

    return metric_diag @ L_uniform
