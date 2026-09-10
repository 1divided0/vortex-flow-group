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

    #für die randzeilen kommen erstmal nur platzhalter mit einseitigen differnzen 1. ordnung 
    D[0, 0] = -1.0 / h
    D[0, 1] = 1.0 / h
    D[-1, -1] = 1.0 / h
    D[-1, -2] = -1.0 / h

    return D.tocsr()

def build_D2_xi(domain):
    #wieder mit zentralen differenzen 2. ordnung

    n = domain.n_xi
    h = domain.dxi

    D2 = sp.diags(
        diagonals=[1.0, -2.0, 1.0],
        offsets=[-1, 0, 1],
        shape=(n, n),
        format='lil'
    ) / (h**2)

    #wieder randzeilen für später ersetzen
    D2[0, :] = 0.0
    D2[0, 0] = 1.0
    D2[-1, :] = 0.0
    D2[-1, -1] = 1.0 
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

    n = domain.n_xi
    h = domain.dxi

    foreward = sp.diags([-1.0, 1.0], [0, 1], shape=(n, n), format='lil') / h
    backward = sp.diags([-1.0, 1.0], [-1, 0], shape=(n, n), format='lil') / h

    return foreward.tocsr(), backward.tocsr()

def build_upwind_theta(domain):

    n = domain.n_theta
    h = domain.dtheta 

    foreward = sp.diags([-1.0, 1.0], [0, 1], shape=(n, n), format='lil') / h
    foreward[-1, 0] = 1.0 / h

    backward = sp.diags([-1.0, 1.0], [-1, 0], shape=(n, n), format='lil') / h
    backward[0, -1] = -1.0 / h

    return foreward.tocsr(), backward.tocsr()       

def build_laplacian(domain):
    #2d laplace operator via kronecker

    n_xi = domain.n_xi
    n_theta = domain.n_theta

    D2_xi = build_D2_xi(domain)
    D2_theta = build_D2_theta(domain)

    I_xi = sp.identity(n_xi, format='csr')
    I_theta = sp.identity(n_theta, format='csr')

    L_uniform = sp.kron(D2_xi, I_theta, format = 'csr') + sp.kron(I_xi, D2_theta, format = 'csr')

    metric = np.repeat(domain.vorfaktor, n_theta)
    metric_diag = sp.diags(metric, format='csr')

    return metric_diag @ L_uniform






#test bereich 
if __name__ == "__main__":
    from config import Config
    from domain import Domain
 
    cfg = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=80, n_theta=160, dt=1e-3)
    dom = Domain(cfg)
 
    # --- Test 1: D_xi, D2_xi an einer Funktion mit bekannter Ableitung ---
    # f(xi) = xi^3, unabhaengig von theta -- Ableitungen analytisch bekannt:
    # f' = 3*xi^2, f'' = 6*xi
    f = np.tile(dom.xi[:, None]**3, (1, dom.n_theta))
 
    # D_xi/D2_xi wirken nur auf die xi-Richtung (Groesse n_xi) -- um sie
    # auf ein volles, geflattetes 2D-Feld (Groesse n_xi*n_theta)
    # anzuwenden, muessen sie erst per Kronecker-Produkt mit der
    # Identitaet in theta-Richtung auf die volle Groesse gebracht werden
    # (exakt dasselbe Prinzip wie beim Laplace-Aufbau weiter unten)
    I_theta = sp.identity(dom.n_theta, format="csr")
    D_xi_full = sp.kron(build_D_xi(dom), I_theta, format="csr")
    D2_xi_full = sp.kron(build_D2_xi(dom), I_theta, format="csr")
 
    df_numeric = dom.unflatten(D_xi_full @ dom.flatten(f))
    df_analytic = 3.0 * dom.xi[:, None]**2
 
    # nur INNERE Punkte vergleichen (Randzeilen sind nur Platzhalter,
    # siehe Docstrings oben)
    interior = slice(1, -1)
    err_D_xi = np.max(np.abs(df_numeric[interior, :] - df_analytic[interior, :]))
    print(f"max Fehler D_xi (innere Punkte, f=xi^3):  {err_D_xi:.6f}")
 
    d2f_numeric = dom.unflatten(D2_xi_full @ dom.flatten(f))
    d2f_analytic = 6.0 * dom.xi[:, None]
    err_D2_xi = np.max(np.abs(d2f_numeric[interior, :] - d2f_analytic[interior, :]))
    print(f"max Fehler D2_xi (innere Punkte, f=xi^3):  {err_D2_xi:.6f}")
 
    # --- Test 2: D_theta, D2_theta an sin(theta), periodisch ---
    g = np.tile(np.sin(dom.theta)[None, :], (dom.n_xi, 1))
 
    I_xi = sp.identity(dom.n_xi, format="csr")
    D_theta_full = sp.kron(I_xi, build_D_theta(dom), format="csr")
    D2_theta_full = sp.kron(I_xi, build_D2_theta(dom), format="csr")
 
    dg_numeric = dom.unflatten(D_theta_full @ dom.flatten(g))
    dg_analytic = np.cos(dom.theta)[None, :]
    err_D_theta = np.max(np.abs(dg_numeric - dg_analytic))
    print(f"max Fehler D_theta (periodisch, g=sin(theta)):  {err_D_theta:.6f}")
 
    d2g_numeric = dom.unflatten(D2_theta_full @ dom.flatten(g))
    d2g_analytic = -np.sin(dom.theta)[None, :]
    err_D2_theta = np.max(np.abs(d2g_numeric - d2g_analytic))
    print(f"max Fehler D2_theta (periodisch, g=sin(theta)):  {err_D2_theta:.6f}")
 
    # --- Test 3: Laplace-Operator an einer Funktion mit bekanntem Laplace ---
    # h(xi,theta) = xi^2 * sin(theta)
    # nabla^2 h = (1/r^2) * (h_xixi + h_thetatheta)
    #           = (1/r^2) * (2*sin(theta) - xi^2*sin(theta))
    h_field = (dom.xi[:, None]**2) * np.sin(dom.theta)[None, :]
 
    L = build_laplacian(dom)
    lap_numeric = dom.unflatten(L @ dom.flatten(h_field))
    lap_analytic = dom.vorfaktor[:, None] * (2.0 * np.sin(dom.theta)[None, :] - (dom.xi[:, None]**2) * np.sin(dom.theta)[None, :])
 
    err_L = np.max(np.abs(lap_numeric[interior, :] - lap_analytic[interior, :]))
    print(f"max Fehler Laplace-Operator (innere Punkte):  {err_L:.6f}")
 
