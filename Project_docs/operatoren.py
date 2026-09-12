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
    #des gesamten operators bei 2 bleibt (an der wand wird daraus u_theta berechnet)
    D[0, 0] = -3.0 / (2.0 * h)
    D[0, 1] = 4.0 / (2.0 * h)
    D[0, 2] = -1.0 / (2.0 * h)

    D[-1, -1] = 3.0 / (2.0 * h)
    D[-1, -2] = -4.0 / (2.0 * h)
    D[-1, -3] = 1.0 / (2.0 * h)

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

    #WICHTIG: hier duerfen keine randbedingungen stehen. frueher stand hier eine
    #identitaetszeile fuer dirichlet - die wurde aber in build_laplacian durch
    #kron(I_xi, D2_theta) wieder ueberschrieben und zusaetzlich mit 1/r^2 skaliert.
    #die randbedingungen setzt jetzt build_poisson_matrix ueber boolesche masken.
    #hier stehen einseitige differenzen 2. ordnung: (2f0 - 5f1 + 4f2 - f3)/h^2
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

    n = domain.n_xi
    h = domain.dxi

    #in xi richtung gibt es keine periodizitaet, also fehlt der vorwaertsdifferenz
    #in der letzten zeile der nachbar (und der rueckwaertsdifferenz in der ersten).
    #ohne korrektur waere die zeilensumme dort ungleich null, d.h. ein konstantes
    #feld haette eine ableitung ungleich null -> kuenstliche konvektion am rand
    forward = sp.diags([-1.0, 1.0], [0, 1], shape=(n, n), format='lil') / h
    forward[-1, -1] = 1.0 / h
    forward[-1, -2] = -1.0 / h

    backward = sp.diags([-1.0, 1.0], [-1, 0], shape=(n, n), format='lil') / h
    backward[0, 0] = -1.0 / h
    backward[0, 1] = 1.0 / h

    return forward.tocsr(), backward.tocsr()

def build_upwind_theta(domain):

    n = domain.n_theta
    h = domain.dtheta

    forward = sp.diags([-1.0, 1.0], [0, 1], shape=(n, n), format='lil') / h
    forward[-1, 0] = 1.0 / h

    backward = sp.diags([-1.0, 1.0], [-1, 0], shape=(n, n), format='lil') / h
    backward[0, -1] = -1.0 / h

    return forward.tocsr(), backward.tocsr()

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


def build_randmasken(domain):
    #boolesche masken der randtypen, schon geflattet.
    #psi wird am gesamten fernfeld per dirichlet gesetzt: psi kodiert den
    #volumenstrom, laesst man es am ausstrom frei schwimmen, ist die
    #durchstroemung nicht mehr festgelegt. die unterscheidung ein-/ausstrom
    #gilt nur fuer omega (siehe Randbedingungen.apply_farfield_bc)
    n_xi, n_theta = domain.n_xi, domain.n_theta

    dirichlet = np.zeros((n_xi, n_theta), dtype=bool)
    dirichlet[domain.i_wall, :] = True      #wand: psi = 0
    dirichlet[domain.i_far, :] = True       #fernfeld: potentialstroemung

    #vorgehalten fuer eine spaetere konvektive ausstrombedingung
    neumann = np.zeros((n_xi, n_theta), dtype=bool)

    return domain.flatten(dirichlet), domain.flatten(neumann)


def build_poisson_matrix(domain, L=None):
    """systemmatrix fuer nabla^2 psi = -omega inklusive randbedingungen.

    innenzeilen: laplace-operator
    dirichlet:   psi_k = psi_rand_k
    neumann:     psi_k - psi_(k-1 in xi) = 0   (nullgradient nach aussen)

    rueckgabe: (A, dirichlet_maske, neumann_maske)
    """
    if L is None:
        L = build_laplacian(domain)

    dirichlet, neumann = build_randmasken(domain)
    rand = dirichlet | neumann

    A = (sp.diags((~rand).astype(float)) @ L
         + sp.diags(dirichlet.astype(float))).tolil()

    #ein schritt in xi-richtung ueberspringt eine ganze theta-zeile,
    #weil flatten zeilenweise ueber (n_xi, n_theta) laeuft
    for k in np.flatnonzero(neumann):
        A[k, k] = 1.0
        A[k, k - domain.n_theta] = -1.0

    return A.tocsr(), dirichlet, neumann






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
 
    # ALLE Punkte inklusive Rand vergleichen. frueher wurden die randzeilen
    # per slice(1,-1) ausgeklammert - genau dort sassen aber die fehler
    err_D_xi = np.max(np.abs(df_numeric - df_analytic))
    print(f"max Fehler D_xi (inkl. Rand, f=xi^3):  {err_D_xi:.6f}")

    d2f_numeric = dom.unflatten(D2_xi_full @ dom.flatten(f))
    d2f_analytic = 6.0 * dom.xi[:, None]
    err_D2_xi = np.max(np.abs(d2f_numeric - d2f_analytic))
    print(f"max Fehler D2_xi (inkl. Rand, f=xi^3):  {err_D2_xi:.6f}")
 
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
 
    err_L = np.max(np.abs(lap_numeric - lap_analytic))
    print(f"max Fehler Laplace-Operator (inkl. Rand):  {err_L:.6f}")

    # --- Test 4: Zeilensummen der Upwind-Matrizen ---
    # eine ableitungsmatrix muss ein konstantes feld auf null abbilden,
    # sonst entsteht am rand kuenstliche konvektion
    for name, (fw, bw) in (("xi", build_upwind_xi(dom)), ("theta", build_upwind_theta(dom))):
        s = max(np.abs(np.asarray(fw.sum(axis=1))).max(),
                np.abs(np.asarray(bw.sum(axis=1))).max())
        print(f"max |Zeilensumme| Upwind {name} (soll 0):  {s:.2e}")

    # --- Test 5: Abnahmetest der Poisson-Matrix ---
    # bei omega = 0 muss exakt die Potentialstroemung herauskommen. das prueft
    # gitter, metrik, operator und randbedingungen in einem durchgang
    import scipy.sparse.linalg as spla
    from Randbedingungen import potential_flow_psi

    A, dirichlet, neumann = build_poisson_matrix(dom)

    psi_rand = np.zeros((dom.n_xi, dom.n_theta))
    psi_rand[dom.i_far] = potential_flow_psi(dom, dom.r[dom.i_far])
    psi_rand[dom.i_wall] = 0.0

    rhs = np.where(dirichlet, dom.flatten(psi_rand), 0.0)
    rhs[neumann] = 0.0
    psi = dom.unflatten(spla.spsolve(A.tocsc(), rhs))

    R_grid, T_grid = np.meshgrid(dom.r, dom.theta, indexing="ij")
    psi_exakt = cfg.U_inf * np.sin(T_grid) * (R_grid - cfg.R**2 / R_grid)

    print(f"psi an der Wand (soll 0):  {np.max(np.abs(psi[dom.i_wall])):.2e}")
    print(f"max Fehler Poisson gegen Potentialstroemung:  {np.max(np.abs(psi - psi_exakt)):.6f}")
