#loeser: alles, was pro zeitschritt gebraucht wird - randbedingungen, poisson-loeser
#und zeitintegration. fasst die frueheren dateien Randbedingungen.py, Poisson.py und
#zeitintegration.py zusammen. die reinen differenzenmatrizen stehen in operatoren.py

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from operatoren import build_D_xi, build_D_theta
from operatoren import build_upwind_xi, build_upwind_theta, build_laplacian
from operatoren import build_D2_xi, build_D2_theta


# ---------------------------------------------------------------------------
# randbedingungen
# ---------------------------------------------------------------------------

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

    #psi: auf dem ganzen fernrand potentialströmung (passend zu build_poisson_matrix),
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


# ---------------------------------------------------------------------------
# poisson-gleichung  nabla^2 psi = -omega
# ---------------------------------------------------------------------------

def build_poisson_matrix(domain):
    #laplace-matrix, deren randzeilen durch dirichlet-zeilen (identitaet) ersetzt werden

    n_theta = domain.n_theta
    A = build_laplacian(domain).tolil()

    i_w = domain.i_wall
    i_f = domain.i_far

    # wand mit dirichlet psi = 0
    for j in range(n_theta):
        idx = i_w * n_theta + j
        A.rows[idx] = [idx]
        A.data[idx] = [1.0]

    # fernfeld: auf dem GANZEN rand dirichlet mit der potentialstroemung.
    # psi gibt den volumenstrom vor - mit nullgradient am ausstrom waere der
    # durchfluss nicht festgelegt (bei omega = 0 kam an der wand u_theta = 1.55
    # statt 2.0 heraus). ein-/ausstrom wird nur fuer omega unterschieden.
    for j in range(n_theta):
        idx = i_f * n_theta + j
        A.rows[idx] = [idx]
        A.data[idx] = [1.0]

    return A.tocsr()


def build_rhs(domain, omega):
    #rechte seite des GLS

    i_f = domain.i_far
    rhs = np.zeros((domain.n_xi, domain.n_theta))

    # innere Punkte: quellterm der poissongleichung
    rhs[1:i_f, :] = -omega[1:i_f, :]

    # wand dirichlet zielwert = 0 -> rhs bleibt 0

    #fernfeld (ganzer rand) -> dirichlet zielwert = potentialströmung
    rhs[i_f, :] = potential_flow_psi(domain, domain.r[i_f])

    return domain.flatten(rhs)


class PoissonSolver:
    # LR-zerlegung nach dem prinzip aus den vorlesungen
    # als klasse, weil die zerlegung nur einmal im konstruktor berechnet und dann
    # fuer jeden aufruf von löse wiederverwendet wird (eine normale funktion
    # koennte sie zwischen den aufrufen nicht behalten)

    def __init__(self, domain):
        self.domain = domain
        A = build_poisson_matrix(domain)
        self.lu = spla.splu(A.tocsc())   #LR zerlegung für die matrix A (splu erwartet csc)

    def löse(self, omega):
        # löst die poissongleichung für ein gegebenes omega.
        # es gehen nur die inneren zeilen von omega ein, die randzeilen nicht
        rhs = build_rhs(self.domain, omega)
        psi_flat = self.lu.solve(rhs)
        return self.domain.unflatten(psi_flat)

class PoissonSolverFFT:
    """
    poisson-loeser ueber eine FFT in theta. loest exakt dasselbe gleichungssystem wie
    PoissonSolver, nur anders zerlegt:
    theta ist periodisch und aequidistant, damit ist D2_theta eine zyklische matrix und
    wird von der diskreten fourier-transformation diagonalisiert. die metrik 1/r^2 haengt
    nur von xi ab, und die dirichlet-randzeilen sind in theta diagonal - deshalb zerfaellt
    das 2d-system in n_theta/2+1 unabhaengige 1d-systeme in xi, eines je fourier-mode.
    das gilt nur fuer ein gleichmaessiges, periodisches theta-gitter; wird daran etwas
    geaendert, ist nur noch PoissonSolver gueltig.
    """

    def __init__(self, domain):
        self.domain = domain
        i_w = domain.i_wall
        i_f = domain.i_far

        #eigenwerte von D2_theta: bei einer zyklischen matrix sind das die
        #fourier-koeffizienten der ersten zeile. rfft liefert nur die haelfte der
        #moden (reelles feld), entsprechend werden auch nur diese gebraucht
        D2_theta = build_D2_theta(domain)
        erste_zeile = np.asarray(D2_theta[0, :].todense()).ravel()
        self.eigenwerte = np.fft.fft(erste_zeile).real[: domain.n_theta // 2 + 1]

        D2_xi = build_D2_xi(domain)
        I_xi = sp.identity(domain.n_xi, format="csr")
        metrisch_diag = sp.diags(domain.vorfaktor, format="csr")

        #je mode eine kleine (n_xi x n_xi)-matrix, mit denselben dirichlet-randzeilen
        #wie in build_poisson_matrix
        self.löser_pro_mode = []
        for lam in self.eigenwerte:
            A_m = (metrisch_diag @ (D2_xi + lam * I_xi)).tolil()

            A_m.rows[i_w] = [i_w]
            A_m.data[i_w] = [1.0]
            A_m.rows[i_f] = [i_f]
            A_m.data[i_f] = [1.0]

            self.löser_pro_mode.append(spla.splu(A_m.tocsc()))

    def löse(self, omega):
        #dieselbe rechte seite wie beim LR-loeser, nur modenweise geloest
        rhs = self.domain.unflatten(build_rhs(self.domain, omega))
        rhs_dach = np.fft.rfft(rhs, axis=1)

        psi_dach = np.empty_like(rhs_dach)
        for m, lu in enumerate(self.löser_pro_mode):
            #die matrizen sind reell: real- und imaginaerteil zusammen in einem
            #aufruf loesen spart gegenueber zwei aufrufen rund 20 % (ergebnis identisch)
            x = lu.solve(np.column_stack([rhs_dach[:, m].real, rhs_dach[:, m].imag]))
            psi_dach[:, m] = x[:, 0] + 1j * x[:, 1]

        return np.fft.irfft(psi_dach, n=self.domain.n_theta, axis=1)


#ab dieser anzahl unbekannter ist der FFT-loeser schneller. gemessen auf einem
#Ryzen 7 9800X3D: bei 81x160 = 12960 unbekannten ist LR noch 20 % schneller
#(der FFT-loeser ruft pro loesung n_theta/2 mal SuperLU auf), bei 121x240 = 29040
#ist FFT 30 % schneller, bei 161x320 = 51520 fast doppelt so schnell
FFT_AB_UNBEKANNTEN = 20_000


def wähle_poisson_löser(domain, schwelle=FFT_AB_UNBEKANNTEN):
    #beide loeser loesen dasselbe gleichungssystem (abweichung ~1e-12), die wahl
    #aendert die ergebnisse also nicht, nur die rechenzeit
    if domain.n_xi * domain.n_theta > schwelle:
        return PoissonSolverFFT(domain)
    return PoissonSolver(domain)


# ---------------------------------------------------------------------------
# zeitintegration
# ---------------------------------------------------------------------------

def build_alle_operatoren(domain):
    #sammelt alle operatoren (auf volle feldgroesse n_xi*n_theta gebracht) in nem dictionary

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


def geschwindigkeit(psi, domain, ops):
    #berechnet aus psi die geschwindigkeiten (zentrale differenzen)
    # u_r     =  1/r * dpsi/dtheta
    # u_theta = -dpsi/dr = -1/r * dpsi/dxi

    dpsi_dtheta = domain.unflatten(ops["D_theta"] @ domain.flatten(psi))
    dpsi_dxi = domain.unflatten(ops["D_xi"] @ domain.flatten(psi))

    r = domain.r[:, None] #format (n_xi, 1) für broadcasting über alle theta spalten (broadcasting ist in numpy ne funktion welche arrays auf die benötigte größe anpasst)
    u_r = dpsi_dtheta / r
    u_theta = -dpsi_dxi / r

    return u_r, u_theta


def upwind_ableitung(field, velocity, forward_op, backward_op, domain):
    #punktweise upwind-ableitung (derzeit ungenutzt, siehe berechne_rhs)
    # bei v >= 0 rueckwaertsdifferenz
    # bei v < 0 vorwaertsdifferenz

    d_forward = domain.unflatten(forward_op @ domain.flatten(field))
    d_backward = domain.unflatten(backward_op @ domain.flatten(field))
    return np.where(velocity >= 0.0, d_backward, d_forward)


def berechne_rhs(omega, psi, domain, cfg, ops):
    #rechte seite der transportgleichung: domega/dt = -advektion + diffusion
    #die werte in den randzeilen sind bedeutungslos, dort setzt apply_bc omega
    # advektion in xi theta mit d/dr = 1/r d/dxi:
    # u.grad(omega) = 1/r * (u_r * domega/dxi + u_theta * domega/dtheta)
    # diffusion mit nu * laplace(omega)

    u_r, u_theta = geschwindigkeit(psi, domain, ops)

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


def rk4(psi, omega, cfg, poisson_solver, ops, dt, domain, randbedingungen=apply_bc):
    # klassisches runge-kutta 4. ordnung fuer omega, psi folgt in jeder stufe aus poisson.
    # randbedingungen ist die funktion, die nach jeder poisson-loesung die raender setzt.
    # voreinstellung ist apply_bc (zylinder); ein periodisches gebiet braucht keine und
    # uebergibt eine funktion, die psi und omega unveraendert zurueckgibt
    # voraussetzung: psi passt zu omega (psi aus poisson, randbedingungen gesetzt).
    # das gilt fuer den anfangszustand aus main.Anfangsbedingungen und fuer jedes
    # ergebnis von rk4. deshalb braucht stufe 1 keine eigene poisson-loesung -
    # sie wuerde bitgenau dasselbe psi liefern und kostete 1 von 5 loesungen pro schritt

    def berechne(omega_stage):
        psi_stage = poisson_solver.löse(omega_stage)
        psi_stage, omega_stage = randbedingungen(psi_stage, omega_stage, domain)
        return berechne_rhs(omega_stage, psi_stage, domain, cfg, ops)

    k1 = berechne_rhs(omega, psi, domain, cfg, ops)
    k2 = berechne(omega + 0.5 * dt * k1)
    k3 = berechne(omega + 0.5 * dt * k2)
    k4 = berechne(omega + dt * k3)

    omega_next = omega + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    psi_next = poisson_solver.löse(omega_next)
    psi_next, omega_next = randbedingungen(psi_next, omega_next, domain)

    return psi_next, omega_next
