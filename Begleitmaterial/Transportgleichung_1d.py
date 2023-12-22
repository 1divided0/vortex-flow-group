import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as ani


""" Hauptprogramm: Transportgleichung_1d """

def Transportgleichung_1d(
    u: float, # Transportgeschwindigkeit
    c: float, # Diffusionskonstante
    r: float, # Quelle (proportional zur Konzentration)
    T: float, # Simulationsendzeit
    gif_export: bool = False # optionaler gif Export
):
    # Diskretisierung des eindimensionalen Rechenraumes
    (x, lbda) = Raumdiskretisierung()

    # Startverteilung
    phi0 = Startverteilung(x)

    # Definition der in der Zeit zu integrierenden Funktion
    # -> Zeitfunktion(phi,u,c,r,lbda)

    # Lösung mit Zeitschrittverfahren
    sol = sp.integrate.solve_ivp(Zeitfunktion,(0,T),phi0,args=(u,c,r,lbda))
    t   = sol.t
    phi = sol.y

    # Plot der zeitlichen Entwicklung
    Plot(u,c,r,T,x,t,phi,gif_export=gif_export)


""" Subfunktion: Raumdiskretisierung """

def Raumdiskretisierung( ):
# Rückgabe
# --------
# x:    äquidistantes 1d-Gitter
# lbda: diskrete Wellenzahlen

    N    = 256 # Anzahl der Stützstellen
    L    = 30  # Gesamtlänge

    x    = np.linspace(0,L,N)
    lbda = (2*np.pi/L)*np.arange(-N/2,N/2)
    lbda = np.fft.fftshift(lbda) # für die FFT benötigte Anordnung
    
    return x, lbda


""" Subfunktion: Startverteilung """

def Startverteilung( x ):
# Rückgabe
# --------
# phi0: Startverteilung der Strömungsgröße als Hut-Funktion

    N = len(x)   # Anzahl der Stützstellen
    a = 1        # linke Position
    b = 2        # rechte Position
    s = 0.8      # Skalierungsfaktor
    phimax = 0.1 # Maximalwert

    # linke Seite
    phi1 = (x-a*np.ones(N)+s*np.ones(N))**3 * (a*np.ones(N)-x+s*np.ones(N))**3
    phi1[np.logical_or(x<a-s,x>a+s)] = 0

    # rechte Seite
    phi2 = (x-b*np.ones(N)+s*np.ones(N))**3 * (b*np.ones(N)-x+s*np.ones(N))**3
    phi2[np.logical_or(x<b-s,x>b+s)] = 0

    # Superposition und Normierung
    phi0 = (phi1+phi2)/np.max(phi1+phi2)*phimax
    
    return phi0


""" Subfunktion: Zeitfunktion """

def Zeitfunktion( t, phi, u, c, r, lbda ):
# Rückgabe
# --------
# phi_dot: Zeitableitung der Strömungsgröße

    return np.real(
            np.fft.ifft((c*(1j*lbda)**2-u*(1j*lbda))*np.fft.fft(phi))
        ) + r*phi


""" Subfunktion: Plot """

def Plot(
    u: float,        # Transportgeschwindigkeit
    c: float,        # Diffusionskonstante
    r: float,        # Quelle (proportional zur Konzentration)
    T: float,        # Simulationsendzeit
    x: np.ndarray,   # äquidistantes 1d-Gitter
    t: np.ndarray,   # Zeitwerte
    phi: np.ndarray, # Funktionswerte
    gif_export: bool = False # optionaler gif Export
):

    frms = 100 # Anzahl der Bilder
    rate = 50  # Bildwiederholungsrate in Millisekunden
    
    # Ploteinstellungen
    fig  = plt.figure(figsize=(4,2))
    line = plt.plot([])[0]
    plt.ylabel(r'$\Phi$')
    plt.xlabel(r'$x$')
    plt.xlim(0,x[-1])
    plt.ylim(-0.01,.12)
    plt.grid(True)
    
    # Animation
    update = lambda i: (
        line.set_data(x,phi[:,i]),
        plt.title(fr'$u={u}$, $c={c}$, $r={r}$, $t={t[i]:.2f}$'),
        fig.tight_layout(pad=0.2)
    )
    seq = ani.FuncAnimation(
        fig,
        lambda frm: update(np.linspace(0,len(t)-1,frms,dtype=int)[frm]),
        frames=frms,
        interval=rate
    )
    
    # optionaler gif Export
    if gif_export:
        seq.save(f'u{u}_c{c}_r{r}_T{T}.gif')
    else:
        plt.show()
