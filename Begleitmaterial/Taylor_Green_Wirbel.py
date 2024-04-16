import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt


""" Visualisierung des Taylor-Green-Wirbels """

def Taylor_Green_Wirbel(
    nu:    float = 0,   # Kinematische Viskosität
    t_end: float = 5,   # Endzeit
    t_nr:  int   = 100, # Anzahl der Zeitschritte
    x_nr:  int   = 100, # Auflösung in x-Richtung
    y_nr:  int   = 20   # Auflösung in y-Richtung
):

    # Gitter
    x    = np.linspace(0,2*np.pi,x_nr)
    y    = np.linspace(0,2*np.pi,y_nr)
    dx   = x[1]-x[0]
    dy   = y[1]-y[0]
    X, Y = np.meshgrid(x,y)

    # Lösung
    U    = lambda t: np.sin(X)*np.cos(Y)*np.exp(-2*nu*t)
    V    = lambda t: -np.cos(X)*np.sin(Y)*np.exp(-2*nu*t)
    Psi  = lambda t: np.sin(X)*np.sin(Y)*np.exp(-2*nu*t)
    
    plt.ion()
    fig = plt.figure(figsize=(6,4))

    # Lagrange Partikel
    ax    = fig.add_subplot(1,2,1)
    X_pos = np.copy(X)
    Y_pos = np.copy(Y)
    sctr  = [ax.scatter(X_pos[row],Y_pos[row],s=2) for row in range(y_nr)]
    ax.set_xlim(0,2*np.pi)
    ax.set_ylim(0,2*np.pi)
    ax.set_title('Lagrange Partikel')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

    # Stromfunktion
    ax   = fig.add_subplot(1,2,2,projection='3d')
    surf = ax.plot_surface(X,Y,Psi(0),cmap=plt.cm.viridis)
    cont = ax.contour(X,Y,Psi(0),zdir='z',offset=-1,cmap=plt.cm.viridis)
    ax.set_xlim(0,2*np.pi)
    ax.set_ylim(0,2*np.pi)
    ax.set_zlim(-1,1)
    ax.set_title('Stromfunktion')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$\Psi$')

    # Zeitschleife
    time  = np.linspace(0,t_end,t_nr)
    t_stp = time[1]-time[0]
    for t in time:

        # Aktuelle Strömungsgrößen
        U_now   = U(t)
        V_now   = V(t)
        Psi_now = Psi(t)
        
        # Plot
        for row in range(y_nr):
            sctr[row].set_offsets(np.stack((X_pos[row],Y_pos[row])).T)
        surf.set_zdata = Psi_now
        cont.set_zdata = Psi_now
        fig.suptitle(fr'$\nu={nu},~t={t:.2f}$')
        fig.tight_layout(pad=2)
        fig.canvas.draw()
        fig.canvas.flush_events()
        
        # Nächster Zeitschritt
        X_pos += t_stp * RectBivariateSpline(x,y,U_now.T)(X_pos,Y_pos,grid=False)
        Y_pos += t_stp * RectBivariateSpline(x,y,V_now.T)(X_pos,Y_pos,grid=False)
