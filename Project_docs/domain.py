import numpy as np 
from log_polar_grid import erzeuge_log_polar_feld, rücktrafo, vorfaktor

class Domain:

    def __init__(self, config): 
        # config wird als Ganzes gespeichert (nicht nur einzelne Werte
        # rauskopiert), damit später jede Funktion, die ein Domain-Objekt
        # bekommt, auch automatisch Zugriff auf z.B. cfg.Re, cfg.dt hat
        self.cfg = config                     

        #trennung von gittererzeugung und gitter arbeit
        self.xi, self.theta, self.r, self.dxi, self.dtheta = erzeuge_log_polar_feld(
            config.R, config.r_max, config.n_xi, config.n_theta
        )
        self.n_xi = config.n_xi
        self.n_theta = config.n_theta
        #metríscher vorfaktor wird einmal berechnet und gespeichert als immer wieder neu für rechenzeit
        self.vorfaktor = vorfaktor(self.r)

        # rand indizes 
        self.i_wall = 0 
        self.i_far = self.n_xi - 1

        #in und outflow
        #die toleranz ist noetig: cos(90°) ergibt +6e-17, cos(270°) aber -1.8e-16.
        #ohne toleranz waere oben ausstrom und unten einstrom, also bekaemen die
        #beiden spiegelpunkte durch reinen rundungszufall verschiedene randbedingungen
        self.is_inflow = np.cos(self.theta) < -1e-12
        self.is_outflow = ~self.is_inflow
    def flatten(self, field2d):
        """flacht ein 2D Feld in ein 1D Feld ab, das die gleiche
        Reihenfolge wie die 2D Matrix hat (d.h. Zeilenweise)"""
        return field2d.reshape(-1)
    
    def unflatten(self, field1d):
        return field1d.reshape(self.n_xi, self.n_theta)

    def kartesisch(self):
        return rücktrafo(self.r, self.theta)


    




        