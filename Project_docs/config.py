#config

from dataclasses import dataclass 

@dataclass 
class Config:
    R: float           #zylinderradius
    r_max: float       #äußerer rand 

    U_inf: float       #Anstromgeschwindigkeit
    Re: float          #reynolds zahl (D =2R)
  
   #Auflösung
    n_xi: int          #anzahl der gitterpunkte in radialer richtung (Xi-Rchtung)
    n_theta: int       #anzahl der gitterpunkte in umfangsrichtung (Theta-Richtung)
    
    dt: float          #maximale zeitschrittweite - den tatsaechlichen schritt bestimmt die cfl-bedingung,
                       #ein zu kleines dt hier bremst die rechnung unnoetig aus

    cfl_target: float = 0.5   #cfl wert (numerischer sicherheitswert)

    @property
    def D(self) -> float:    #zylinderdurchmesser für Re
      return 2.0 * self.R 
    
    @property 
    def nu(self) -> float:   #kinematische viskosität aus Re = U_inf * D / nu
      return self.U_inf * self.D / self.Re 


    
