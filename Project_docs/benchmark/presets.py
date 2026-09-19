#voreinstellungen der benchmark-serien, beschrieben in README_Presets.md
#(die einzellaeufe von main.py stehen dort in EINZELLAEUFE)

import math
import os
import sys
from dataclasses import replace

#der loeser liegt eine ebene hoeher (Project_docs). der eintrag macht ihn importierbar,
#egal aus welchem ordner heraus benchmark/benchmark.py aufgerufen wird
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gitter import Config


# ---------------------------------------------------------------------------
# benchmark-serien fuer benchmark.py:  python benchmark.py <serie> --stufe <stufe>
# ---------------------------------------------------------------------------

#basis aller serien. n_xi - 1 = 80 statt 79, damit sich die gitterweite bei
#41 -> 81 -> 161 punkten exakt halbiert (noetig fuer die richardson-extrapolation)
BENCHMARK_BASIS = Config(R=0.5, r_max=20.0, U_inf=1.0, Re=100.0, n_xi=81, n_theta=160, dt=0.05, cfl_target=0.5)

#lang genug, damit auch grobe gitter nach dem anlauf noch viele perioden liefern
BENCHMARK_T_END = 150.0

#eine stufe enthaelt immer auch alle billigeren stufen
STUFEN = ("schnell", "mittel", "voll")


def _gebiet(n_xi):
    #r_max so waehlen, dass dxi exakt wie in der basis bleibt: dxi = ln(r_max/R) / (n_xi - 1)
    b = BENCHMARK_BASIS
    #gerundet, damit fuer n_xi = 81 exakt die basis (r_max = 20.0) herauskommt und
    #der lauf mit dem gleichen lauf der anderen serien zusammenfaellt
    dxi = math.log(b.r_max / b.R) / (b.n_xi - 1)
    return dict(n_xi=n_xi, r_max=round(b.R * math.exp(dxi * (n_xi - 1)), 10))


SERIEN = {
    "gitter": dict(
        beschreibung="ortsaufloesung: gitter verfeinern, seitenverhaeltnis n_theta = 2*(n_xi - 1) bleibt fest",
        x=lambda cfg: cfg.n_theta, x_name="n_theta",
        laeufe=[
            dict(stufe="schnell", aenderung=dict(n_xi=41, n_theta=80)),
            dict(stufe="schnell", aenderung=dict(n_xi=61, n_theta=120)),
            dict(stufe="mittel", aenderung=dict(n_xi=81, n_theta=160)),
            dict(stufe="voll", aenderung=dict(n_xi=121, n_theta=240)),
            dict(stufe="voll", aenderung=dict(n_xi=161, n_theta=320)),
        ],
    ),
    "zeit": dict(
        #oberhalb von etwa cfl_target = 2 bricht die rechnung nicht ab, sondern verfaelscht:
        #fehler erhoehen die geschwindigkeiten, die cfl-bedingung verkleinert dt wieder.
        #erkennbar daran, dass dt_mittel nicht mehr mit cfl_target waechst
        beschreibung="zeitaufloesung: cfl_target auf dem 81x160-gitter, bis die loesung instabil oder verfaelscht wird",
        x=lambda cfg: cfg.cfl_target, x_name="cfl_target",
        laeufe=[
            dict(stufe="voll", aenderung=dict(cfl_target=0.125)),
            dict(stufe="mittel", aenderung=dict(cfl_target=0.25)),
            dict(stufe="schnell", aenderung=dict(cfl_target=0.5)),
            dict(stufe="schnell", aenderung=dict(cfl_target=1.0)),
            dict(stufe="mittel", aenderung=dict(cfl_target=1.5)),
            dict(stufe="schnell", aenderung=dict(cfl_target=2.0)),
            dict(stufe="mittel", aenderung=dict(cfl_target=2.5)),
            dict(stufe="mittel", aenderung=dict(cfl_target=3.0)),
            dict(stufe="voll", aenderung=dict(cfl_target=4.0)),
        ],
    ),
    "gebiet": dict(
        beschreibung="gebietsgroesse: r_max bei festem dxi und dtheta (gitterweiten wie 81x160)",
        x=lambda cfg: cfg.r_max / cfg.D, x_name="r_max / D",
        laeufe=[
            dict(stufe="mittel", aenderung=_gebiet(41)),     #r_max =   3.2 D
            dict(stufe="schnell", aenderung=_gebiet(61)),    #r_max =   8.0 D
            dict(stufe="schnell", aenderung=_gebiet(81)),    #r_max =  20.0 D (basis)
            dict(stufe="mittel", aenderung=_gebiet(101)),    #r_max =  50.3 D
            dict(stufe="voll", aenderung=_gebiet(121)),      #r_max = 126.5 D
        ],
    ),
}


def laeufe_der_serie(serie, stufe):
    #liefert die configs aller laeufe der serie bis einschliesslich der stufe
    grenze = STUFEN.index(stufe)
    return [(lauf["stufe"], replace(BENCHMARK_BASIS, **lauf["aenderung"]))
            for lauf in SERIEN[serie]["laeufe"]
            if STUFEN.index(lauf["stufe"]) <= grenze]
