"""
gemeinsames geruest der pruefprogramme in diesem ordner.

die tests liegen bewusst ausserhalb des loesers: gitter.py, operatoren.py, loeser.py und
main.py enthalten nur rechnenden code, jedes pruefprogramm hier ruft ihn nur auf. wer den
loeser aendert, sieht mit einem lauf von selbsttest_alle.py, ob er etwas kaputt gemacht hat.

dieses modul zu importieren macht den loeser importierbar, egal aus welchem ordner heraus
ein pruefprogramm gestartet wird.
"""

import os
import sys

#der loeser liegt eine ebene hoeher (Project_docs)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

#grad- und pfeilzeichen sollen auch dann funktionieren, wenn die ausgabe in eine datei
#umgeleitet wird - windows faellt sonst auf cp1252 zurueck und bricht mit UnicodeEncodeError ab
try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, OSError):
    pass


class Pruefung:
    """
    sammelt einzelne pruefungen und merkt sich, ob alle bestanden sind.
    jede pruefung druckt ihren messwert mit - eine zahl zu sehen ist oft nuetzlicher
    als nur "ok", und bei einem fehlschlag steht der vergleich direkt daneben
    """

    def __init__(self, titel=None):
        self.alles_ok = True
        if titel:
            print(titel)

    def abschnitt(self, name):
        print(f"{name}:")

    def pruefe(self, bedingung, text):
        bedingung = bool(bedingung)
        self.alles_ok &= bedingung
        print(f"  [{'ok' if bedingung else 'FEHLER'}] {text}")
        return bedingung

    def fazit(self, name):
        print(f"{name} " + ("bestanden" if self.alles_ok else "FEHLGESCHLAGEN"))
        return self.alles_ok
