"""
fuehrt alle schnellen pruefprogramme dieses ordners nacheinander aus.

  python Validierung/selbsttest_alle.py

rueckgabewert 0, wenn alles bestanden ist, sonst 1 - der aufruf laesst sich damit auch
automatisiert verwenden. die benchmark-SERIEN gehoeren nicht dazu, die rechnen minuten bis
stunden (python Validierung/benchmark.py); hier laeuft nur, was in sekunden fertig ist.
"""

import sys
import time

from pruefung import Pruefung      #setzt den suchpfad auf Project_docs

import benchmark
import selbsttest_loeser
import selbsttest_operatoren
import taylor_green

PROGRAMME = (
    ("operatoren", selbsttest_operatoren.main),
    ("loeser", selbsttest_loeser.main),
    ("taylor-green", taylor_green.main),
    ("benchmark (auswertung)", benchmark.selbsttest),
)


def main():
    ergebnisse = []
    for name, ausfuehren in PROGRAMME:
        print("=" * 72)
        print(f"  {name}")
        print("=" * 72)
        uhr = time.perf_counter()
        bestanden = bool(ausfuehren())
        ergebnisse.append((name, bestanden, time.perf_counter() - uhr))
        print()

    print("=" * 72)
    for name, bestanden, dauer in ergebnisse:
        print(f"  [{'ok' if bestanden else 'FEHLER'}] {name:24s} {dauer:5.1f} s")
    alles_ok = all(b for _, b, _ in ergebnisse)
    print("=" * 72)
    print("alle selbsttests " + ("bestanden" if alles_ok else "FEHLGESCHLAGEN"))
    return alles_ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
