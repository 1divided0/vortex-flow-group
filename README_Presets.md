# Zylinderumströmung – Programmcode

Finite-Differenzen-Löser für die 2D-inkompressible Umströmung eines Kreiszylinders
(Wirbelstärke–Stromfunktion, logarithmische Polarkoordinaten, RK4) mit Kármán'scher
Wirbelstraße bei Re = 100.

Voraussetzungen: Python 3 mit `numpy` (≥ 2.0), `scipy`, `matplotlib`
(getestet mit Python 3.14.4, numpy 2.4.6, scipy 1.17.1, matplotlib 3.10.9).
Alle Befehle werden im Ordner `Project_docs` ausgeführt.

## Dateien

| Datei | Inhalt |
|---|---|
| `gitter.py` | `Config` (Geometrie, Re, Auflösung, `cfl_target`), log-polares Gitter, `Domain` (Randindizes, Ein-/Ausströmrand, flatten/unflatten) |
| `operatoren.py` | Sparse-Differenzenmatrizen und Laplace-Operator, Selbsttest |
| `loeser.py` | Randbedingungen (Wand/Thom, Fernfeld), Poisson-Löser ∇²ψ = −ω (LR-Zerlegung und FFT, siehe unten), Geschwindigkeiten, rechte Seite, CFL-Zeitschritt, RK4 |
| `main.py` | Anfangszustand, Zeitschleife `eine_Schleife`, Presets `EINZELLAEUFE` |
| `Animation.py` | Auswertung eines Einzellaufs: Wirbelstärke, Strouhal-Zahl, FTLE |
| `benchmark/benchmark.py` | Benchmark-Serien (Orts-/Zeitauflösung, Gebietsgröße) |
| `benchmark/presets.py` | Voreinstellungen der Benchmark-Serien |

Aufrufreihenfolge der Module: `gitter` → `operatoren` → `loeser` → `main` → `Animation` bzw.
`benchmark`. `benchmark/benchmark.py` lässt sich aus jedem Verzeichnis starten; die
Ergebnisse landen immer in `Project_docs/ergebnisse/benchmark/`.

## Poisson-Löser

Für ∇²ψ = −ω gibt es zwei Löser, die **dasselbe** Gleichungssystem lösen (Abweichung ~10⁻¹²).
`wähle_poisson_löser` entscheidet allein nach der Gittergröße; die Ergebnisse ändern sich dadurch
nicht, nur die Rechenzeit:

| Löser | verwendet | Prinzip |
|---|---|---|
| `PoissonSolver` | bis 20 000 Unbekannte | eine LR-Zerlegung des ganzen Systems (`splu`) |
| `PoissonSolverFFT` | darüber | FFT in θ (D2_theta ist zyklisch) → ein kleines 1D-System je Fourier-Mode |

Gemessen auf einem Ryzen 7 9800X3D, Zeit für **eine** Poisson-Lösung:

| Gitter | LR | FFT | |
|---|---|---|---|
| 81×160 | 0,62 ms | 0,76 ms | LR schneller |
| 121×240 | 1,65 ms | 1,26 ms | FFT schneller |
| 161×320 | 3,36 ms | 1,90 ms | FFT fast doppelt so schnell |

Der FFT-Löser braucht zudem deutlich weniger Speicher (0,7 statt 14 MB bei 121×240) und ist
schneller aufgebaut. Er setzt ein **gleichmäßiges, periodisches θ-Gitter** voraus – wird daran
etwas geändert, gilt nur noch `PoissonSolver`. `python loeser.py` prüft beide gegeneinander.

## Einzelläufe (`main.py`)

```
python main.py            # = python main.py schnell
python main.py lang
```

| Preset | Gitter | Zeit | Snapshots | Dauer | Zweck |
|---|---|---|---|---|---|
| `schnell` | 50×100 | t = 0 … 12 | alle 10 Schritte | ≈ 1 s | Funktionstest |
| `lang` | 160×320 | t = 0 … 100 | alle 20 Schritte ab t = 60 | ≈ 15 min | Produktionslauf für `Animation.py` |

Beide schreiben nach `simulation_snapshots.npz` – der Schnelltest überschreibt also einen
vorhandenen Produktionslauf.

**Achtung bei `lang`:** Das Gitter 160×320 erzwingt dt ≈ 1,2·10⁻³, also rund 83 000 Zeitschritte
und etwa 1 660 Snapshots. Die Snapshot-Datei wird damit **ca. 680 MB** groß, und `Animation.py`
lädt sie als float64 – dafür sind mehrere GB Arbeitsspeicher nötig. Wer das vermeiden will, setzt
`Snapshotrange` höher (z. B. 100 statt 20 → ca. 140 MB) oder rechnet auf 80×160.

Danach:

```
python Animation.py       # schreibt ergebnisse/wirbelstaerke.png/.gif, ftle.png/.gif
```

## Selbsttests

```
python operatoren.py                  # Ableitungsmatrizen, Laplace, Poisson gegen Potentialströmung
python loeser.py                      # FFT-Löser gegen LR-Löser (lösen sie dasselbe System?)
python benchmark/benchmark.py selbsttest   # Sonde und Auswertung des Benchmarks an bekannten Fällen
```

## Benchmark (`benchmark/benchmark.py`)

Die Serien stehen in `benchmark/presets.py`, die Auswertung in `benchmark/benchmark.py`.
Der Benchmark rechnet Serien von Läufen, in denen jeweils **ein** Parameter verändert wird,
und bestimmt für jeden Lauf Genauigkeits- und Aufwandsgrößen. Die Serien `gitter`, `zeit` und
`gebiet` verändern nur die **Numerik** (bei Re = 100) und prüfen damit, ob die Lösung
gitter-, zeitschritt- und gebietsunabhängig ist; die Serie `reynolds` verändert die **Physik**
und vergleicht mit Literaturwerten. Alle Läufe gehen bis t = 150. Basis ist das Gitter 81×160, r_max = 20 D,
`cfl_target` = 0,5 (81 statt 80 Punkte, damit sich die Gitterweite bei 41 → 81 → 161 exakt
halbiert).

### Aufruf

```
python benchmark/benchmark.py                          # Übersicht über alle Serien und Stufen
python benchmark/benchmark.py gitter                   # Serie "gitter", Stufe "schnell"
python benchmark/benchmark.py gitter --stufe voll      # alle Läufe der Serie
python benchmark/benchmark.py alle --stufe mittel      # alle Serien bis Stufe "mittel"
python benchmark/benchmark.py zeit --parallel 4        # 4 Läufe gleichzeitig (siehe Hinweise)
python benchmark/benchmark.py gitter --nur-auswerten   # nur Tabelle, CSV und Diagramm neu erstellen
python benchmark/benchmark.py gitter --neu             # vorhandene Ergebnisse neu rechnen
```

Weitere Optionen: `--t-end` (Simulationsende, Standard 150), `--ausgabe` (Ergebnisordner).

### Stufen

Jede Stufe enthält auch alle billigeren Stufen:
`schnell` ⊂ `mittel` ⊂ `voll`. Fertige Läufe werden gespeichert und beim nächsten Aufruf
übersprungen; ein abgebrochener Benchmark setzt also dort fort, wo er aufgehört hat. Läufe mit
gleicher Konfiguration (z. B. 81×160 mit `cfl_target` = 0,5 in allen drei Serien) werden nur
einmal gerechnet.

### Serien

Laufzeiten seriell auf AMD Ryzen 7 9800X3D (gemessen bzw. aus der Schrittzeit hochgerechnet).

**`gitter` – Ortsauflösung.** Gitter verfeinern bei festem Seitenverhältnis
n_theta = 2·(n_xi − 1), sodass die Zellen ihre Form behalten.

| Stufe | Gitter | Zellen an der Wand in 0,1 D | Dauer |
|---|---|---|---|
| schnell | 41×80 | 2,2 | ≈ 6 s |
| schnell | 61×120 | 3,3 | ≈ 30 s |
| mittel | 81×160 | 4,3 | ≈ 1,5 min |
| voll | 121×240 | 6,5 | ≈ 10 min |
| voll | 161×320 | 8,7 | ≈ 35 min |

Aus 41×80, 81×160 und 161×320 (Gitterweite jeweils halbiert) berechnet die Auswertung die
beobachtete Konvergenzordnung p und den nach Richardson auf h → 0 extrapolierten Wert.

**`zeit` – Zeitauflösung.** Der Zeitschritt wird über `cfl_target` gesteuert (`cfg.dt` ist nur
eine Obergrenze und greift nicht). Gitter 81×160.

| Stufe | `cfl_target` | Dauer |
|---|---|---|
| voll | 0,125 | ≈ 6 min |
| mittel | 0,25 | ≈ 3 min |
| schnell | 0,5 · 1,0 · 2,0 | ≈ 1,5 min · 50 s · 25 s |
| mittel | 1,5 · 2,5 · 3,0 | ≈ 1 min zusammen |
| voll | 4,0 | < 1 min |

Zu große Werte führen nicht zwingend zum Abbruch, sondern **verfälschen** die Lösung: Bei
`cfl_target` = 2,0 läuft die Rechnung durch, aber der Ablösewinkel fällt von 116,7° auf 99°
und die Amplitude steigt um 5 %. Das passt zur Stabilitätsgrenze von RK4 für die Diffusion
(λ·dt ≤ 2,785): mit der Zeitschrittformel aus `cfl_zeitschritt` ist an der Wand
λ·dt ≈ 1,73·`cfl_target`, grob abgeschätzt also stabil bis `cfl_target` ≈ 1,6 … 1,8 (die Stufe
`mittel` mit 1,5 und 2,5 grenzt das ein). Bei noch größeren Werten
erhöhen die Fehler die Geschwindigkeiten, und die CFL-Bedingung verkleinert den Zeitschritt
wieder – `dt_mittel` wächst dann nicht mehr mit `cfl_target`. Ein verfälschter Lauf kann also
den Status `periodisch` haben; maßgeblich ist der Vergleich mit den kleineren `cfl_target`.

**`gebiet` – Gebietsgröße.** Äußerer Rand r_max bei festen Gitterweiten dξ, dθ (wie 81×160);
nur n_xi ändert sich. Zeigt den Einfluss der Fernfeld-Randbedingung.

| Stufe | n_xi × n_theta | r_max | Dauer |
|---|---|---|---|
| mittel | 41×160 | 3,2 D | ≈ 45 s |
| schnell | 61×160 | 8,0 D | ≈ 70 s |
| schnell | 81×160 | 20 D | (wie oben) |
| mittel | 101×160 | 50,3 D | ≈ 2 min |
| voll | 121×160 | 126,5 D | ≈ 2,5 min |

**`reynolds` – Reynolds-Zahl.** Einzige physikalische Serie: Gitter 81×160, r_max = 20 D,
`cfl_target` = 0,5, nur Re ändert sich. Unterhalb von Re ≈ 47 bildet sich ein **stationäres**
Wirbelpaar, darüber löst der Nachlauf periodisch ab. Kleines Re bedeutet großes ν, der
Zeitschritt ist diffusionsbegrenzt und wird dadurch kleiner – Re = 20 kostet etwa das
Fünffache von Re = 100.

| Stufe | Re | Erwartung | Dauer |
|---|---|---|---|
| mittel | 20 | stationär, L/D ≈ 0,93, Ablösewinkel ≈ 136° | ≈ 8 min |
| schnell | 40 | stationär, L/D ≈ 2,24, Ablösewinkel ≈ 126,5° | ≈ 4 min |
| schnell | 60 | Wirbelstraße, St ≈ 0,136 | ≈ 3 min |
| schnell | 100 | Wirbelstraße, St ≈ 0,164 | (wie oben) |
| mittel | 150 | Wirbelstraße, St ≈ 0,183 | ≈ 1,5 min |
| voll | 200 | St ≈ 0,197; Auflösung wird knapp | ≈ 1,5 min |

Vergleichswerte: die Kurve St(Re) = −3,3265/Re + 0,1816 + 1,6·10⁻⁴·Re (Williamson 1988,
gültig 47 < Re < 180) und für das stationäre Wirbelpaar Rückströmlänge und Ablösewinkel nach
Coutanceau & Bouard (1977) bzw. Fornberg (1980). Beides wird im Diagramm automatisch als
Literaturvergleich eingezeichnet (gestrichelte Kurve bzw. Sterne). Die Grenzschichtdicke
skaliert mit Re^(−1/2): das Gitter 81×160 ist für Re = 200 nur noch knapp ausreichend, der
Lauf gehört mit zur Aussage.

**Gesamtdauer** (alle Serien, gemeinsame Läufe nur einmal): Stufe `schnell` ≈ 4,5 min (gemessen),
`mittel` ≈ 12 min, `voll` ≈ 65 min.

### Ergebnisse

In `ergebnisse/benchmark/`:
- `laeufe/<lauf>.json` – alle Kennzahlen eines Laufs inkl. Konfiguration und Rechnerdaten,
  `laeufe/<lauf>.npz` – Sondensignal (jeder Zeitschritt), gemittelte Wandwirbelstärke und
  gemittelte Längsgeschwindigkeit auf der Nachlaufachse (mit `theta` bzw. `xi` dazu)
- `<serie>.csv` – Tabelle aller vorhandenen Läufe der Serie
- `<serie>.png` – Diagramm: drei Genauigkeitsgrößen und der Rechenaufwand über dem Parameter.
  Standard sind St, Amplitude und Ablösewinkel; eine Serie kann in `presets.py` über `panels`
  andere wählen (`reynolds` zeigt statt der Amplitude die Rückströmlänge, weil es die auch ohne
  Ablösung gibt) und über `referenz` Literaturwerte eintragen

### Messgrößen

| Größe | Bedeutung |
|---|---|
| `status` | `periodisch`, `nicht_periodisch` (Anlauf nicht abgeschlossen), `keine_abloesung`, `instabil` (Abbruch) |
| `St` | Strouhal-Zahl D/(U·T) aus der mittleren Periode T im Auswertefenster |
| `periode_streuung` | relative Standardabweichung der Perioden im Fenster |
| `amplitude` | Amplitude der Quergeschwindigkeit u_θ an der Sonde (r = 2 D, θ = 0, x = 1,5 D hinter dem Zylinder) |
| `rueckstroemlaenge` | Länge L/D des Rückströmgebiets ab der Zylinderrückseite: erster Vorzeichenwechsel der zeitgemittelten Längsgeschwindigkeit auf der Nachlaufachse (θ = 0), linear in ξ interpoliert. `nan`, wenn die Strömung direkt hinter dem Zylinder schon nach außen zeigt |
| `restschwankung` | größte Abweichung des Sondensignals vom Endwert in den letzten 10 Zeiteinheiten, bezogen auf U. Klein (< 10⁻³) heißt: der Lauf ist stationär eingelaufen, `keine_abloesung` ist dann das physikalische Ergebnis und kein unfertiger Anlauf |
| `fenster_art` | `perioden` (Mittelung über das Periodenfenster) oder `ende` (über die letzten 10 Zeiteinheiten, wenn es keine Perioden gibt) |
| `abloesewinkel` | zeitlich gemittelter Ablösewinkel, gemessen vom vorderen Staupunkt (Vorzeichenwechsel der mittleren Wandwirbelstärke); `_oben`/`_unten` einzeln |
| `t_einsatz` | Zeitpunkt, an dem das Sondensignal die halbe Endamplitude erreicht |
| `n_perioden`, `t_fenster_start/-ende` | Auswertefenster |
| `schritte`, `dt_mittel` | Anzahl und mittlere Größe der Zeitschritte |
| `t_aufbau_s`, `t_schleife_s`, `ms_pro_schritt` | Rechenzeit für Aufbau (inkl. LR-Zerlegung) und Zeitschleife |
| `t_loese_ms`, `poisson_anteil` | Zeit einer Poisson-Lösung, Anteil der 4 Lösungen pro Schritt an der Schrittzeit |
| `poisson_loeser` | `PoissonSolver` (LR) oder `PoissonSolverFFT`; in der Tabelle als `LR` / `FFT` |
| `rechenzeit_pro_periode_s` | Rechenzeit für eine Ablöseperiode |
| `nnz_LU`, `speicher_LU_MB` | Größe der Zerlegung (beim FFT-Löser Summe über alle Moden) |

**Auswertefenster:** Ablösewinkel und Rückströmlänge sind Zeitmittel. Aus den Aufwärts-Nulldurchgängen des Sondensignals werden Periode und
Amplitude jeder Schwingung bestimmt. Das Fenster umfasst die letzten Perioden, deren Dauer um
höchstens 1 % und deren Amplitude um höchstens 2 % vom Median der letzten drei Perioden
abweicht – der Anlauf fällt damit automatisch heraus. `periodisch` erfordert mindestens
4 Perioden im Fenster.

### Hinweise

- **Parallelbetrieb:** Ein einzelner Lauf nutzt genau einen CPU-Kern. Mit `--parallel N`
  laufen N Rechnungen gleichzeitig (4 parallele Läufe ≈ 3,2-facher, 8 ≈ 4,2-facher Durchsatz).
  Die Läufe bremsen sich dabei gegenseitig, die **Zeitmessungen sind dann nicht vergleichbar**
  und werden als `zeitmessung_vergleichbar = false` markiert. Genauigkeitsgrößen (St usw.) sind
  davon nicht betroffen. Für die Aufwandsauswertung seriell rechnen und den Rechner währenddessen
  nicht anderweitig belasten.
- Die Sonde sitzt auf jedem Gitter an exakt derselben Stelle (lineare Interpolation in ξ), die
  Werte verschiedener Gitter sind also direkt vergleichbar.
- Die Messung ändert die Rechnung nicht (Endzustand bitgleich mit und ohne Messung) und kostet
  weniger als 1 % Rechenzeit.
