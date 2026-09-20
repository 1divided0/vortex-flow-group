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
| `signalauswertung.py` | Periodenfenster, Strouhal-Zahl, Amplitude, Restschwankung und Sondenlage – **gemeinsam** von `Animation.py` und dem Benchmark benutzt |
| `Animation.py` | Auswertung eines Einzellaufs: Wirbelstärke, Strouhal-Zahl, FTLE |
| `Validierung/benchmark.py` | Benchmark-Serien (Orts-/Zeit-/Gebiets-/Kreuz-/Reynolds-Serie) |
| `Validierung/presets.py` | Voreinstellungen der Benchmark-Serien |
| `Validierung/taylor_green.py` | Validierung am Taylor-Green-Wirbel (exakte Lösung) |
| `Validierung/mms_logpolar.py` | Manufactured Solution auf dem log-polaren Gitter: Metrik, Advektion, Diffusion, Thom-Formel |
| `Validierung/selbsttest_operatoren.py` | Differenzenmatrizen gegen analytische Ableitungen |
| `Validierung/selbsttest_loeser.py` | Poisson-Löser gegen die Potentialströmung, FFT gegen LR |
| `Validierung/selbsttest_alle.py` | führt alle schnellen Prüfprogramme aus |
| `Validierung/pruefung.py` | gemeinsames Gerüst der Prüfprogramme |

**Löser und Tests sind getrennt:** `gitter.py`, `operatoren.py`, `loeser.py` und `main.py`
enthalten nur rechnenden Code, alle Prüfprogramme liegen in `Validierung/` und rufen den
Löser von außen auf.

Aufrufreihenfolge der Module: `gitter` → `operatoren` → `loeser` → `main` → `Animation` bzw.
`Validierung`. Alles in `Validierung/` lässt sich aus jedem Verzeichnis starten; die
Ergebnisse der Benchmark-Serien landen immer in `Project_docs/ergebnisse/benchmark/`.

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
etwas geändert, gilt nur noch `PoissonSolver`. `Validierung/selbsttest_loeser.py` prüft beide
gegeneinander.

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

## Selbsttests und Validierung (`Validierung/`)

Alle Prüfprogramme zusammen laufen in unter zwei Sekunden und liefern den Rückgabewert 0,
wenn alles bestanden ist:

```
python Validierung/selbsttest_alle.py
```

Einzeln:

```
python Validierung/selbsttest_operatoren.py   # Differenzenmatrizen gegen analytische Ableitungen
python Validierung/selbsttest_loeser.py       # Poisson gegen Potentialströmung, Ordnung, FFT gegen LR
python Validierung/taylor_green.py            # Taylor-Green-Wirbel: exakte Lösung, Ordnung, nu = 0
python Validierung/mms_logpolar.py            # Manufactured Solution auf dem log-polaren Gitter
python Validierung/benchmark.py selbsttest    # Sonde und Auswertung des Benchmarks
```

### Taylor-Green-Wirbel

Das Pflichtbeispiel der Aufgabenstellung und der einzige Test, der Advektion,
Poisson-Lösung und RK4 **gemeinsam** gegen eine exakte Lösung der nichtlinearen Gleichungen
stellt. Die log-polaren Gleichungen enthalten den kartesischen Fall als exakten Spezialfall:
mit r = 1 (metrischer Vorfaktor 1) und ξ = x, θ = y werden Laplace, Geschwindigkeiten,
Advektion und Zeitschrittformel Zeile für Zeile zu den kartesischen. `geschwindigkeit`,
`berechne_rhs`, `cfl_zeitschritt` und `rk4` laufen deshalb **unverändert** – es ist derselbe
Code, der auch die Zylinderumströmung rechnet. Ersetzt werden nur der Poisson-Löser (periodisch
statt Dirichlet) und die Randbedingungen (`rk4` bekommt dafür eine leere Funktion übergeben).

| Prüfung | Ergebnis |
|---|---|
| Poisson-Löser (periodisch), Residuum | 7·10⁻¹⁵ |
| Advektionsterm gegen analytische Ableitung, Ordnung | 1,91 / 1,98 |
| Taylor-Green ν = 0,1, rel. L2-Fehler ω (16/32/64) | 1,3·10⁻³ / 3,2·10⁻⁴ / 8,0·10⁻⁵, Ordnung 2,00 |
| Abklingrate gegen 2ν | −0,08 % |
| ν = 0: stationäre Lösung nach 100 Schritten | 1·10⁻¹⁷ |
| ν = 0: Enstrophie-Drift | exakt 0 |

**Was der Test nicht abdeckt:** Mit r = 1 und periodischen Rändern sind die Metrik 1/r²,
die Thom-Formel, die Fernfeld-Randbedingung und die einseitigen Randzeilen von
`build_D_xi`/`build_D2_xi` abgeschaltet. Und weil die Nichtlinearität im Taylor-Green-Fall
identisch verschwindet, sagt er über die Genauigkeit des Advektionsterms nichts aus – dafür
steht der separate Test mit einer Lösung, bei der ω kein Vielfaches von ψ ist. Die
zylinderspezifischen Teile decken `selbsttest_loeser.py` (Abnahmetest gegen die
Potentialströmung), `mms_logpolar.py` (siehe unten) und die Benchmark-Serie `reynolds`
(Literaturvergleich) ab.

### Manufactured Solution auf dem log-polaren Gitter

Schließt genau die Lücke, die der Taylor-Green-Test bauartbedingt lässt: Statt eine Lösung zu
suchen, wird eine beliebige glatte Funktion **vorgegeben**, in die diskreten Operatoren
eingesetzt und mit der analytisch differenzierten verglichen. Das läuft auf dem Gitter der
Zylinderrechnung, die Metrik 1/r² ist also aktiv. ω ist bewusst kein Vielfaches von ψ, damit
der Advektionsterm nicht wie beim Taylor-Green-Wirbel identisch verschwindet.

| Prüfung | Ergebnis |
|---|---|
| `berechne_rhs` (Metrik + Advektion + Diffusion), n_xi = 41/81/161 | 7,9·10⁻³ / 2,0·10⁻³ / 5,0·10⁻⁴ → **Ordnung 2,00** |
| nur Advektion (ν = 0) | Ordnung 2,00 |
| nur Diffusion (ψ = 0) | Ordnung 2,00 |
| Thom-Formel an der Wand | 8,8·10⁻² / 4,5·10⁻² / 2,3·10⁻² → **Ordnung 0,97** |

Die letzte Zeile ist kein Fehler, sondern der Nachweis der bekannten Eigenschaft: Die
Thom-Formel ist erster Ordnung und damit die einzige Stelle des Lösers, die nicht zweiter
Ordnung ist. Ergänzend prüft `selbsttest_operatoren.py` die einseitigen Randzeilen von
`build_D_xi`/`build_D2_xi` an einer nicht-polynomialen Funktion (exp) – an ξ³ allein sind
die einseitigen Formeln exakt und ihre Ordnung bliebe ungeprüft.

## Benchmark (`Validierung/benchmark.py`)

Die Serien stehen in `Validierung/presets.py`, die Auswertung in `Validierung/benchmark.py`.
Der Benchmark rechnet Serien von Läufen, in denen jeweils **ein** Parameter verändert wird,
und bestimmt für jeden Lauf Genauigkeits- und Aufwandsgrößen. Die Serien `gitter`, `zeit`,
`gebiet` und `gitter_fern` verändern nur die **Numerik** (bei Re = 100) und prüfen damit, ob die
Lösung gitter-, zeitschritt- und gebietsunabhängig ist; die Serie `reynolds` verändert die
**Physik** und vergleicht mit Literaturwerten. Alle Läufe gehen bis t = 150. Basis ist das Gitter 81×160, r_max = 20 D,
`cfl_target` = 0,5 (81 statt 80 Punkte, damit sich die Gitterweite bei 41 → 81 → 161 exakt
halbiert).

### Aufruf

```
python Validierung/benchmark.py                          # Übersicht über alle Serien und Stufen
python Validierung/benchmark.py gitter                   # Serie "gitter", Stufe "schnell"
python Validierung/benchmark.py gitter --stufe voll      # alle Läufe der Serie
python Validierung/benchmark.py alle --stufe mittel      # alle Serien bis Stufe "mittel"
python Validierung/benchmark.py zeit --parallel 4        # 4 Läufe gleichzeitig (siehe Hinweise)
python Validierung/benchmark.py gitter --nur-auswerten   # nur Tabelle, CSV und Diagramm neu erstellen
python Validierung/benchmark.py gitter --neu             # vorhandene Ergebnisse neu rechnen
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

**`gitter_fern` – Kreuzprobe Gitter × Gebiet.** Dieselbe Verfeinerung wie `gitter`, aber auf
dem größeren Gebiet r_max = 50,3 D (dξ und dθ jeweils wie im zugehörigen 20-D-Lauf, das
20-D-Gitter ist genau der innere Teil des größeren).

| Stufe | n_xi × n_theta | entspricht in `gitter` | Dauer |
|---|---|---|---|
| schnell | 51×80 | 41×80 | ≈ 10 s |
| schnell | 76×120 | 61×120 | ≈ 40 s |
| mittel | 101×160 | 81×160 | (= `gebiet` 50,3 D) |
| voll | 151×240 | 121×240 | ≈ 10 min |
| voll | 201×320 | 161×320 | ≈ 40 min |

**Wozu die Serie da ist.** `gitter` und `gebiet` allein sagen nichts darüber, ob sich Orts- und
Gebietsfehler **addieren**. Täten sie es, wäre der Abstand zwischen dieser Kurve und der von
`gitter` konstant. Gemessen wächst er:

| n_theta | St bei 20 D | St bei 50,3 D | Differenz |
|---|---|---|---|
| 80 | 0,156894 | 0,155171 | −0,00172 |
| 120 | 0,162064 | 0,160251 | −0,00181 |
| 160 | 0,164418 | 0,162007 | −0,00241 |
| 240 | 0,166464 | 0,163060 | −0,00340 |
| 320 | 0,166926 | 0,163361 | −0,00357 |

Die Fehler sind also **nicht separierbar**, und dass sich bei 81×160 / 20 D ein Wert nahe der
Literatur ergibt, ist kein Argument, sondern auflösungsabhängige Kompensation. Entscheidend ist
die Gitterkonvergenz auf dem großen Gebiet (Richardson 80/160/320):

| | p | St(h → 0) | vs. Williamson 0,16434 |
|---|---|---|---|
| `gitter` (20 D) | 1,58 | 0,1682 | +2,4 % |
| `gitter_fern` (50,3 D) | **2,34** | **0,1637** | **−0,4 %** |

Die dominante Fehlerquelle ist damit die Fernfeld-Randbedingung, nicht das Gitter – und das
p = 1,58 der Serie `gitter` entsteht dadurch, dass der Gebietsfehler selbst gitterabhängig ist
und die Richardson-Extrapolation verunreinigt. Auf dem großen Gebiet erreicht der Löser
praktisch die nominelle Ordnung 2.

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

**Gesamtdauer** (alle Serien, gemeinsame Läufe nur einmal): Stufe `schnell` ≈ 5 min,
`mittel` ≈ 13 min, `voll` ≈ 2 h (davon rund 75 min für die beiden feinsten Läufe
161×320 und 201×320).

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
  Werte verschiedener Gitter sind also direkt vergleichbar. `Animation.py` benutzt über
  `signalauswertung.py` dieselbe Sondenlage **und** dasselbe Periodenfenster; Bild und Benchmark
  liefern damit dieselbe Messvorschrift (Differenz für 81×160: 0,02 %, Periodenstreuung 0,1 %).
- **Ergebnisformat:** `format_version` in der json steht auf 2. Ältere Läufe (Version 1) werden
  nicht übersprungen, sondern neu gerechnet – ihnen fehlen `rueckstroemlaenge`, `fenster_art`,
  `restschwankung` und `poisson_loeser`, und die Rückströmlänge lässt sich aus der alten npz
  nicht nachtragen (sie enthält `xi` und `achse_mittel` noch nicht).
- Die Messung ändert die Rechnung nicht (Endzustand bitgleich mit und ohne Messung) und kostet
  weniger als 1 % Rechenzeit.
