# Der Lösungsweg: Schritt für Schritt zum eigenen Strömungslöser

Dieses Dokument geht den kompletten Weg von der Aufgabenstellung bis zur
laufenden Zylinderumströmung durch — aber nicht als Ergebnisbericht, sondern
als Gedankengang. Bei jedem Schritt steht im Vordergrund: **Welche Frage
stellt sich hier, warum ist die gewählte Antwort naheliegend, und woran hätte
man selbst darauf kommen können?** Die Theorie dazu steht in den Themen 1–6;
hier geht es um den roten Faden, der sie verbindet.

**Der rote Faden in einem Satz:** Wir formen die Navier-Stokes-Gleichungen so
um, dass ihre unangenehmen Eigenschaften (Druck, Kontinuitätszwang)
verschwinden, übersetzen die verbleibenden Operatoren in Matrizen, zerlegen
jeden Zeitschritt in drei einfache Teilprobleme, stecken die Physik in die
Randbedingungen — und beweisen an einem analytisch lösbaren Fall, dass alles
stimmt, **bevor** wir uns an die Geometrie des Zylinders wagen.

[TOC]


---
## Schritt 0: Das Problem so lange vereinfachen, bis es lösbar wird

**Die Leitfrage:** *Was ist die einfachste Version meines Problems, die den
interessanten Effekt noch enthält?*

Die Wirbelstraße hinter einem Zylinder ist ein Phänomen, das schon in zwei
Dimensionen, bei laminarer, inkompressibler Strömung mit konstanten
Stoffwerten auftritt. Also darf man alles andere weglassen:

- **Kontinuum statt Moleküle** (Knudsen-Zahl klein, Thema 1) → Feldgrößen
  und partielle Differentialgleichungen sind überhaupt erst zulässig.
- **Inkompressibel** → die Dichte ist konstant, die Energiegleichung
  entkoppelt, es bleiben Masse- und Impulserhaltung.
- **2D** → aus drei Geschwindigkeitskomponenten werden zwei; entscheidender
  noch: die Wirbelstärke wird vom Vektor zum Skalar (Schritt 1).

> **Merksatz:** Jede Vereinfachung, die man *vor* dem Programmieren macht,
> spart ein Vielfaches an Arbeit *beim* Programmieren. Die Kunst ist, genau
> so weit zu vereinfachen, dass das Zielphänomen überlebt.

Das Ausgangssystem ist damit (Thema 1):

$$
\boldsymbol{\nabla}\cdot\boldsymbol{u} = 0, \qquad
\partial_t\boldsymbol{u} + (\boldsymbol{u}\cdot\boldsymbol{\nabla})\boldsymbol{u}
= -\frac{1}{\rho}\boldsymbol{\nabla}p + \nu\nabla^2\boldsymbol{u}
$$


---
## Schritt 1: Die Formulierung wählen — Druck und Kontinuität eliminieren

**Die Leitfrage:** *Welche Bestandteile der Gleichungen machen die numerische
Lösung schwer — und gibt es eine Umformung, die sie beseitigt?*

Wer die Gleichungen direkt in $u, v, p$ („primitiven Variablen") lösen will,
stößt auf zwei strukturelle Ärgernisse:

1. **Der Druck hat keine eigene Zeitentwicklungsgleichung.** Er ist kein
   transportiertes Feld, sondern stellt sich instantan so ein, dass die
   Kontinuität erfüllt bleibt. Man müsste ihn in jedem Schritt implizit
   mitbestimmen (Projektionsverfahren & Co.).
2. **Die Kontinuitätsgleichung ist eine Zwangsbedingung**, keine
   Bestimmungsgleichung: Sie sagt nicht, *wie* sich etwas entwickelt, sondern
   *was verboten ist*.

Für beide Probleme gibt es je einen klassischen Handgriff — und beide sollte
man als allgemeine Denkwerkzeuge abspeichern:

**Handgriff 1 — Störendes eliminieren durch einen Differentialoperator:**
Der Druck taucht nur als Gradient auf, und die Rotation eines Gradienten
verschwindet identisch: $\boldsymbol{\nabla}\times\boldsymbol{\nabla}p = 0$.
Wendet man also die Rotation auf die Impulsgleichung an, fällt der Druck
komplett heraus. Übrig bleibt eine Transportgleichung für die Wirbelstärke
$\omega_z = \partial_x v - \partial_y u$ (Thema 2):

$$
\partial_t\omega_z + \boldsymbol{u}\cdot\boldsymbol{\nabla}\omega_z = \nu\nabla^2\omega_z
$$

In 2D ist das ein **Skalarfeld** — und die Gleichung hat die Form der
allgemeinen Konvektions-Diffusions-Gleichung aus Thema 1. Genau deshalb lohnt
es sich, die 1D-Transportgleichung vorher verstanden zu haben: Die
Wirbeltransportgleichung ist ihr zweidimensionaler Bruder.

**Handgriff 2 — Zwangsbedingungen „by construction" erfüllen:** Statt die
Kontinuität in jedem Schritt zu erzwingen, wählt man Variablen, die sie
automatisch erfüllen. Der Ansatz der Stromfunktion $\Psi$ mit

$$
u = \partial_y\Psi, \qquad v = -\partial_x\Psi
$$

liefert $\boldsymbol{\nabla}\cdot\boldsymbol{u} = \partial_x\partial_y\Psi -
\partial_y\partial_x\Psi \equiv 0$ — die Kontinuität ist keine Bedingung
mehr, sondern eine Eigenschaft der Konstruktion (Cauchy-Riemann, Thema 2).

**Die Verbindung beider Größen** liefert die Definition der Wirbelstärke,
in die man den Stromfunktionsansatz einsetzt:

$$
\omega_z = \partial_x v - \partial_y u = -\nabla^2\Psi
$$

Das ist eine **Poisson-Gleichung** — und damit ist das System geschlossen:
zwei Skalarfelder ($\omega_z$, $\Psi$), drei Gleichungen (Poisson,
Cauchy-Riemann, Transport), kein Druck, keine Zwangsbedingung.

> **Wie kommt man selbst darauf?** Man schreibt auf, was einen stört
> (Druck ohne Gleichung, Kontinuität als Zwang), und fragt für jedes Ärgernis:
> *Gibt es einen Operator, unter dem der Störterm verschwindet? Gibt es eine
> Variablenwahl, die den Zwang automatisch erfüllt?* Diese zwei Fragen führen
> fast zwangsläufig zu Rotation und Stromfunktion.


---
## Schritt 2: Vom Kontinuum zum Gitter — Ableitungen werden Matrizen

**Die Leitfrage:** *Wie mache ich aus einem Differentialoperator etwas, das
ein Computer ausführen kann?*

**2a — Taylor-Reihe als Fundament.** Jede finite Differenz entsteht aus
Taylor-Entwicklungen benachbarter Gitterpunkte, die man so kombiniert, dass
die gewünschte Ableitung übrig bleibt und möglichst viele Fehlerterme
wegfallen (Thema 3). Die Zentraldifferenz

$$
\partial_x f \approx \frac{f_{i+1} - f_{i-1}}{2h}
$$

hat Fehlerordnung 2, weil sich die geraden Terme der beiden Entwicklungen
gegenseitig auslöschen. Wichtig für später: Am Gebietsrand fehlen Nachbarn,
dort braucht man **einseitige Differenzen gleicher Ordnung** — sonst
verschmutzt der Randfehler die Konvergenzordnung des Gesamtverfahrens.

**2b — Das Gitter wird ein Vektor, die Ableitung eine Matrix.** Die
entscheidende Umdeutung: Legt man alle Gitterwerte in einen langen Vektor
(zeilenweise „Vektorisierung", Thema 3), dann ist *jede* finite Differenz
eine **lineare Abbildung** — also eine Matrix $\boldsymbol{D}$. Aus
$\partial_x f$ wird $\boldsymbol{D}_x\boldsymbol{f}$, aus $\nabla^2$ wird
$\boldsymbol{D}_x^{(2)} + \boldsymbol{D}_y^{(2)}$. Der ganze PDE-Apparat
schrumpft auf Matrixprodukte und das Lösen linearer Gleichungssysteme.

**2c — Von 1D nach 2D per Kronecker-Produkt.** Man baut nur eindimensionale
Ableitungsmatrizen und hebt sie mit dem Kronecker-Produkt auf das 2D-Gitter:

$$
\boldsymbol{D}_x = I_y \otimes \boldsymbol{d}_x, \qquad
\boldsymbol{D}_y = \boldsymbol{d}_y \otimes I_x
$$

(Reihenfolge passend zur Vektorisierungskonvention!). Das ist der Grund,
warum sich der Baukasten in `src/operators.py` auf den 1D-Fall konzentrieren
kann.

**2d — Dünnbesetzt speichern.** Eine Ableitungsmatrix für ein 320×160-Gitter
hat 51200² ≈ 2,6 Milliarden Einträge, von denen nur ein paar hunderttausend
ungleich null sind. Sparse-Formate machen aus einem unmöglichen Problem ein
triviales (`Sparse_Performance` im Begleitmaterial zeigt die Größenordnungen).

> **Fertig-Kriterium dieses Schritts:** Die Matrizen angewandt auf bekannte
> Funktionen (z. B. $\sin$) reproduzieren deren Ableitungen, und der Fehler
> fällt bei Gitterverfeinerung mit der nominellen Ordnung
> (`tests/test_operators.py`). Erst wenn das steht, weiterbauen — jeder
> spätere Fehler ließe sich sonst nicht mehr zuordnen.


---
## Schritt 3: Der Lösungsalgorithmus — ein nichtlineares Problem in drei lineare Teilschritte zerlegen

**Die Leitfrage:** *In welcher Reihenfolge löse ich die drei gekoppelten
Gleichungen, sodass jede Einzelne einfach wird?*

Zuerst die grundsätzliche Strategie (Thema 4): **Semi-Diskretisierung** —
erst den Ort auflösen, dann bleibt ein System gewöhnlicher
Differentialgleichungen in der Zeit („Linienmethode"). Orts- und
Zeitdiskretisierung lassen sich damit getrennt entwerfen und getrennt
analysieren.

Dann der Kern: Die Nichtlinearität der Navier-Stokes-Gleichungen steckt
allein im Konvektionsterm $\boldsymbol{u}\cdot\boldsymbol{\nabla}\omega_z$ —
Geschwindigkeit mal Wirbelstärkegradient, und beide hängen voneinander ab.
Der Trick: Innerhalb eines Zeitschritts das **eingefrorene** Feld verwenden.
Dann zerfällt jeder Schritt in drei *lineare* Teilprobleme:

1. **Poisson:** $\nabla^2\Psi = -\omega_z$ — aus der aktuellen Wirbelstärke
   die Stromfunktion berechnen. (Ein lineares Gleichungssystem.)
2. **Cauchy-Riemann:** $u = \partial_y\Psi$, $v = -\partial_x\Psi$ — reine
   Matrixprodukte.
3. **Transport:** $\dot{\omega}_z = \nu\nabla^2\omega_z -
   (\boldsymbol{u}\cdot\boldsymbol{\nabla})\omega_z$ mit dem nun *bekannten*
   $\boldsymbol{u}$ — ein linearer Ausdruck in $\omega_z$, den das
   Zeitschrittverfahren integriert. Danach zurück zu 1.

Warum diese Reihenfolge? Weil die Kausalkette der Physik folgt: Die
Wirbelstärke *ist* der Zustand (sie wird transportiert), Stromfunktion und
Geschwindigkeit sind nur abgeleitete Momentaufnahmen davon.

**Der Performance-Blick:** Das teuerste Stück ist das Poisson-System. Aber
seine Matrix ändert sich über die Zeit nie — nur die rechte Seite. Also
zerlegt man sie **einmal** (LR-Zerlegung vor der Zeitschleife) und macht pro
Schritt nur noch Vorwärts-/Rückwärtseinsetzen. In `src/solver.py` ist das die
eine Zeile `splu(...)` im Konstruktor; sie macht den Unterschied zwischen
Sekunden und Stunden.

> **Wie kommt man selbst darauf?** Zwei Standardfragen der numerischen
> Praxis: *„Was ändert sich pro Iteration wirklich?"* (nur die rechte Seite →
> Zerlegung cachen) und *„Was macht das Problem nichtlinear, und darf ich es
> pro Schritt einfrieren?"* (ja — der Fehler ist von derselben Ordnung wie
> der Zeitschrittfehler selbst).


---
## Schritt 4: Randbedingungen — hier steckt die eigentliche Physik

**Die Leitfrage:** *Die Gleichungen sind überall dieselben — was
unterscheidet dann eine Kanalströmung von einer Zylinderumströmung? Antwort:
nur die Ränder.*

Deshalb lohnt es, die Bedeutung der Randwerte wirklich zu verstehen
(Thema 4):

- **Dirichlet ($\Psi$ vorgeben):** Die Differenz der Stromfunktion zwischen
  zwei Randpunkten ist der Volumenstrom, der zwischen ihnen hindurchtritt.
  $\Psi$ auf dem Rand festzulegen heißt also: **die Normalgeschwindigkeit
  vorschreiben.** Konstantes $\Psi$ entlang eines Randstücks bedeutet „hier
  fließt nichts hindurch" — jede Wand ist eine Stromlinie. Und weil die
  Strömung inkompressibel ist, muss $\oint d\Psi = 0$ gelten: Was hineinfließt,
  muss hinaus (Kontinuitätsvoraussetzung).
- **Neumann ($\partial_n\Psi$ vorgeben):** Wegen Cauchy-Riemann ist die
  Normalableitung von $\Psi$ die **Tangentialgeschwindigkeit** — das ist die
  Haftbedingung an Wänden.
- **Periodisch:** Ränder werden miteinander identifiziert (das Gebiet wird
  topologisch ein Torus). Numerisch heißt das: zirkulante Matrizen, die
  Differenzensterne laufen per Modulo über den Rand, und der doppelte
  Endpunkt wird weggelassen.

**Das Wand-Problem — die Thom-Formel.** An einer Haftwand sind *beide*
Bedingungen gleichzeitig zu erfüllen ($\Psi$ konstant *und*
$\partial_n\Psi = 0$) — aber die Poisson-Gleichung kann nur eine davon
aufnehmen. Wohin mit der zweiten? Die Antwort findet man durch eine
Taylor-Entwicklung der Stromfunktion von der Wand aus ins Fluid:

$$
\Psi_1 = \Psi_w + h\underbrace{\partial_n\Psi_w}_{=0}
       + \frac{h^2}{2}\partial_n^2\Psi_w + \mathcal{O}(h^3)
\quad\Rightarrow\quad
\omega_w = -\nabla^2\Psi_w \approx \frac{2(\Psi_w - \Psi_1)}{h^2}
$$

Die überzählige Bedingung wird also zu einer **Vorschrift für die
Wirbelstärke an der Wand**. Das ist physikalisch tiefsinnig: Die Wand ist
genau der Ort, an dem Wirbelstärke *erzeugt* wird — die Formel beschreibt die
Quelle der ganzen späteren Wirbelstraße.

> **Wie kommt man selbst darauf?** Wenn eine Größe am Rand überbestimmt und
> eine andere unterbestimmt ist, ist das ein Hinweis, dass die überzählige
> Bedingung die fehlende *ersetzen* muss. Das Werkzeug, um solche lokalen
> Zusammenhänge herzuleiten, ist immer dasselbe wie in Schritt 2: die
> Taylor-Entwicklung.


---
## Schritt 5: Zeitintegration und Stabilität — warum der Zeitschritt nicht frei wählbar ist

**Die Leitfrage:** *Warum explodiert meine Simulation, obwohl jede einzelne
Formel richtig ist?*

Nach der Semi-Diskretisierung steht ein ODE-System
$\dot{\boldsymbol{\omega}} = f(\boldsymbol{\omega})$, und man braucht ein
Zeitschrittverfahren. Die Auswahlkriterien versteht man über zwei Analysen
aus Thema 5:

- **Das numerische Übertragungsverhalten** (modifizierte Wellenzahl):
  Zentraldifferenzen sind **dispersiv, aber nicht dissipativ** — sie
  verfälschen Phasen, fressen aber keine Amplitude. Einseitige
  (Aufwind-)Differenzen 1. Ordnung sind dissipativ: Sie wirken wie künstliche
  Viskosität („numerische Diffusion").
- **Die von-Neumann-Stabilitätsanalyse:** Man setzt eine Fourier-Mode in das
  volldiskretisierte Schema ein und prüft, ob ihr Verstärkungsfaktor
  betragsmäßig ≤ 1 bleibt. Daraus fallen die beiden zentralen Grenzen heraus:
  - **CFL-Bedingung** (Konvektion): $|u|\,\Delta t \lesssim h$ — die
    Information darf pro Zeitschritt höchstens etwa eine Zelle weit laufen.
  - **Diffusionsbedingung:** $\nu\,\Delta t \lesssim h^2/4$ — beachte das
    *Quadrat*: Gitterverfeinerung verteuert explizite Verfahren doppelt.

Eine klassische Falle, die man mit diesen Werkzeugen sofort durchschaut: Das
**explizite Euler-Verfahren mit Zentraldifferenzen ist für reine Konvektion
immer instabil** — sein Stabilitätsgebiet enthält kein Stück der imaginären
Achse, auf der die Eigenwerte der Zentraldifferenzen-Konvektion liegen.
Klassisches RK4 dagegen schneidet die imaginäre Achse — deshalb ist es hier
das Standardverfahren. (Mit Viskosität rutschen die Eigenwerte in die linke
Halbebene, dann geht Euler mit kleinem $\Delta t$ gerade noch — aber „geht
gerade noch" ist keine Verfahrensgrundlage.)

> **Merksatz:** Räumliches Schema und Zeitschrittverfahren sind ein *Paar*.
> Stabilität ist keine Eigenschaft des einen oder des anderen, sondern der
> Kombination.


---
## Schritt 6: Validierung — erst beweisen, dann bauen

**Die Leitfrage:** *Woher weiß ich, dass mein Löser rechnet, was er soll —
bevor ich ihn auf ein Problem loslasse, dessen Lösung niemand kennt?*

**Warum der Taylor-Green-Wirbel das perfekte Testproblem ist:** Er ist einer
der ganz wenigen Fälle, in denen die *nichtlinearen* Gleichungen eine exakte
Lösung haben. Der Grund ist eine gezielte Annahme (Thema 6): Wenn die
Wirbelstärke überall ein festes Vielfaches der Stromfunktion ist
($\omega_z = \eta\Psi$, Beltrami-Strömung), dann steht der Gradient von
$\omega_z$ senkrecht auf $\boldsymbol{u}$ — **der Konvektionsterm, die einzige
Nichtlinearität, verschwindet identisch.** Übrig bleibt lineare Diffusion:

$$
\Psi = \sin x \sin y\, \mathrm{e}^{-2\nu t}, \qquad \omega_z = 2\Psi
$$

(Der Ansatz $\Psi \sim \sin x \sin y$ in der Poisson-Gleichung liefert
$\eta = 2$, der Ansatz $\Psi \sim \mathrm{e}^{-\nu\eta t}$ im
Diffusionsrest den Zeitabfall.)

**Was man damit prüft — drei Ebenen:**

1. **Konvergenzordnung messen, nicht behaupten:** Fehler gegen die analytische
   Lösung bei 16/32/64 Gitterpunkten; der Logarithmus der Fehlerverhältnisse
   muss die nominelle Ordnung ergeben (hier: 1,98 und 2,00 — Punktlandung
   auf Ordnung 2). Eine falsche Randbehandlung würde genau hier auffliegen.
2. **Erhaltungsgrößen als Diffusionsdetektor:** Die Enstrophie
   $\mathcal{E} = \tfrac12\int\omega_z^2\,dA$ muss im reibungsfreien Fall
   konstant bleiben. Jede Abweichung ist *numerische* Diffusion — dieser eine
   Zahlenwert misst die Qualität des Schemas. (RK4 + Zentraldifferenzen:
   Drift exakt 0.)
3. **Grenzfälle:** $\nu = 0$ (stationär), $\nu > 0$ (exponentieller Abfall
   mit bekannter Rate $4\nu$).

> **Die wichtigste Regel des ganzen Projekts:** Niemals Numerik und Geometrie
> gleichzeitig debuggen. Der Taylor-Green-Wirbel validiert die komplette
> Numerik auf dem einfachsten Gebiet (periodisches Quadrat). Erst wenn hier
> alles bewiesen ist, kommt die schwierige Geometrie dazu — dann kann jeder
> neue Fehler nur noch von der Geometrie stammen.


---
## Schritt 7: Die Zylinderumströmung — Geometrie als Koordinatenproblem

**Die Leitfrage:** *Der Löser kann Rechtecke. Der Zylinder ist rund. Wer von
beiden passt sich an?*

**7a — Die Koordinatenwahl.** Auf einem kartesischen Gitter liegt die
Zylinderwand quer zu allen Gitterlinien — Randbedingungen werden dann
schmutzig (Treppenstufen). Die elegante Antwort: **Koordinaten wählen, in
denen der Rand eine Koordinatenlinie ist.** In Polarkoordinaten $(r,\theta)$
ist der Zylinder einfach die Linie $r = r_0$, und das Rechengebiet
$[r_0, r_\infty] \times [0, 2\pi)$ ist wieder ein Rechteck — der komplette
validierte Apparat aus den Schritten 2–6 bleibt verwendbar. Die
$\theta$-Richtung ist dabei von Natur aus periodisch (Schritt 4 liefert die
Matrizen dafür).

**7b — Die logarithmische Verfeinerung.** Physikalisch spielt die Musik in
der dünnen Grenzschicht an der Wand; weit draußen passiert wenig. Die
Substitution $r_{\ln} = \ln r$ verdichtet das Gitter automatisch an der Wand
(gleiche $r_{\ln}$-Schritte sind nahe der Wand kleine, außen große
$r$-Schritte) — und als Bonus wird der Laplace-Operator *einfacher*:

$$
\nabla^2 = \mathrm{e}^{-2r_{\ln}}\left(\partial_{r_{\ln}}^2 + \partial_\theta^2\right)
$$

Der Kern ist derselbe wie im kartesischen Fall, nur mit einem metrischen
Vorfaktor. (Herleitung über die Kettenregel: $\partial_{r_{\ln}} =
r\,\partial_r$ — Thema 6.)

**7c — Entdimensionalisierung.** Statt $r_0$, $u_\infty$ und $\nu$ einzeln zu
wählen, skaliert man alle Größen mit $r_0$ und $u_\infty$ (Buckingham,
Thema 6). Übrig bleibt **ein einziger Parameter**:

$$
\mathrm{Re} = \frac{2 r_0 u_\infty}{\nu}
$$

Das ist nicht Kosmetik: Es macht Ergebnisse vergleichbar (mit Literatur, mit
anderen Codes) und reduziert den Parameterraum von drei Dimensionen auf eine.

**7d — Die Randbedingungen zusammensetzen.** Jetzt zahlt sich Schritt 4 aus,
denn alles ist schon da:
- **Wand** ($r_{\ln} = 0$): $\Psi = 0$ (Stromlinie) + Thom-Formel für
  $\omega_w$ (Haftung).
- **Außenrand:** die bekannte Potentialströmung um den Zylinder
  ($\Psi = \sin\theta\,(r - 1/r)$) als Fernfeld, $\omega_z = 0$.
- **$\theta$:** periodisch.

**7e — Die Symmetrie brechen.** Eine Überraschung mit Lerneffekt: Die
Wirbelstraße kommt von allein *nicht*. Gleichungen, Gitter und
Randbedingungen sind perfekt symmetrisch, und anders als im Experiment gibt
es kein Rauschen, das die Instabilität anstößt — die Simulation verharrt in
der (instabilen!) symmetrischen Lösung. Man muss die Symmetrie gezielt
stören, z. B. durch eine kurze Drehung des Zylinders oder ein kurzes Kippen
der Anströmung. Danach schaukelt sich die Instabilität von selbst zur
periodischen Ablösung auf.

**7f — Validieren ohne analytische Lösung.** Für den Zylinder gibt es keine
exakte Lösung — aber Jahrzehnte an Messwerten. Der Standard-Check ist die
**Strouhal-Zahl** (dimensionslose Ablösefrequenz aus einem
Geschwindigkeitssignal im Nachlauf): Bei $\mathrm{Re} = 100$ muss
$\mathrm{St} \approx 0{,}16$–$0{,}17$ herauskommen. Unser Löser liefert
0,168 — das ist der Beleg, dass die ganze Kette stimmt.

> **Wie kommt man selbst darauf?** Die Schlüsselidee ist 7a: *Nicht der Löser
> muss die Geometrie lernen, sondern die Koordinaten müssen sich der
> Geometrie anpassen.* Danach ist der Zylinder „nur" eine Wiederverwendung
> aller vorherigen Schritte auf einem transformierten Rechteck.


---
## Schritt 8: Sichtbar machen — der Wechsel in die Lagrange-Perspektive

Die Felder $\omega_z$ und $\Psi$ sind Euler'sche Größen (feste Orte). Was man
im Experiment *sieht* (Rauch, Farbstoff), sind aber mitschwimmende Teilchen —
die Lagrange-Perspektive (Thema 1). Für die Visualisierung wechselt man also
den Blickwinkel: masselose Partikel werden mit dem interpolierten
Geschwindigkeitsfeld integriert,

$$
\boldsymbol{p}(t+\Delta t) \approx \boldsymbol{p}(t) + \Delta t\,
\boldsymbol{u}(\boldsymbol{p}(t), t)
$$

Die Steigerung davon sind **Lagrange-kohärente Strukturen** (FTLE): Man
verfolgt ein ganzes Partikelgitter, misst über die Jacobimatrix der
Flussabbildung, wie stark benachbarte Partikel auseinandergezogen werden
(größter Eigenwert des Deformationstensors $J^\top J$), und erhält die
verborgenen Transportbarrieren der Strömung — vorwärts in der Zeit die
instabilen (rot), rückwärts die stabilen Strukturen (blau).


---
## Der Bauplan zum Selbermachen

Die Reihenfolge ist der halbe Erfolg. Jede Stufe hat ein hartes
Fertig-Kriterium — erst wenn es erfüllt ist, geht es weiter:

| # | Stufe | Fertig-Kriterium |
|---|-------|------------------|
| 1 | 1D-Ableitungsmatrizen (zentral, einseitig, periodisch) | Konvergenzordnung an $\sin x$ gemessen = nominelle Ordnung |
| 2 | 2D-Matrizen per Kronecker-Produkt, sparse | $\nabla^2 \sin x \sin y = -2\sin x\sin y$ auf dem Gitter |
| 3 | Poisson-Löser mit Dirichlet-Maske + gecachter LR-Zerlegung | bekannte $(\omega, \Psi)$-Paare werden reproduziert |
| 4 | Dreischritt-Schleife (Poisson → Cauchy-Riemann → Transport) mit RK4 | läuft stabil bei eingehaltener CFL-/Diffusionsgrenze |
| 5 | Taylor-Green-Validierung | Ordnung ≈ 2 gemessen; Enstrophie-Drift ≈ 0 bei $\nu=0$ |
| 6 | Log-Polargitter + Metrikfaktoren | Potentialströmung ($\omega=0$) wird stationär gehalten |
| 7 | Wand- und Fernfeld-Randbedingungen, Thom-Formel | Grenzschicht baut sich auf; stationäres Wirbelpaar bei Re ≈ 20–40 |
| 8 | Störung + lange Laufzeit bei Re ≈ 100 | periodische Ablösung; $\mathrm{St} \approx 0{,}16$ |
| 9 | Partikel/FTLE-Visualisierung | Strukturen wie in Thema 6 sichtbar |

Zuordnung zum Programmcode: Stufen 1–2 → `src/operators.py`, 3–4 →
`src/solver.py`, 5 → `src/taylor_green.py`, 6–8 → `src/cylinder.py`,
9 → `src/visualization.py`; die Messungen dazu in `tests/` und
`src/benchmark.py`.


---
## Die Meta-Prinzipien (das eigentlich Übertragbare)

Wer diese sieben Denkmuster mitnimmt, kann nicht nur diesen Löser
nachbauen, sondern auch das nächste, ganz andere Problem angehen:

1. **Vereinfache zuerst das Problem, nicht den Code.** (2D, inkompressibel,
   laminar — Schritt 0)
2. **Eliminiere Störendes durch Operatoren, erfülle Zwänge durch
   Konstruktion.** (Rotation gegen den Druck, Stromfunktion gegen die
   Kontinuität — Schritt 1)
3. **Mache Ableitungen zu Matrizen — dann ist Analysis nur noch lineare
   Algebra.** (Schritt 2)
4. **Zerlege Nichtlineares in eine Folge linearer Teilprobleme und cache,
   was sich nicht ändert.** (Dreischritt-Schleife, LR-Zerlegung — Schritt 3)
5. **Die Physik eines Problems steckt fast vollständig in den
   Randbedingungen.** (Thom-Formel: die Wand als Wirbelquelle — Schritt 4)
6. **Validiere am analytischen Grenzfall, bevor du Geometrie hinzufügst —
   und miss Konvergenzordnungen, statt sie zu glauben.** (Taylor-Green —
   Schritt 6)
7. **Passe die Koordinaten der Geometrie an, nicht den Löser.**
   (Log-Polar — Schritt 7)
