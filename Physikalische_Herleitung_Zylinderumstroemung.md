# Die Zylinderumströmung — vom physikalischen Grundgesetz zur Kármán'schen Wirbelstraße


---

## Inhalt

1. [Das physikalische Problem](#1-das-physikalische-problem)
2. [Erhaltungssätze: Kontinuität und Impuls](#2-erhaltungssätze-kontinuität-und-impuls)
3. [Die allgemeine Transportgleichung als Blaupause](#3-die-allgemeine-transportgleichung-als-blaupause)
4. [Formulierungswechsel: Wirbelstärke und Stromfunktion](#4-formulierungswechsel-wirbelstärke-und-stromfunktion)
5. [Das geschlossene Gleichungssystem](#5-das-geschlossene-gleichungssystem)
6. [Entdimensionalisierung: die Reynolds-Zahl](#6-entdimensionalisierung-die-reynolds-zahl)
7. [Koordinatenwahl: logarithmische Polarkoordinaten](#7-koordinatenwahl-logarithmische-polarkoordinaten)
8. [Randbedingungen — hier steckt die Physik](#8-randbedingungen--hier-steckt-die-physik)
9. [Diskretisierung: von der Analysis zur linearen Algebra](#9-diskretisierung-von-der-analysis-zur-linearen-algebra)
10. [Der Lösungsalgorithmus](#10-der-lösungsalgorithmus)
11. [Zeitintegration und Stabilität](#11-zeitintegration-und-stabilität)
12. [Symmetriebrechung: warum die Wirbelstraße nicht von allein kommt](#12-symmetriebrechung-warum-die-wirbelstraße-nicht-von-allein-kommt)
13. [Physikalische Auswertung](#13-physikalische-auswertung)
14. [Validierung](#14-validierung)
15. [Die durchgehende Kausalkette im Überblick](#15-die-durchgehende-kausalkette-im-überblick)

---

## 1. Das physikalische Problem

Ein starrer Kreiszylinder mit Radius $r_0$ steht senkrecht in einer Strömung, die
weit stromauf die gleichförmige Geschwindigkeit $u_\infty$ besitzt. Das Fluid ist
newtonsch, inkompressibel und isotherm mit konstanter Dichte $\rho$ und kinematischer
Viskosität $\nu$.

**Die zu erklärende Beobachtung:** Oberhalb einer kritischen Anströmgeschwindigkeit
löst die Strömung nicht mehr symmetrisch ab, sondern es lösen sich abwechselnd von
Ober- und Unterseite Wirbel entgegengesetzten Drehsinns — die *Kármán'sche
Wirbelstraße*. Diese Instabilität ist kein Rechenartefakt, sondern eine echte
Eigenschaft der Navier-Stokes-Gleichungen: die symmetrische Lösung existiert
weiterhin, ist aber instabil.

**Vereinfachende Annahmen und ihre Rechtfertigung:**

| Annahme | Begründung |
|---|---|
| Zweidimensional | Für $\mathrm{Re} \lesssim 190$ ist der Nachlauf des Kreiszylinders experimentell nachweislich spannweitig kohärent — 3D-Instabilitäten (Mode A/B) setzen erst darüber ein. |
| Inkompressibel | Machzahl $\ll 0{,}3$; Dichteänderungen sind vernachlässigbar. |
| Newtonsch | Linearer Zusammenhang zwischen Spannung und Deformationsgeschwindigkeit (Luft, Wasser). |
| Isotherm | Keine Auftriebs- oder Temperaturkopplung; $\rho, \nu$ konstant. |
| Laminar | Bei $\mathrm{Re} \approx 100$ ist die Strömung noch nicht turbulent; kein Turbulenzmodell nötig. |

Diese Vereinfachungen betreffen das *Problem*, nicht die Numerik. Das ist eine
methodische Grundregel: erst das Problem so weit reduzieren, dass es genau die
gesuchte Physik noch enthält, dann erst rechnen.

---

## 2. Erhaltungssätze: Kontinuität und Impuls

### 2.1 Massenerhaltung

Bilanziert man den Massenstrom über ein infinitesimales, ortsfestes Kontrollvolumen
(Euler'sche Betrachtungsweise), so ergibt sich mit der Taylor-Entwicklung der
austretenden Ströme und dem Grenzübergang $V \to 0$ die Kontinuitätsgleichung:

$$
\partial_t \rho + \boldsymbol{\nabla}\cdot(\rho\boldsymbol{u}) = 0
$$

Unter der Annahme der Inkompressibilität ($\rho = \text{konst.}$) reduziert sie sich auf
eine reine **Zwangsbedingung an das Geschwindigkeitsfeld**:

$$
\boxed{\;\boldsymbol{\nabla}\cdot\boldsymbol{u} = 0\;}
$$

Diese Gleichung beschreibt keine Entwicklung in der Zeit — sie sagt nicht, *wie* sich
etwas verhält, sondern *was verboten ist*. Dieser strukturelle Unterschied wird in
Abschnitt 4 entscheidend.

### 2.2 Impulserhaltung — die Navier-Stokes-Gleichung

Die Impulsbilanz am selben Kontrollvolumen entspricht dem 2. Newton'schen Gesetz.
Weil Newtons Kraftbegriff einer *Punktmasse* — also der Lagrange'schen Sicht —
entspringt, muss die Geschwindigkeit substantiell abgeleitet werden, bevor im
Euler'schen Sinne bilanziert werden kann:

$$
\boldsymbol{F} = m\boldsymbol{a} = (\rho V)(D_t\boldsymbol{u})
= \left[\boldsymbol{F}_\text{ein} - \boldsymbol{F}_\text{aus}\right] + \boldsymbol{F}_\text{quell}
$$

Mit den Oberflächenkräften (Druck zentrisch, Spannung exzentrisch) und den
Volumenkräften im Quellterm folgt

$$
D_t\boldsymbol{u} = -\frac{\boldsymbol{\nabla}p}{\rho} + \frac{\boldsymbol{\nabla}\cdot\boldsymbol{\tau}}{\rho} + \boldsymbol{g}
$$

Der Stokes'sche Spannungsansatz für newtonsche Fluide

$$
\boldsymbol{\tau} = \mu\left[\boldsymbol{\nabla}\boldsymbol{u} + (\boldsymbol{\nabla}\boldsymbol{u})^\top\right]
$$

liefert zusammen mit den Identitäten

$$
\boldsymbol{\nabla}\cdot(\boldsymbol{\nabla}\boldsymbol{u}) = \nabla^2\boldsymbol{u},
\qquad
\boldsymbol{\nabla}\cdot(\boldsymbol{\nabla}\boldsymbol{u})^\top
= \boldsymbol{\nabla}(\underbrace{\boldsymbol{\nabla}\cdot\boldsymbol{u}}_{=\,0}) = 0,
\qquad \mu = \nu\rho
$$

und der Zerlegung der substantiellen Ableitung

$$
D_t\boldsymbol{u} = \partial_t\boldsymbol{u} + (\boldsymbol{u}\cdot\boldsymbol{\nabla})\boldsymbol{u}
$$

die inkompressible Navier-Stokes-Gleichung:

$$
\boxed{\;\partial_t\boldsymbol{u} + (\boldsymbol{u}\cdot\boldsymbol{\nabla})\boldsymbol{u}
= -\frac{\boldsymbol{\nabla}p}{\rho} + \nu\nabla^2\boldsymbol{u} + \boldsymbol{g}\;}
$$

Damit ist das physikalische Modell vollständig: zwei Gleichungen (Kontinuität,
Impuls) für drei Unbekannte ($u$, $v$, $p$) in 2D.

---

## 3. Die allgemeine Transportgleichung als Blaupause

Beide Herleitungen folgen demselben Muster, und dieses Muster lässt sich
verallgemeinern. Für eine beliebige Feldgröße $\Phi$ mit Quellterm $R$ lautet die
Bilanz

$$
\partial_t\Phi + \boldsymbol{\nabla}\cdot\boldsymbol{J} = R
$$

wobei sich die Stromdichte aus einem Diffusions- und einem Konvektionsanteil
zusammensetzt:

$$
\boldsymbol{J} = \underbrace{-c\,\boldsymbol{\nabla}\Phi}_{\text{1. Fick'sches Gesetz}} + \underbrace{\boldsymbol{u}\Phi}_{\text{Konvektion}}
$$

Deren Divergenz ergibt bei konstantem $c$ und mit der Produktregel

$$
\boldsymbol{\nabla}\cdot\boldsymbol{J} = -c\nabla^2\Phi + \underbrace{(\boldsymbol{\nabla}\cdot\boldsymbol{u})\Phi}_{=\,0\ \text{(inkompressibel)}} + (\boldsymbol{u}\cdot\boldsymbol{\nabla})\Phi
$$

und damit die **Konvektions-Diffusions-Gleichung**:

$$
\boxed{\;
\underbrace{\partial_t\Phi}_{\text{Instationarität}}
= \underbrace{c\nabla^2\Phi}_{\text{Diffusion}} - \underbrace{(\boldsymbol{u}\cdot\boldsymbol{\nabla})\Phi}_{\text{Konvektion}} + \underbrace{R}_{\text{Quelle}}
\;}
$$

Das ist mehr als eine Formalie: Sobald die Zielgleichung diese Struktur hat, ist
bekannt, welche numerischen Werkzeuge greifen — Aufwind für die Konvektion,
Diffusionsstabilitätsgrenze für den Laplace-Term, CFL-Bedingung für die Zeitschrittweite.
Die eindimensionale Version dieser Gleichung ist deshalb der sinnvolle Testfall,
bevor die zweidimensionale Simulation überhaupt beginnt.

---

## 4. Formulierungswechsel: Wirbelstärke und Stromfunktion

Die primitiven Variablen $(u, v, p)$ direkt zu lösen, bringt zwei strukturelle
Schwierigkeiten mit sich:

1. **Der Druck hat keine eigene Zeitentwicklungsgleichung.** Er wird nicht
   transportiert, sondern stellt sich instantan so ein, dass die Kontinuität erhalten
   bleibt — er müsste in jedem Schritt implizit mitbestimmt werden.
2. **Die Kontinuitätsgleichung ist eine Zwangsbedingung**, keine
   Bestimmungsgleichung (siehe 2.1).

Für beide Probleme gibt es je einen klassischen Handgriff.

### 4.1 Handgriff 1 — Störendes durch einen Operator eliminieren

Der Druck erscheint ausschließlich als Gradient, und Gradientenfelder sind wirbelfrei:
$\boldsymbol{\nabla}\times\boldsymbol{\nabla}p = \boldsymbol{0}$. Wendet man also die
Rotation auf die Impulsgleichung an, verschwindet der Druck vollständig — ebenso die
Schwerkraft, denn auch sie ist konservativ ($\boldsymbol{F}_\text{grav} = m\boldsymbol{g}
= -\boldsymbol{\nabla}W_\text{pot}$).

Mit der **Wirbelstärke** $\boldsymbol{\omega} \coloneqq \boldsymbol{\nabla}\times\boldsymbol{u}$
und den Identitäten

$$
\begin{aligned}
(\boldsymbol{u}\cdot\boldsymbol{\nabla})\boldsymbol{u}
&= \tfrac{1}{2}\boldsymbol{\nabla}\boldsymbol{u}^2 - \boldsymbol{u}\times\boldsymbol{\omega} \\
\boldsymbol{\nabla}\times(\boldsymbol{u}\times\boldsymbol{\omega})
&= \boldsymbol{u}(\underbrace{\boldsymbol{\nabla}\cdot\boldsymbol{\omega}}_{=\,0}) - \boldsymbol{\omega}(\underbrace{\boldsymbol{\nabla}\cdot\boldsymbol{u}}_{=\,0}) + (\boldsymbol{\omega}\cdot\boldsymbol{\nabla})\boldsymbol{u} - (\boldsymbol{u}\cdot\boldsymbol{\nabla})\boldsymbol{\omega}
\end{aligned}
$$

ergibt sich die **Wirbeltransportgleichung**

$$
\partial_t\boldsymbol{\omega} + (\boldsymbol{u}\cdot\boldsymbol{\nabla})\boldsymbol{\omega}
= \underbrace{(\boldsymbol{\omega}\cdot\boldsymbol{\nabla})\boldsymbol{u}}_{\text{Wirbelstreckung}} + \nu\nabla^2\boldsymbol{\omega}
$$

**Der 2D-Sonderfall.** Mit $u_z = 0$ und $\partial_z = 0$ hat die Wirbelstärke nur
eine Komponente,

$$
\omega_z = \partial_x u_y - \partial_y u_x,
$$

die senkrecht auf der Strömungsebene steht. Damit gilt
$(\boldsymbol{\omega}\cdot\boldsymbol{\nabla})\boldsymbol{u} = \omega_z \partial_z \boldsymbol{u} = \boldsymbol{0}$ —
**der Wirbelstreckungsterm entfällt**. Übrig bleibt eine skalare Gleichung:

$$
\boxed{\;\partial_t\omega_z + (\boldsymbol{u}\cdot\boldsymbol{\nabla})\omega_z
= \nu\nabla^2\omega_z\;}
$$

Dies ist exakt die allgemeine Transportgleichung aus Abschnitt 3 mit $\Phi = \omega_z$,
$c = \nu$, $R = 0$. Physikalisch: **Wirbelstärke wird mit der Strömung mitgeführt und
diffundiert viskos, wird aber im Inneren des Fluids weder erzeugt noch vernichtet.**
Die einzige Quelle liegt am Rand (Abschnitt 8) — das ist der Kern des ganzen Problems.

Dass die Wirbelstreckung in 2D wegfällt, ist zugleich der Grund, warum die
2D-Turbulenz sich qualitativ von der 3D-Turbulenz unterscheidet (inverse
Energiekaskade statt Kolmogorov-Kaskade) — und warum die 2D-Simulation überhaupt mit
moderatem Aufwand rechenbar ist.

### 4.2 Handgriff 2 — Zwangsbedingungen durch Konstruktion erfüllen

Statt $\boldsymbol{\nabla}\cdot\boldsymbol{u} = 0$ in jedem Schritt zu erzwingen, wählt
man Variablen, die sie automatisch erfüllen.

**Stromlinien** sind Kurven, deren Tangenten kollinear zum Geschwindigkeitsfeld
verlaufen: $d\boldsymbol{x}\times\boldsymbol{u} \overset{!}{=} \boldsymbol{0}$, in 2D
also $dx/u = dy/v$.

Die **Stromfunktion** $\Psi$ wird so definiert, dass die Stromlinien ihre Niveaulinien
sind und die Differenz zweier Niveaus dem dazwischen hindurchtretenden Volumenstrom
entspricht:

$$
\Psi(B) - \Psi(A) \coloneqq \int_A^B \underbrace{(d\dot V_x - d\dot V_y)}_{d\dot V_\perp}
= \int_A^B (u\,dy - v\,dx)
$$

Für ein infinitesimales Wegelement entfällt das Integral, $d\Psi = u\,dy - v\,dx$.
Ein Vergleich mit dem exakten Ortsdifferential
$d\Psi = (\partial_x\Psi)dx + (\partial_y\Psi)dy$ liefert die
**Cauchy-Riemann-Gleichungen**:

$$
\boxed{\;u = \partial_y\Psi, \qquad v = -\partial_x\Psi\;}
$$

Damit ist die Kontinuität automatisch erfüllt — nach dem Satz von Schwarz:

$$
\boldsymbol{\nabla}\cdot\boldsymbol{u} = \partial_x u + \partial_y v
= \partial_x\partial_y\Psi - \partial_y\partial_x\Psi \equiv 0
$$

Die Zwangsbedingung ist keine Bedingung mehr, sondern eine Eigenschaft der
Konstruktion. Nebenbei reduziert die Stromfunktion das Vektorfeld
$\boldsymbol{u}$ auf ein Skalarfeld.

### 4.3 Die Verbindung — die Poisson-Gleichung

Setzt man den Stromfunktionsansatz in die Definition der Wirbelstärke ein:

$$
\omega_z = \partial_x v - \partial_y u
= -\partial_x(\partial_x\Psi) - \partial_y(\partial_y\Psi)
= -(\partial_x^2 + \partial_y^2)\Psi
$$

folgt die **Poisson-Gleichung**

$$
\boxed{\;\omega_z = -\nabla^2\Psi\;}
$$

Sie schließt das System.

---

## 5. Das geschlossene Gleichungssystem

Zwei Skalarfelder, drei Gleichungen, kein Druck, keine Zwangsbedingung:

$$
\begin{aligned}
\text{(I) Poisson:} \qquad & \omega_z = -\nabla^2\Psi \\[4pt]
\text{(II) Cauchy-Riemann:} \qquad & \boldsymbol{u} =
\begin{bmatrix} \partial_y\Psi \\ -\partial_x\Psi \end{bmatrix} \\[4pt]
\text{(III) Transport:} \qquad & \partial_t\omega_z
= \nu\nabla^2\omega_z - (\boldsymbol{u}\cdot\boldsymbol{\nabla})\omega_z
\end{aligned}
$$

Die Rollenverteilung ist dabei nicht symmetrisch und für den Algorithmus zentral:
**$\omega_z$ ist der eigentliche Zustand des Systems** — es ist die einzige Größe mit
einer Zeitableitung. $\Psi$ und $\boldsymbol{u}$ sind daraus abgeleitete
Momentaufnahmen, die instantan mitfolgen. Diese Kausalkette bestimmt die Reihenfolge
im Lösungsalgorithmus (Abschnitt 10).

### Nebenprodukt: die Enstrophie

Aus der zeitlichen Ableitung der kinetischen Energie folgt mit der
Navier-Stokes-Gleichung und partieller Integration ein bemerkenswerter Zusammenhang:

$$
\partial_t \underbrace{\left(\frac{\rho}{2}\int_\Omega \boldsymbol{u}^2\,dV\right)}_{E_\text{kin}}
= -\mu \underbrace{\left(\int_\Omega \boldsymbol{\omega}^2\,dV\right)}_{\text{Enstrophie }\mathcal{E}}
$$

Die kinetische Energie wird also **über die Wirbelstärke und die Viskosität** in Wärme
dissipiert. Für $\nu = 0$ ist die Enstrophie eine Erhaltungsgröße — was sie zum
idealen Detektor für *numerische* Diffusion macht (Abschnitt 14).

---

## 6. Entdimensionalisierung: die Reynolds-Zahl

Das System enthält neun dimensionsbehaftete Größen:

$$
t,\ \omega_z,\ \Psi,\ x,\ y,\ u,\ v,\ r_0,\ u_\infty
$$

Nach dem Buckingham'schen $\Pi$-Theorem reduziert die Wahl von zwei
Bezugsgrößen die Anzahl unabhängiger Parameter entsprechend. Zulässig ist jede Wahl,
deren Einheiten zusammen alle vorkommenden Einheiten aufspannen. Mit $r_0$ und
$u_\infty$:

$$
[x] = [y] = [r_0], \quad [u] = [v] = [u_\infty], \quad
[\Psi] = [r_0][u_\infty], \quad [\omega_z] = [u_\infty]/[r_0]
$$

Das genügt. Skaliert man jede Größe entsprechend um, bleibt genau **ein** Parameter
übrig — die Reynolds-Zahl:

$$
\boxed{\;\mathrm{Re} = \frac{2r_0 u_\infty}{\nu} = \frac{D\,u_\infty}{\nu}\;}
$$

und das System lautet

$$
\begin{aligned}
\omega_z &= -\nabla^2\Psi \\
\boldsymbol{u} &= \begin{bmatrix} \partial_y\Psi \\ -\partial_x\Psi \end{bmatrix} \\
\partial_t\omega_z &= \left(\frac{2}{\mathrm{Re}}\nabla^2 - \boldsymbol{u}\cdot\boldsymbol{\nabla}\right)\omega_z
\end{aligned}
$$

Der Faktor 2 entsteht dadurch, dass mit dem *Radius* skaliert, die Reynolds-Zahl aber
über den *Durchmesser* definiert wird. Skaliert man stattdessen mit $D$ — was bei
kartesischen Gittern naheliegt, wie in `src/immersed.py` —, steht dort $1/\mathrm{Re}$.

**Warum das mehr als Kosmetik ist:** Die Reynolds-Zahl ist das Verhältnis von
Trägheits- zu Zähigkeitskräften. Zwei Strömungen mit gleichem $\mathrm{Re}$ sind
dynamisch ähnlich, unabhängig von der absoluten Größe. Der Parameterraum schrumpft
von drei Dimensionen auf eine, und die Ergebnisse werden direkt mit Literaturwerten
und anderen Codes vergleichbar.

---

## 7. Koordinatenwahl: logarithmische Polarkoordinaten

**Die Leitfrage:** Der Löser arbeitet auf Rechtecken. Der Zylinder ist rund. Wer von
beiden passt sich an?

### 7.1 Polarkoordinaten

Auf einem kartesischen Gitter läuft die Zylinderwand quer zu allen Gitterlinien; die
Randbedingungen werden dann treppenförmig approximiert. Die elegantere Antwort:
**Koordinaten wählen, in denen der Rand eine Koordinatenlinie ist.** In
Polarkoordinaten $(r,\theta)$ ist der Zylinder die Linie $r = r_0$, und das
Rechengebiet $[r_0, r_\infty] \times [0, 2\pi)$ ist wieder ein Rechteck. Die
$\theta$-Richtung ist von Natur aus periodisch.

Mit

$$
\boldsymbol{\nabla} = \frac{1}{r}\begin{bmatrix} r\partial_r \\ \partial_\theta \end{bmatrix},
\qquad
\nabla^2 = \frac{1}{r^2}\left[(r\partial_r)(r\partial_r) + \partial_\theta^2\right]
$$

lautet das System

$$
\omega_z = -\nabla^2\Psi, \qquad
\boldsymbol{u} = \frac{1}{r}\begin{bmatrix} \partial_\theta\Psi \\ -r\partial_r\Psi \end{bmatrix}
$$

### 7.2 Die logarithmische Verfeinerung

Physikalisch spielt sich das Wesentliche in der dünnen Grenzschicht an der Wand ab;
weit draußen ist die Strömung nahezu potentialförmig. Ein Gitter mit gleichmäßigem
$\Delta r$ würde Rechenpunkte dort verschwenden, wo nichts passiert.

Die Substitution

$$
r_{\ln} \coloneqq \ln r
$$

löst beide Probleme auf einmal. Aus

$$
\frac{dr}{dr_{\ln}} = \left(\frac{d}{dr}\ln r\right)^{-1} = r
\qquad\Longrightarrow\qquad
\partial_{r_{\ln}} = r\,\partial_r
$$

folgt

$$
\boxed{\;
\boldsymbol{\nabla} = \mathrm{e}^{-r_{\ln}}\begin{bmatrix}\partial_{r_{\ln}}\\\partial_\theta\end{bmatrix},
\qquad
\nabla^2 = \mathrm{e}^{-2r_{\ln}}\left(\partial_{r_{\ln}}^2 + \partial_\theta^2\right)
\;}
$$

Zwei Vorteile gleichzeitig:

- **Automatische Grenzschichtauflösung:** Gleiche Schritte in $r_{\ln}$ bedeuten
  kleine $\Delta r$ nahe der Wand und große weit außen — genau die gewünschte
  Verteilung, ohne manuelles Gitterstretching.
- **Einfacherer Operator:** Der Laplace-Operator ist wieder der kartesische Kern, nur
  mit einem skalaren metrischen Vorfaktor $\mathrm{e}^{-2r_{\ln}}$. Alle bereits
  validierten Ableitungsmatrizen bleiben unverändert nutzbar.

### 7.3 Das Endsystem der Simulation

Mit $r_0 = 1$ liegt die Wand bei $r_{\ln} = 0$.

$$
\begin{aligned}
\omega &= -\mathrm{e}^{-2r_{\ln}}\left(\partial_\theta^2 + \partial_{r_{\ln}}^2\right)\Psi \\
u_r &= \mathrm{e}^{-r_{\ln}}\,\partial_\theta\Psi,
\qquad u_\theta = -\mathrm{e}^{-r_{\ln}}\,\partial_{r_{\ln}}\Psi \\
\dot\omega &= \mathrm{e}^{-2r_{\ln}}\left[\frac{2}{\mathrm{Re}}
\left(\partial_\theta^2 + \partial_{r_{\ln}}^2\right)\omega - \partial_\theta\Psi\,\partial_{r_{\ln}}\omega + \partial_{r_{\ln}}\Psi\,\partial_\theta\omega\right]
\end{aligned}
$$

Der Konvektionsterm ist hier ausgeschrieben die Jacobi-Determinante
$J(\Psi, \omega)$ — Ausdruck dessen, dass Wirbelstärke entlang von Stromlinien
transportiert wird.

**Rechengebiet der Implementierung:** $128 \times 96$ Punkte in
$(\theta, r_{\ln})$, $r_\text{max} = 30$. Der Außenrand muss weit genug entfernt
sein, damit die dort angesetzte Fernfeldbedingung den Nachlauf nicht künstlich
beeinflusst.

---

## 8. Randbedingungen — hier steckt die Physik

Die Gleichungen sind überall dieselben. Was eine Kanalströmung von einer
Zylinderumströmung unterscheidet, sind **ausschließlich die Ränder**.

### 8.1 Bedeutung der Randtypen

Die Stromfunktionsformulierung gibt den Randbedingungen eine direkte physikalische
Lesart:

- **Dirichlet ($\Psi$ vorgegeben):** Die Differenz von $\Psi$ zwischen zwei
  Randpunkten ist der dazwischen hindurchtretende Volumenstrom. $\Psi$ am Rand
  festzulegen heißt also, die **Normalgeschwindigkeit** vorzuschreiben. Konstantes
  $\Psi$ entlang eines Randstücks bedeutet „hier fließt nichts hindurch" — jede
  undurchlässige Wand ist eine Stromlinie. Aus der Kontinuität folgt zwingend
  $\oint d\Psi = 0$.
- **Neumann ($\partial_n\Psi$ vorgegeben):** Wegen Cauchy-Riemann ist die
  Normalableitung von $\Psi$ die **Tangentialgeschwindigkeit** — das ist die
  Haftbedingung.
- **Periodisch:** Gegenüberliegende Ränder werden identifiziert. Numerisch bedeutet
  das zirkulante Matrizen (Modulo-Indizierung der Differenzensterne) und das
  Weglassen des doppelten Endpunkts.

### 8.2 Das Wandproblem und die Thom-Formel

An einer Haftwand müssen **beide** Bedingungen gleichzeitig gelten:

$$
\Psi_w = \text{konst.} \quad (\text{keine Durchströmung}), \qquad
\partial_n\Psi_w = 0 \quad (\text{Haftung})
$$

Die Poisson-Gleichung kann pro Rand aber nur *eine* Bedingung aufnehmen. Wohin also
mit der zweiten? Zugleich fehlt für die Transportgleichung ein Randwert für $\omega$
— es gibt keine physikalische Vorschrift, die die Wandwirbelstärke direkt festlegt.

Eine Größe ist überbestimmt, eine andere unterbestimmt. Das ist der Hinweis: **die
überzählige Bedingung muss die fehlende ersetzen.** Die Taylor-Entwicklung der
Stromfunktion von der Wand aus ins Fluid leistet die Übersetzung:

$$
\Psi_1 = \Psi_w + h\underbrace{\partial_n\Psi_w}_{=\,0} + \frac{h^2}{2}\partial_n^2\Psi_w + \mathcal{O}(h^3)
$$

Mit $\omega_w = -\nabla^2\Psi_w \approx -\partial_n^2\Psi_w$ (die tangentiale
Krümmung verschwindet, da $\Psi$ entlang der Wand konstant ist) folgt die
**Thom-Formel**:

$$
\boxed{\;\omega_w \approx \frac{2(\Psi_w - \Psi_1)}{h^2}
\;\overset{\Psi_w = 0}{=}\; -\frac{2\Psi_1}{h^2}\;}
$$

**Das ist der physikalisch tiefsinnigste Punkt der ganzen Herleitung.** In Abschnitt
4.1 wurde festgestellt, dass die 2D-Wirbeltransportgleichung keinen Quellterm hat —
Wirbelstärke entsteht im Fluidinneren nicht. Die Thom-Formel zeigt, wo sie
stattdessen entsteht: **an der Wand.** Die Haftbedingung erzwingt einen
Geschwindigkeitsgradienten, und dieser Gradient *ist* Wirbelstärke. Die Wand ist die
Quelle der Grenzschicht, der Ablösung und letztlich der gesamten Wirbelstraße.

Die Formel ist von erster Ordnung genau in $h$ — hier zahlt sich die logarithmische
Gitterverfeinerung ein zweites Mal aus, weil $h$ genau dort klein ist, wo diese
Formel angewendet wird.

### 8.3 Der Außenrand

Weit vom Zylinder entfernt ist die Strömung nahezu wirbelfrei. Angesetzt wird die
**Potentialströmung um den Zylinder**:

$$
\Psi_\infty = \sin\theta\left(r - \frac{1}{r}\right), \qquad \omega = 0
$$

Das ist die exakte Lösung der reibungsfreien Umströmung (Superposition aus
Parallelströmung und Dipol) und erfüllt $\nabla^2\Psi = 0$ — konsistent mit
$\omega = 0$.

Die Einschränkung: Diese Bedingung ist am Nachlaufrand physikalisch nicht korrekt,
denn dort verlässt Wirbelstärke das Gebiet. Bei $r_\text{max} = 30$ ist die
Wirbelstärke dort so stark abgeklungen, dass der Fehler tolerierbar bleibt. Sauberer
wäre eine konvektive Ausströmbedingung
$\partial_t\omega + u_\infty\partial_x\omega = 0$.

### 8.4 Zusammenfassung der Ränder

| Rand | Stromfunktion | Wirbelstärke | Physik |
|---|---|---|---|
| Wand ($r_{\ln}=0$) | $\Psi = 0$ (Dirichlet) | Thom-Formel | Haftung; Wirbelquelle |
| Außen ($r=30$) | $\Psi = \sin\theta(r-1/r)$ | $\omega = 0$ | ungestörte Fernströmung |
| $\theta$ | periodisch | periodisch | geschlossener Umfang |

---

## 9. Diskretisierung: von der Analysis zur linearen Algebra

### 9.1 Taylor-Reihe als Fundament

Jedes Finite-Differenzen-Schema folgt aus der Taylor-Entwicklung. Für die zentrale
zweite Ableitung:

$$
\partial_x^2 \Phi_j \approx \frac{\Phi_{j-1} - 2\Phi_j + \Phi_{j+1}}{h^2} + \mathcal{O}(h^2)
$$

Die Fehlerordnung ist kein Nebenaspekt, sondern die Größe, an der später die
Korrektheit gemessen wird (Abschnitt 14).

### 9.2 Operatoren als Matrizen

Wendet man ein Schema auf alle Gitterpunkte an und ordnet die Feldwerte in einem
Vektor an, wird jeder Differentialoperator zu einer dünnbesetzten Matrix:

$$
\partial_x \;\longrightarrow\; \boldsymbol{D}_x, \qquad
\nabla^2 \;\longrightarrow\; \boldsymbol{L} = \boldsymbol{D}_{xx} + \boldsymbol{D}_{yy}
$$

Die 2D-Matrizen entstehen per **Kronecker-Produkt** aus den 1D-Matrizen — was die
Trennbarkeit der Ableitungsrichtungen ausnutzt und die Konstruktion sowohl kompakt
als auch fehlerarm macht. Ab hier ist die gesamte Analysis lineare Algebra.

In `src/operators.py`: Zentraldifferenzen der Fehlerordnung 2 und 4, einseitige
Randschemata, periodische (zirkulante) Varianten.

### 9.3 Randknoten über boolesche Masken

Ränder werden als boolesche Vektoren markiert (`wand`, `aussen`). Die
Poisson-Systemmatrix ersetzt in den Randzeilen die Laplace-Zeile durch die Identität:

$$
\boldsymbol{A} = \mathrm{diag}(\neg\boldsymbol{B})\,\boldsymbol{L} + \mathrm{diag}(\boldsymbol{B})
$$

Damit steht in Randzeilen schlicht $\Psi_i = \Psi_{i,\text{Rand}}$, in Innenzeilen die
Differenzengleichung. Ein einheitliches lineares System für beides.

### 9.4 Semi-Diskretisierung (Linienmethode)

Wird nur der Ort diskretisiert, bleibt ein System **gewöhnlicher**
Differentialgleichungen in der Zeit:

$$
\dot{\boldsymbol{\omega}} = f(\boldsymbol{\omega})
$$

Der Vorteil ist methodisch: Orts- und Zeitdiskretisierung lassen sich getrennt
entwerfen, getrennt analysieren und getrennt austauschen.

---

## 10. Der Lösungsalgorithmus

Die einzige Nichtlinearität des Systems steckt im Konvektionsterm
$\boldsymbol{u}\cdot\boldsymbol{\nabla}\omega_z$ — Geschwindigkeit mal
Wirbelstärkegradient, und beide hängen voneinander ab. Der entscheidende Kunstgriff:
Innerhalb einer Auswertung das Geschwindigkeitsfeld als **eingefroren** behandeln.
Dann zerfällt jeder Zeitschritt in drei *lineare* Teilprobleme, in genau der
Reihenfolge der physikalischen Kausalkette aus Abschnitt 5:

```
        ┌─────────────────────────────────────────────────┐
        │  ω(t)  — der Zustand                            │
        └────────────────────┬────────────────────────────┘
                             ↓
        (1) POISSON:   ∇²Ψ = −ω,  mit Randwerten
                       → lineares Gleichungssystem
                             ↓
        (2) CAUCHY-RIEMANN:  u = ∂_yΨ,  v = −∂_xΨ
                       → reine Matrix-Vektor-Produkte
                             ↓
        (3) THOM:      ω_w aus Ψ₁ an der Wand setzen
                       → Wirbelerzeugung am Rand
                             ↓
        (4) TRANSPORT: ω̇ = ν∇²ω − (u·∇)ω  mit bekanntem u
                       → linear in ω
                             ↓
        (5) ZEITSCHRITT (RK4) → ω(t+Δt)
                             ↓
                    ┌────────┘
                    ↓  zurück zu (1)
```

**Der entscheidende Performance-Aspekt:** Das teuerste Element ist Schritt (1) — das
Lösen des Poisson-Systems. Aber dessen Matrix ist *zeitlich konstant*; nur die rechte
Seite ändert sich. Deshalb wird sie **einmal** vor der Zeitschleife LR-zerlegt
(`scipy.sparse.linalg.splu`), und pro Schritt verbleiben nur Vorwärts- und
Rückwärtseinsetzen. Das ist der Unterschied zwischen Sekunden und Stunden Rechenzeit.

Der durch das Einfrieren entstehende Fehler ist von derselben Ordnung wie der
Zeitschrittfehler selbst — also ohne Verlust an Gesamtgenauigkeit.

---

## 11. Zeitintegration und Stabilität

Nach der Semi-Diskretisierung steht ein ODE-System, und es braucht ein
Zeitschrittverfahren. Zwei Analysen bestimmen die Wahl.

### 11.1 Übertragungsverhalten (modifizierte Wellenzahl)

- **Zentraldifferenzen sind dispersiv, aber nicht dissipativ.** Sie verfälschen die
  Phasengeschwindigkeit kurzwelliger Moden, dämpfen aber keine Amplitude — sie
  erzeugen keine numerische Diffusion.
- **Einseitige (Aufwind-)Differenzen 1. Ordnung sind dissipativ.** Sie wirken wie eine
  künstliche Viskosität. Das stabilisiert, verfälscht aber die effektive
  Reynolds-Zahl — die Simulation rechnet dann implizit ein zäheres Fluid als
  vorgegeben.

### 11.2 Von-Neumann-Stabilitätsanalyse

Man setzt eine räumliche Fourier-Mode

$$
\Phi(x_j, t_l) = \xi^l\,\mathrm{e}^{\mathrm{i}\lambda j h_x}
$$

in das volldiskretisierte Schema ein und fordert, dass der Verstärkungsfaktor
betragsmäßig nicht wächst:

$$
\left|\frac{\Phi(x_j, t_{l+1})}{\Phi(x_j, t_l)}\right| = |\xi| \overset{!}{\leq} 1
$$

Daraus folgen die beiden zentralen Grenzen:

$$
\mathrm{CFL}_x = |u|\frac{h_t}{h_x} \lesssim 1
\qquad\text{(Konvektion)}
$$

$$
\mathrm{CFL}_{xx} = \nu\frac{h_t}{h_x^2} \lesssim \tfrac{1}{4}
\qquad\text{(Diffusion)}
$$

Die CFL-Bedingung besagt physikalisch: Information darf pro Zeitschritt höchstens
etwa eine Gitterzelle weit laufen. Das **Quadrat** in der Diffusionsbedingung ist
folgenreich — Gitterverfeinerung verteuert explizite Verfahren doppelt.

Auf dem logarithmischen Polargitter ist die kleinste Zelle an der Wand; sie bestimmt
die Zeitschrittweite für das gesamte Gebiet. Das ist der Preis für die
Grenzschichtauflösung.

### 11.3 Warum RK4

Eine klassische Falle: **Das explizite Euler-Verfahren mit Zentraldifferenzen ist für
reine Konvektion immer instabil.** Sein Stabilitätsgebiet enthält kein Stück der
imaginären Achse — und genau dort liegen die Eigenwerte der
Zentraldifferenzen-Konvektion. Klassisches RK4 dagegen schneidet die imaginäre Achse
und ist damit das passende Verfahren. (Mit Viskosität rutschen die Eigenwerte in die
linke Halbebene, dann geht Euler mit kleinem $\Delta t$ gerade noch — aber „geht
gerade noch" ist keine tragfähige Verfahrensgrundlage.)

> **Merksatz:** Räumliches Schema und Zeitschrittverfahren sind ein *Paar*.
> Stabilität ist keine Eigenschaft des einen oder des anderen, sondern der
> Kombination.

Die Implementierung nutzt RK4 mit $\Delta t = 0{,}02$ bei $\mathrm{Re} = 100$.

### 11.4 Hybridschema für höhere Reynolds-Zahlen

Bei größerem $\mathrm{Re}$ oder gröberem Gitter wird die **Gitter-Reynolds-Zahl**

$$
\mathrm{Re}_\text{Zelle} = |u|\,h\,\mathrm{Re}
$$

zum Problem: Überschreitet sie 2, destabilisieren Zentraldifferenzen die Lösung
(oszillatorische Moden). Der Löser `src/immersed.py` mischt deshalb knotenweise
zwischen Zentral- und Aufwindverfahren:

$$
\theta = \min\left(1, \frac{2}{\mathrm{Re}_\text{Zelle}}\right),
\qquad
\text{Konvektion} = \theta \cdot \text{zentral} + (1-\theta)\cdot\text{aufwind}
$$

Wo das Gitter fein genug auflöst, wird dispersionsarm zentral gerechnet, sonst stabil
aufwind. Der messbare Effekt: Die Strouhal-Zahl bei $\mathrm{Re} = 100$ steigt von
0,142 (reines Aufwind) auf 0,163 — die künstliche Viskosität des Aufwindverfahrens
verfälschte die Ablösefrequenz spürbar.

---

## 12. Symmetriebrechung: warum die Wirbelstraße nicht von allein kommt

Ein physikalisch lehrreicher Befund: Startet man die Simulation impulsiv aus der Ruhe
($\omega = 0$), stellt sich **kein** Wirbelablösen ein. Stattdessen wächst ein
symmetrisches, stationäres Wirbelpaar hinter dem Zylinder — und bleibt.

Der Grund: Gleichungen, Gitter und Randbedingungen sind exakt symmetrisch zur
$x$-Achse. Die symmetrische Lösung existiert bei $\mathrm{Re} = 100$ mathematisch
weiterhin; sie ist nur **instabil**. Im Experiment sorgen Rauschen, Rauheit und
Turbulenz der Anströmung dafür, dass die Instabilität sofort angeregt wird. Die
Simulation hat dieses Rauschen nicht — bis auf Rundungsfehler, die zu klein und zu
symmetrisch verteilt sind.

Die Symmetrie muss also gezielt gestört werden. Die Implementierung dreht dazu kurz
den Zylinder ($u_\text{Wand} = 0{,}5$ für $1 < t < 3$); alternativ funktioniert ein
kurzes Kippen der Anströmung. Danach schaukelt sich die Instabilität selbständig zur
periodischen Ablösung auf — ab $t \approx 15$ steht die Kármán'sche Wirbelstraße.

Wichtig zum Verständnis: Die Störung *erzeugt* die Wirbelstraße nicht, sie
**triggert** sie nur. Form und Frequenz des sich einstellenden Grenzzyklus sind von
der Art der Störung unabhängig — das ist gerade das Kennzeichen einer
Hopf-Verzweigung, die beim Kreiszylinder bei $\mathrm{Re}_\text{krit} \approx 47$ liegt.

---

## 13. Physikalische Auswertung

### 13.1 Strouhal-Zahl

Die dimensionslose Ablösefrequenz ist die zentrale Vergleichsgröße:

$$
\mathrm{St} = \frac{f \cdot 2r_0}{u_\infty} = \frac{f\,D}{u_\infty}
$$

Gemessen wird sie an einer Sonde im Nachlauf ($r = 3$, $\theta = 0$): Die
Quergeschwindigkeit $u_\theta$ oszilliert dort periodisch; aus den mittleren
Abständen der Nulldurchgänge in aufsteigender Richtung folgt die Periode und daraus
$f$. Die erste Hälfte des Signals wird verworfen, um Einschwingvorgänge auszuschließen.

**Ergebnis:** $\mathrm{St} = 0{,}168$ gegenüber $\approx 0{,}16$–$0{,}17$ in der
Literatur für $\mathrm{Re} = 100$.

### 13.2 Kraftbeiwerte

Der kartesische Löser berechnet zusätzlich Widerstand und Auftrieb. Weil die
Wirbel-Stromfunktions-Formulierung den Druck eliminiert hat, muss er für die
Kraftberechnung rekonstruiert werden — aus der Divergenz der Impulsgleichung:

$$
\nabla^2 p = 2\left(\partial_x u\,\partial_y v - \partial_y u\,\partial_x v\right)
$$

Die Wandschubspannung folgt direkt aus der Wandwirbelstärke:
$\tau = \omega_w/\mathrm{Re}$. Integration von Druck und Schubspannung über die
Körperoberfläche liefert $C_W$ und $C_A$.

Ein instruktives Detail: Der Körper muss druckseitig entkoppelt werden
($\partial p/\partial n = 0$ an der Wand). Sonst „leckt" der Staudruck durch das
Körperinnere zur Rückseite und halbiert den berechneten Formwiderstand.

### 13.3 FTLE — Lagrange-kohärente Strukturen

Wirbelstärkefelder zeigen, *wo* Rotation ist. Die Transportstrukturen der Strömung —
welche Fluidpakete sich trennen oder zusammenfinden — werden dagegen über den
zeitlich abgeschätzten Ljapunow-Exponenten sichtbar.

Ein Partikelgitter wird durch das zeitabhängige Geschwindigkeitsfeld integriert:

$$
\boldsymbol{p}(t) = \boldsymbol{p}(t_0) + \int_{t_0}^{t} \dot{\boldsymbol{p}}(\tau)\,d\tau
$$

Die Jacobimatrix $J_{\boldsymbol{p}}$ dieser Flussabbildung wird über finite
Differenzen benachbarter Partikel approximiert, und der maximale Eigenwert des
rechten Cauchy-Green-Deformationstensors $C = J^\top J$ liefert das
Streckungsmaß:

$$
\sigma(t) = \frac{1}{t - t_0}\ln\left(\sqrt{\lambda_\text{max}\{J^\top J\}}\right)
$$

Da $C$ eine $2\times 2$-Matrix ist, wird $\lambda_\text{max}$ analytisch als
Nullstelle des quadratischen charakteristischen Polynoms berechnet — deutlich
schneller als ein numerischer Eigenwertlöser bei Zehntausenden von Partikeln.

**Vorwärts** integriert markieren die Maxima instabile Mannigfaltigkeiten (rot,
abstoßend), **rückwärts** die stabilen (blau, anziehend). Ihre Schnittpunkte sind die
Transportbarrieren, die die Wirbel der Straße gegeneinander abgrenzen.

---

## 14. Validierung

Für die Zylinderumströmung existiert keine analytische Lösung. Deshalb wird die
Numerik **vorher** an einem Problem bewiesen, dessen Lösung exakt bekannt ist.

### 14.1 Der Taylor-Green-Wirbel

Er gehört zu den wenigen Fällen mit exakter Lösung der *nichtlinearen* Gleichungen.
Der Grund ist eine gezielte Annahme: Ist die Wirbelstärke überall ein festes
Vielfaches der Stromfunktion ($\omega = \eta\Psi$, Beltrami-Strömung), dann gilt

$$
\boldsymbol{u}\cdot\boldsymbol{\nabla}\Psi
= \partial_y\Psi\,\partial_x\Psi - \partial_x\Psi\,\partial_y\Psi = 0,
$$

der Gradient der Stromfunktion steht also senkrecht auf der Geschwindigkeit — **die
einzige Nichtlinearität verschwindet identisch.** Der Ansatz $\Psi \sim \sin x \sin y$
liefert über die Poisson-Gleichung $\eta = 2$, der Ansatz
$\Psi \sim \mathrm{e}^{-\nu\eta t}$ den Zeitabfall:

$$
\Psi = \sin x \sin y\,\mathrm{e}^{-2\nu t}, \qquad \omega = 2\Psi
$$

Mit periodischen Randbedingungen ist das zugrundeliegende Gebiet topologisch ein
**flacher Torus** $T^2 = S^1 \times S^1$.

### 14.2 Drei Prüfebenen

**1. Konvergenzordnung messen, nicht behaupten.** Der Fehler gegen die analytische
Lösung bei 16/32/64 Punkten muss die nominelle Ordnung ergeben. Gemessen: 1,98 und
2,00 — Punktlandung auf Ordnung 2. Eine falsche Randbehandlung würde genau hier
auffliegen.

**2. Erhaltungsgrößen als Diffusionsdetektor.** Die Enstrophie
$\mathcal{E} = \tfrac{1}{2}\int\omega^2\,dA$ muss im reibungsfreien Fall konstant
bleiben (Abschnitt 5). Analytisch: $\mathcal{E}(t) = 2\pi^2\mathrm{e}^{-4\nu t}$. Jede
Abweichung bei $\nu = 0$ ist *numerische* Diffusion. Gemessen (RK4 +
Zentraldifferenzen): Drift exakt 0.

**3. Grenzfälle.** $\nu = 0$ liefert eine stationäre Lösung (relativer L2-Fehler
$3\cdot 10^{-14}$, Maschinengenauigkeit); $\nu > 0$ den exponentiellen Abfall mit der
bekannten Rate $4\nu$.

| Messung | Ergebnis |
|---|---|
| Rel. L2-Fehler $\omega$ ($\nu = 0{,}1$, 64×64, RK4) | $6{,}4\cdot10^{-4}$ |
| Rel. L2-Fehler ($\nu = 0$, stationär) | $3\cdot10^{-14}$ |
| Enstrophie-Drift bei $\nu = 0$ | 0 |
| Konvergenzordnung (16→32→64) | 1,98 / 2,00 |

### 14.3 Validierung des Zylinders

Für den Zylinder tritt an die Stelle der analytischen Lösung die
**Literaturvalidierung** über die Strouhal-Zahl. $\mathrm{St} = 0{,}168$ gegenüber
$\approx 0{,}16$–$0{,}17$ ist der Beleg, dass die gesamte Kette — Formulierung,
Koordinatentransformation, Randbedingungen, Zeitintegration — konsistent ist.

Zusätzlich prüfbar: Auf dem Log-Polargitter muss die Potentialströmung ($\omega = 0$)
stationär gehalten werden, und bei $\mathrm{Re} \approx 20$–$40$ muss sich das
bekannte stationäre Wirbelpaar einstellen (unterhalb der kritischen Reynolds-Zahl).

> **Die methodisch wichtigste Regel:** Niemals Numerik und Geometrie gleichzeitig
> debuggen. Der Taylor-Green-Wirbel validiert die vollständige Numerik auf dem
> einfachsten Gebiet. Erst wenn das bewiesen ist, kommt die schwierige Geometrie
> hinzu — dann kann jeder neue Fehler nur noch von der Geometrie stammen.

---

## 15. Die durchgehende Kausalkette im Überblick

```
  Massenerhaltung                    Impulserhaltung (Newton II)
        │                                      │
        ↓                                      ↓
  ∇·u = 0                          Navier-Stokes-Gleichung
  (Zwangsbedingung)                (mit Druck ohne eigene Gleichung)
        │                                      │
        │  Handgriff 2:                        │  Handgriff 1:
        │  Variablenwahl                       │  Rotation anwenden
        ↓                                      ↓
  Stromfunktion Ψ                    Wirbeltransportgleichung
  u = ∂_yΨ, v = −∂_xΨ                ∂_tω + u·∇ω = ν∇²ω
  (Kontinuität automatisch)          (2D: kein Streckungsterm,
        │                             kein Druck, kein Quellterm)
        └──────────────┬───────────────────────┘
                       ↓
              ω = −∇²Ψ  (Poisson)
              → System geschlossen: 2 Felder, 3 Gleichungen
                       ↓
        Entdimensionalisierung (Buckingham)
              → ein Parameter: Re = D·u∞/ν
                       ↓
        Koordinaten: r_ln = ln r
              → Rand ist Koordinatenlinie
              → Grenzschicht automatisch aufgelöst
              → ∇² = e^(−2r_ln)(∂²_r_ln + ∂²_θ)
                       ↓
        Randbedingungen
              → Wand: Ψ = 0 + Thom-Formel ω_w = −2Ψ₁/h²
                ↑ HIER entsteht alle Wirbelstärke
              → Außen: Potentialströmung, ω = 0
              → θ: periodisch
                       ↓
        Diskretisierung (Taylor → Matrizen → Kronecker)
              → ω̇ = f(ω), ODE-System
                       ↓
        Dreischritt-Schleife + RK4 (CFL- und Diffusionsgrenze)
              + LR-Zerlegung der Poisson-Matrix gecacht
                       ↓
        Symmetriebrechung (Störung triggert Instabilität)
                       ↓
        ═══════════════════════════════════════════
         Kármán'sche Wirbelstraße,  St = 0,168
        ═══════════════════════════════════════════
                       ↓
        Auswertung: Strouhal-Zahl, C_W/C_A, FTLE
```

### Die sieben übertragbaren Prinzipien

1. **Vereinfache zuerst das Problem, nicht den Code.** (2D, inkompressibel, laminar)
2. **Eliminiere Störendes durch Operatoren, erfülle Zwänge durch Konstruktion.**
   (Rotation gegen den Druck, Stromfunktion gegen die Kontinuität)
3. **Mache Ableitungen zu Matrizen — dann ist Analysis nur noch lineare Algebra.**
4. **Zerlege Nichtlineares in eine Folge linearer Teilprobleme und cache, was sich
   nicht ändert.** (Dreischritt-Schleife, LR-Zerlegung)
5. **Die Physik eines Problems steckt fast vollständig in den Randbedingungen.**
   (Thom-Formel: die Wand als Wirbelquelle)
6. **Validiere am analytischen Grenzfall, bevor Geometrie hinzukommt — und miss
   Konvergenzordnungen, statt sie zu glauben.**
7. **Passe die Koordinaten der Geometrie an, nicht den Löser.**

