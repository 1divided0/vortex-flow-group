# Thema 6: Anwendungsbeispiele

Nach ein paar Hinweisen zur Programmierung, wird hier der Taylor-Green-Wirbel als Beispiel für die numerische Untersuchung vorgestellt.

[TOC]


<!------------------------------------------------------------------------------
Bewährte Programmierpraktiken
------------------------------------------------------------------------------->
## Bewährte Programmierpraktiken

In dem [Begleitmaterial](Begleitmaterial) befinden sich bereits Beispiele für bewährte Programmierpraktiken. Denn das Einhalten von Konvention macht die Entwicklung von Programmcode, für alle Beteiligten, sehr viel leichter. Genauso wie das wissenschaftliche Arbeiten zum Verfassen von fundierten Berichten gehört, so ist auch der Programmierstil für die Entwicklung von Programmcode maßgeblich. Dazu folgende Anmerkungen.

1. **Dokumentation:** Der Programmcode lässt sich durch Kommentare in Sinnabschnitte unterteilen. Außerdem bietet es sich an, mit einem sog. Docstring, kurze Funktionsbeschreibungen vorzunehmen.

2. **Vermeidung von Redundanz:** Werte die mehrmals vorkommen, sollten als Variablen verwendet werden. Routinen welche mehrmals vorkommen, sollten als Funktionen ausgelagert werden.

3. **Verwendung von Parametern:** Bei einer Simulation gibt es immer Parameter. Es bietet sich dafür an, die Hauptfunktion mit sog. Keyword-Arguments zu versehen, um dadurch Voreinstellungen für die möglichen Parameter zu treffen.

4. **Typisierung:** Durch das sog. Type-Hinting kann bei Programmiersprachen mit impliziter Typisierung explizit auf den Datentyp hingewiesen werden. Dies ist gerade bei Funktionsargumenten und Rückgabewerten sehr nützlich, da so Missverständnisse vermieden werden.

<!--
5. **Handhabung von großen Datenmengen:** Bei Skriptsprachen wie Matlab oder Python erfolgt der Funktionsaufruf meistens durch sog. Wertparameter (call by value) und nicht durch sog. Referenzparameter (call by reference), wodurch das gesamte Funktionsargument bei dem Funktionsaufruf kopiert wird. Da Speicheroperationen aber in der Regel sehr langsam sind, ist es in diesem Fall günstiger, große Objekte – wie z. B. die Differenzenmatrizen – global zu definieren, um sie nicht ständig neu zu initialisieren.
-->


<!------------------------------------------------------------------------------
Taylor-Green-Wirbel
------------------------------------------------------------------------------->
## Taylor-Green-Wirbel

Für den hier betrachteten zweidimensionalen Wirbel hat erstmals Geoffrey I. Taylor allein eine analytische Lösung hergeleitet (On the decay of vortices in a viscous fluid – 1923). Fälschlicherweise wird dieses Resultat oft auf seine Zusammenarbeit mit Albert E. Green von 1937 zurückgeführt, obwohl diese Veröffentlichung gar nicht die besagte Herleitung enthält. Der dort betrachtete Wirbel ist nämlich dreidimensional und baut insofern auf den vorherigen Erkenntnissen auf. Trotzdem werden in der Literatur beide Fälle meist als »Taylor-Green-Wirbel« bezeichnet.

### Analytische Lösung

Die Herleitung ist so simpel wie genial. Ausgehend von der Wirbeltransportgleichung für zweidimensionale inkompressible Strömungen

$$
\partial_t \omega_z + \boldsymbol{u} \cdot \nabla\omega_z = \nu\nabla^2\omega_z
$$

wird die Annahme einer allgemeinen Beltrami-Strömung getroffen, dass sich die Wirbelstärke entlang von Stromlinien konstant verhält, also in diesem Fall ein $\eta$-Faches der Stromfunktion darstellt:

$$
\omega_z = \eta\Psi
$$

Die Poisson-Gleichung hängt dann nur noch von der Stromfunktion ab und in der Wirbeltransportgleichung entfällt der Konvektionsterm, da der Gradient der Stromfunktion senkrecht zur Geschwindigkeit ist. Die beiden Ansätze

$$
\begin{alignat*}{2}
\Psi &\sim \mathrm{e}^{-\nu\eta t}, &&\quad\text{für die Wirbeltransportgleichung, und}\\
\Psi &\sim \sin(x)\sin(y), &&\quad\text{für die Poisson-Gleichung,}
\end{alignat*}
$$

lassen sich kombiniert einsetzten, was zu einem bestimmten Wert für $\eta$ und somit zur analytischen Lösung führt.

---
> **Aufgabe (Herleitung der analytischen Lösung des Taylor-Green-Wirbels)**
>
> Leitet die analytische Lösung des Taylor-Green-Wirbels her.

### Visualisierungsbeispiel

Die Lösung lässt sich durch das Verfolgen sog. Lagrange-Partikel visualisieren. Hierzu wird auf dem Rechengebiet eine Startverteilung masseloser Punkte gewählt und deren Position nach der Lagrange'schen Betrachtungsweise in der Zeit durch Interpolation des Geschwindigkeitsvektorfeldes entwickelt.

$$
\boldsymbol{p}(t+h_t) \approx \boldsymbol{p}(t) + h_t \cdot \dot{\boldsymbol{p}}(t) = \boldsymbol{p}(t) + h_t \cdot \boldsymbol{u}(\boldsymbol{p}(t),t)
$$

> **Begleitmaterial (Visualisierung des Taylor-Green-Wirbels)**
>
> [![Matlab Beispiel](.Dateien/Plaketten/Matlab_Beispiel.svg)](Begleitmaterial/Taylor_Green_Wirbel.m) [![Python Beispiel](.Dateien/Plaketten/Python_Beispiel.svg)](Begleitmaterial/Taylor_Green_Wirbel.py)
>
> ![Taylor-Green-Wirbel ohne Viskosität](.Dateien/Bilder/Taylor_Green_Wirbel.gif)
>
> ![Taylor-Green-Wirbel mit Viskosität](.Dateien/Bilder/Taylor_Green_Wirbel_nu.gif)

---
> **Aufgabe (Topologie des Taylor-Green-Wirbels)**
>
> Welche Topologie liegt dem Taylor-Green-Wirbel zugrunde? Welche Randbedingungen lassen sich hier setzen?

---
> **Aufgabe (Vergleich der Enstrophie)**
>
> Vergleicht bei eurer Simulation die Enstrophie mit ihrem analytischen Ergebnis, indem ihr das zweidimensionale Volumenintegral über die Riemann-Summe approximiert. Wie hängt der Verlauf von dem gewählten Differenzenverfahren ab? Lässt sich numerische Diffusion erkennen?

---
> **Aufgabe (Vergleich der Laufzeit)**
>
> Wie verändert sich die Laufzeit eures Lösungsalgorithmus mit der Ordnung der gewählten Differenzenschemata?
