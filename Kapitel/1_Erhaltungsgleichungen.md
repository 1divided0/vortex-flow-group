# Teil 1: Erhaltungsgleichungen

...

[TOC]


<!------------------------------------------------------------------------------
Voraussetzungen
------------------------------------------------------------------------------->
## Voraussetzungen

Die Formulierung der Gleichungen stützt sich auf der Kontinuumshypothese. Demnach befinden sich in einem hinreichend kleinen Fluidvolumen ausreichend viele Moleküle, sodass von einem Kontinuum ausgegangen werden kann. Etwas genauer lässt sich die Bedingung mit der Knudsen-Zahl $\mathrm{Kn}$ definieren, welche das Verhältnis der mittleren freien molekularen Weglänge $l$ zur charakteristischen Länge des Strömungsfeldes $L$ (z. B. dem Durchmesser eines durchströmten Rohres) beschreibt.

$$
    \mathrm{Kn} = \frac{l}{L}
$$

Für eine Kontinuumsströmung darf die Knudsen-Zahl nicht größer als ein Hundertstel sein.

> **Tabelle (Strömungsart nach Knudsen-Zahl)**
>
> ![Knudsen-Zahl](../Dateien/Tabellen/Knudsen_Zahl.svg)

Außerdem lässt sich eine Fluidströmung auf zwei unterschiedliche Weisen betrachten. Einerseits gibt es die Euler'sche Betrachtungsweise, bei der das Koordinatensystem ortsfest ist, und andererseits die Lagrange'sche Betrachtungsweise, bei der das Koordinatensystem mit der Strömung eines einzelnen Fluidpartikels mitgeführt wird, wobei dessen Hülle zwar beliebig flexibel, aber dennoch undurchlässig ist. Nach dem klassischen Relativitätsprinzip kann man sich leicht davon überzeugen, dass sich beide Betrachtungsweisen sowohl mathematisch, als auch physikalisch ineinander überführen lassen. In der folgenden Abbildung ist das Bezugssystem in schwarz und das Fluidpartikel in blau angedeutet.

> **Abbildung (Euler'sche und Lagrange'sche Betrachtungsweise)**
>
> ![Betrachtungsweise nach Euler und nach Lagrange](../Dateien/Bilder/Euler_vs_Lagrange.svg)

Da die Euler'sche Betrachtungsweise (links im Bild) der äußeren Perspektive entspricht, ist sie womöglich etwas intuitiver als die Lagrange'sche. Darum wird im Folgenden auf diese Weise formuliert.


<!------------------------------------------------------------------------------
Massenerhaltung (alias Kontinuitätsgleichung)
------------------------------------------------------------------------------->
## Massenerhaltung (alias Kontinuitätsgleichung)

Für die Herleitung der Massenerhaltungsgleichung wird von einem infinitesimalen ortsfesten Kontrollvolumen $\Omega$ ausgegangen. Dafür wird zunächst der Massenstrom über die einzelnen Raumachsen $i$ bilanziert.

$$
    \partial_t m = \partial_t (\rho \Omega) = \sum_i (\dot m_{\mathrm{ein}} - \dot m_{\mathrm{aus}})_i = \dot m
$$

> **Abbildung (Massenstrombilanz)
>
> ![Massenstrombilanz](../Dateien/Bilder/Massenstrombilanz.svg)

Der an der vorderen Querschnittsfläche einer Achse eintretende Massenstrom ergibt sich durch Multiplikation dieser Querschnittsfläche mit der senkrecht zu ihr stehenden Eintrittsgeschwindigkeit und der entsprechenden Fluiddichte.

$$
    (\dot m_{\mathrm{ein}})_i = \rho \dot x_i \Omega/dx_i
$$

Der an der hinteren Querschnittsfläche einer Achse austretende Massenstrom ergibt sich wiederum durch die Taylorreihe des eintretenden Massenstroms, entwickelt an der vorderen und ausgewertet an der hinteren Stelle.

$$
    (\dot m_{\mathrm{aus}})_i = (\rho \dot x_i /dx_i + \partial_{x_i}(\rho \dot x_i) + \mathcal{O}(dx_i)) \Omega
$$

Dieser Zusammenhang kann nun in die Bilanzgleichung für den Massenstrom eingesetzt werden.

$$
    \partial_t(\rho V) = -\sum_i (\partial_{x_i}(\rho \dot x_i) + \mathcal{O}(dx_i)) \Omega
$$

Da das Kontrollvolumen nach der Euler'schen Betrachtungsweise nicht von der Zeit abhängt, kann es aus der Gleichung gekürzt werden. Aufgrund seiner Infinitesimalität verschwinden außerdem die Terme höherer Ordnung.

$$
    \lim_{\Omega\to 0} \partial_t \rho = -\sum_i \partial_{x_i}(\rho \dot x_i)
$$

Übrig bleibt die Kontinuitätsgleichung in ihrer vollen Pracht

$$
    \partial_t\, \rho + \nabla\cdot (\rho \mathbf{u}) = 0,
$$

oder unter der Annahme von Inkompressibilität

$$
    \rho \nabla\cdot (\mathbf{u}) = 0.
$$
