# Teil 1: Erhaltungsgleichungen

...

[TOC]


<!------------------------------------------------------------------------------
Voraussetzungen
------------------------------------------------------------------------------->
## Voraussetzungen

Die Formulierung der Gleichungen stützt sich auf der Kontinuumshypothese. Demnach befinden sich in einem hinreichend kleinen Fluidvolumen ausreichend Moleküle, sodass von einem Kontinuum ausgegangen werden kann. Etwas genauer lässt sich die Bedingung mit der Knudsen-Zahl $\mathrm{Kn}$ definieren, welche das Verhältnis der mittleren freien molekularen Weglänge $l$ zur charakteristischen Länge des Strömungsfeldes $L$ (z. B. dem Durchmesser eines durchströmten Rohres) beschreibt.

$$
    \mathrm{Kn} = \frac{l}{L}
$$

> **Tabelle (Strömungsart nach Knudsen-Zahl)**
>
> Für eine Kontinuumsströmung darf die Knudsen-Zahl nicht größer als ein Hundertstel sein.
>
> | Bereich                         | Art der Strömung        |
> |:-------------------------------:|:------------------------|
> | $10 < \mathrm{Kn}$              | freie Molekularströmung |
> | $1/10 < \mathrm{Kn} \leq 10$    | Übergangsströmung       |
> | $1/100 < \mathrm{Kn} \leq 1/10$ | Gleitströmung           |
> | $\mathrm{Kn} \leq 1/100$        | Kontinuumsströmung      |

Außerdem lässt sich eine Fluid Strömung auf zwei unterschiedliche Weisen betrachten.

> **Abbildung (Euler'sche und Lagrange'sche Betrachtungsweise)**
>
> Bei der Euler'schen Betrachtungsweise ist das Koordinatensystem ortsfest. Bei der Lagrange'schen Betrachtungsweise hingegen, wird das Koordinatensystem mit der Strömung eines einzelnen Fluidpartikels mitgeführt, wobei dessen Hülle zwar beliebig flexibel, aber dennoch undurchlässig ist. Nach dem klassischen Relativitätsprinzip kann man sich leicht davon überzeugen, dass sich beide Betrachtungsweisen sowohl mathematisch, als auch physikalisch ineinander überführen lassen. In der Abbildung ist das Bezugssystem in schwarz und das Fluidpartikel in blau angedeutet.
>
> ![Betrachtungsweise nach Euler und nach Lagrange](../Dateien/Bilder/Euler_vs_Lagrange.svg)

Da die Euler'sche Betrachtungsweise der eines außenstehenden Betrachters entspricht, ist sie womöglich etwas intuitiver als die Lagrange'sche. Darum wird im Folgenden auf diese Weise verfahren, um die strömungsmechanischen Erhaltungsgleichungen herzuleiten.


<!------------------------------------------------------------------------------
Massenerhaltung (alias Kontinuitätsgleichung)
------------------------------------------------------------------------------->
## Massenerhaltung (alias Kontinuitätsgleichung)

Für die Herleitung der Massenerhaltungsgleichung wird nun also von einem infinitesimalen, raum- und ortsfesten Kontrollvolumen $\Omega$ ausgegangen. Dafür wird zunächst der Massenstrom über die einzelnen Raumachsen $i$ bilanziert.

> **Abbildung (Massenstrombilanz)
>
> $$ \partial_t\, m = \partial_t (\rho \Omega) = \sum_i (\dot m_{\mathrm{ein}} - \dot m_{\mathrm{aus}})_i = \dot m $$
>
> ![Massenstrombilanz](../Dateien/Bilder/Massenstrombilanz.svg)

Der an der vorderen Querschnittsfläche einer Achse $i$ eintretende Massenstrom ergibt sich durch Multiplikation dieser Querschnittsfläche $\Omega/dx_i$ mit der senkrecht zu ihr stehenden Eintrittsgeschwindigkeit $\dot x_i$ und der entsprechenden Fluiddichte $\rho$.

$$
    (\dot m_{\mathrm{ein}})_i = \rho \dot x_i \Omega/dx_i
$$

Der an der hinteren Querschnittsfläche einer Achse $i$ austretende Massenstrom ergibt sich wiederum durch die Taylorreihe des eintretenden Massenstroms, entwickelt an der Stelle $x_i$ und ausgewertet an der Stelle $x_i + dx_i$.

$$
    (\dot m_{\mathrm{aus}})_i = (\rho \dot x_i /dx_i + \partial_{x_i}(\rho \dot x_i) + \mathcal{O}(dx_i)) \Omega
$$

Wird das in die Bilanzgleichung für den Massenstrom eingesetzt, so ergibt sich

$$
    \partial_t(\rho V) = -\sum_i (\partial_{x_i}(\rho \dot x_i) + \mathcal{O}(dx_i)) \Omega.
$$

Da das Kontrollvolumen $\Omega$ nach der Euler'schen Betrachtungsweise nicht von der Zeit $t$ abhängt, kann es aus der Gleichung gekürzt werden. Aufgrund seiner Infinitesimalität verschwinden außerdem die Terme höherer Ordnung $\mathcal{O}(dx_i)$.

$$
    \lim_{\Omega\to 0} \partial_t\, \rho = -\sum_i \lim_{dx_i\to 0} (\partial_{x_i}(\rho \dot x_i) + \underbrace{\mathcal{O}(dx_i)}_{\to 0})
$$

Übrig bleibt die Kontinuitätsgleichung in ihrer vollen Pracht

$$
\begin{equation} \tag{Kontinuitätsgleichung}
    \partial_t\, \rho + \nabla\cdot (\rho \mathbf{u}) = 0,
\end{equation}
$$

oder unter der Annahme von Inkompressibilität

$$
    \rho \nabla\cdot (\mathbf{u}) = 0.
$$
