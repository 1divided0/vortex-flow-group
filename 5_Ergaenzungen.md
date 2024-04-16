# Thema 5: Ergänzungen

Es wird zunächst eine analytische Lösung der eindimensionalen Transportgleichung hergeleitet, wie sie hier zum Einsatz kommt und damit als Ansatz dient, um die numerische Einflussnahme auf das Systemverhalten hinsichtlich der Stabilität zu untersuchen. Daraus ergeben sich diverse Möglichkeiten, die Approximationsverfahren zu kombinieren und auf die Problemstellung anzupassen, sodass die numerische Berechnung sicher und stabil vonstattengeht.

[TOC]


<!------------------------------------------------------------------------------
Analytische Lösung der eindimensionalen Transportgleichung
------------------------------------------------------------------------------->
## Analytische Lösung der eindimensionalen Transportgleichung

Für die eindimensionale Transportgleichung ohne Quellterm, lässt sich zu Vergleichszwecken eine analytische Lösung herleiten. Diese Gleichung ist gegeben durch

$$
\partial_t \Phi(x,t) = c\partial_x^2 \Phi(x,t) - u\partial_x \Phi(x,t),
$$

wobei _`c`_ die Dämpfungskonstante, _`u`_ die Transportgeschwindigkeit und _`Φ`_ die Strömungsgröße geschreibt. Mit der Wellenzahl _`λ:=2π/L`_, Frequenz _`f:=2π/T`_ und imaginären Einheit `i` lässt sich der folgende Ansatz wählen:

$$
\Phi(x,t) \sim \mathrm{e}^{\mathrm{i}(\lambda x - f t)}
$$

Wird dieser nun in die obige Transportgleichung eingesetzt, dann ergibt sich

$$
-\mathrm{i}f \mathrm{e}^{\mathrm{i}(\lambda x - f t)} = -c\lambda^2 \mathrm{e}^{\mathrm{i}(\lambda x - f t)} -iu\lambda \mathrm{e}^{\mathrm{i}(\lambda x - f t)}
$$

und gekürzt

$$
f = -\mathrm{i}c\lambda^2 + u\lambda,
$$

sodass die Frequenz in dem Ansatz substituiert werden kann, was durch Multiplikation mit einem Vorfaktor _`ξ`_ zur allgemeinen Lösung führt.

$$
\Phi(x,t) = \xi \mathrm{e}^{\mathrm{i}\lambda(x - u t)}\mathrm{e}^{-c\lambda^2 t}
$$


<!------------------------------------------------------------------------------
Numerisches Übertragungsverhalten
------------------------------------------------------------------------------->
## Numerisches Übertragungsverhalten

Das numerische Übertragungsverhalten gibt an, wie sich das gewählte Approximationsverfahren auf die Systemdynamik auswirkt. Wird die Transportgleichung zunächst ohne Diffusion (_`c=0`_) betrachtet, dann ergibt sich die Lösung durch die Exponentialfunktion mit komplexem Argument ohne Realteil

$$
\Phi(x,t) = \xi \mathrm{e}^{\mathrm{i}(\lambda x - f t)},
$$

sodass

$$
|\Phi(x,t)| = |\xi|.
$$

Wenn jedoch die Approximation der Ableitung die Wellenzahl im Raum bzw. die Frequenz in der Zeit derart modifiziert, dass im Argument der Exponentialfunktion ein Realteil entsteht, dann schlägt sich dies auf die Konvektion und Diffusion nieder, sodass

$$
|\Phi(x,t)| \approx |\xi\mathrm{e}^{\mathrm{i}(\lambda^\prime x - f^\prime t)}| \neq |\xi|.
$$

Da die Form der Lösung bekannt ist, kann in diesem Fall heuristisch davon ausgegangen werden, dass durch

$$
f = u\lambda \quad\Rightarrow\quad f^\prime \approx u\lambda^\prime
$$

ein Zusammenhang zwischen der räumlichen (konvektiven) und zeitlichen (diffusiven) Konzentrationsänderung besteht. So sollte die Wahl des Approximationsverfahrens mit Bedacht dieser Auswirkungen getroffen werden.

### Numerische Konvektion

Wenn sich die Approximation der Ableitung auf die räumliche Änderung der Konzentration auswirkt, dann ist von numerischer Konvektion die Rede. Um diesen Effekt zu untersuchen wird für die räumliche Ableitung

$$
t = 0
$$

gesetzt und für die zeitliche Ableitung mit der Heuristik

$$
\lambda = f / u \quad\Rightarrow\quad \lambda^\prime \approx f^\prime / u
$$

substituiert.

$$
\begin{alignat*}{3}
& \partial_x \Phi &&= \mathrm{i}\lambda\Phi &&\approx \mathrm{i}\lambda^\prime\Phi \\
\Rightarrow\quad& &&= \mathrm{i}(f / u)\Phi &&\approx \mathrm{i}(f^\prime / u)\Phi
\end{alignat*}
$$

Für die Konzentration folgt somit:

$$
\begin{alignat*}{5}
& |\Phi| &&= |\xi\mathrm{e}^{\mathrm{i}\lambda x}| &&\approx |\xi\mathrm{e}^{\mathrm{i}\lambda^\prime x}| &&= |\xi\mathrm{e}^{(-\Im(\lambda^\prime) + \mathrm{i}\Re(\lambda^\prime))x}| &&= |\xi\mathrm{e}^{-\Im(\lambda^\prime)x}| \\
\Rightarrow\quad& &&= |\xi\mathrm{e}^{\mathrm{i}(f / u)x}| &&\approx |\xi\mathrm{e}^{\mathrm{i}(f^\prime / u)x}| &&= |\xi\mathrm{e}^{(-\Im(f^\prime / u) + \mathrm{i}\Re(f^\prime / u))x}| &&= |\xi\mathrm{e}^{-\Im(f^\prime / u)x}|
\end{alignat*}
$$

> **Abbildung (Numerische Konvektion)**
>
> ![negative numerische Konvektion](.Dateien/Bilder/Uebertragungsverhalten_exp_pos_x.svg) ![null numerische Konvektion](.Dateien/Bilder/Uebertragungsverhalten_exp_zero_x.svg) ![positive numerische Konvektion](.Dateien/Bilder/Uebertragungsverhalten_exp_neg_x.svg)

Dabei bestimmt der Imaginärteil der modifizierten Wellenzahl bzw. modifizierten Frequenz das numerische Konvektionsverhalten, wobei die räumliche Konzentrationsänderung (als Ursache) der Konvektion ebenjener Strömungsgröße entgegenwirkt.

### Numerische Diffusion

Wenn sich die Approximation der Ableitung auf die zeitliche Änderung der Konzentration auswirkt, dann ist von numerischer Diffusion die Rede. Um diesen Effekt zu untersuchen wird für die zeitliche Ableitung

$$
x = 0
$$

gesetzt und für die räumliche Ableitung mit der Heuristik

$$
f = u\lambda \quad\Rightarrow\quad f^\prime \approx u\lambda^\prime
$$

substituiert.

$$
\begin{alignat*}{3}
& \partial_t \Phi &&= -\mathrm{i}f\Phi &&\approx -\mathrm{i}f^\prime\Phi \\
\Rightarrow\quad& &&= -\mathrm{i}(u\lambda)\Phi &&\approx -\mathrm{i}(u\lambda^\prime)\Phi
\end{alignat*}
$$

Für die Konzentration folgt somit:

$$
\begin{alignat*}{5}
& |\Phi| &&= |\xi\mathrm{e}^{-\mathrm{i}f t}| &&\approx |\xi\mathrm{e}^{-\mathrm{i}f^\prime t}| &&= |\xi\mathrm{e}^{(\Im(f^\prime) - \mathrm{i}\Re(f^\prime))t}| &&= |\xi\mathrm{e}^{\Im(f^\prime)t}| \\
\Rightarrow\quad& &&= |\xi\mathrm{e}^{-\mathrm{i}(u\lambda)t}| &&\approx |\xi\mathrm{e}^{-\mathrm{i}(u\lambda^\prime)t}| &&= |\xi\mathrm{e}^{(\Im(u\lambda^\prime) - \mathrm{i}\Re(u\lambda^\prime))t}| &&= |\xi\mathrm{e}^{\Im(u\lambda^\prime)t}|
\end{alignat*}
$$

> **Abbildung (Numerische Diffusion)**
>
> ![positive numerische Diffusion](.Dateien/Bilder/Uebertragungsverhalten_exp_neg_t.svg) ![null numerische Diffusion](.Dateien/Bilder/Uebertragungsverhalten_exp_zero_t.svg) ![negative numerische Diffusion](.Dateien/Bilder/Uebertragungsverhalten_exp_pos_t.svg)

Dabei bestimmt der Imaginärteil der modifizierten Frequenz bzw. modifizierten Wellenzahl das numerische Diffusionsverhalten, wobei die zeitliche Konzentrationsänderung (als Ursache) der Diffusion ebenjener Strömungsgröße entgegenwirkt.

### Beispiele

---
<details>
<summary markdown="span"><b>Räumliche Ableitung: Zentraldifferenz 2. Ordnung</b></summary>
<br>

$$
\begin{align*}
\partial_x \Phi(x_j,0) &\approx \frac{1}{2 h_x}\left[ \Phi(x_{j+1},0)-\Phi(x_{j-1},0) \right] \\
{\color{blue}\mathrm{i}{\color{red}\lambda} \mathrm{e}^{\mathrm{i}\lambda j h_x}} &\approx \frac{1}{2 h_x}\left[ \mathrm{e}^{\mathrm{i}\lambda (j+1) h_x}-\mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
&\approx \frac{1}{2 h_x}\left[ \mathrm{e}^{\mathrm{i}\lambda h_x}-\mathrm{e}^{-\mathrm{i}\lambda h_x} \right] \mathrm{e}^{\mathrm{i}\lambda j h_x} \\
&\approx {\color{blue}\mathrm{i}{\color{red}\underbrace{\left[\frac{\sin(\lambda h_x)}{h_x}\right]}_{\eqqcolon \lambda^\prime}} \mathrm{e}^{\mathrm{i}\lambda j h_x}}
\end{align*}
$$

Die modifizierte Wellenzahl hat keinen Imaginärteil. Es wird somit keine numerische Konvektion und auch keine numerische Diffusion verursacht.

> **Abbildung (Modifizierte Wellenzahl)**
>
> ![Modifizierte Wellenzahl](.Dateien/Bilder/Modifizierte_Wellenzahl_Zentraldifferenz_O2.svg)

</details>

---
<details>
<summary markdown="span"><b>Räumliche Ableitung: Rückwärtsdifferenz 1. Ordnung</b></summary>
<br>

$$
\begin{align*}
\partial_x \Phi(x_j,0) &\approx \frac{1}{h_x}\left[ \Phi(x_{j},0)-\Phi(x_{j-1},0) \right] \\
{\color{blue}\mathrm{i}{\color{red}\lambda} \mathrm{e}^{\mathrm{i}\lambda j h_x}} &\approx \frac{1}{h_x}\left[ \mathrm{e}^{\mathrm{i}\lambda j h_x}-\mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
&\approx \frac{1}{h_x}\left[ 1-\mathrm{e}^{-\mathrm{i}\lambda h_x} \right] \mathrm{e}^{\mathrm{i}\lambda j h_x} \\
&\approx \frac{1}{h_x}\left[ 1-\cos(\lambda h_x)+\mathrm{i}\sin(\lambda h_x) \right] \mathrm{e}^{\mathrm{i}\lambda j h_x} \\
&\approx {\color{blue}\mathrm{i}{\color{red}\underbrace{\left[ \frac{\sin(\lambda h_x)}{h_x}+\mathrm{i}\frac{\cos(\lambda h_x)-1}{h_x} \right]}_{\eqqcolon \lambda^\prime}} \mathrm{e}^{\mathrm{i}\lambda j h_x}}
\end{align*}
$$

Die modifizierte Wellenzahl hat einen negativen Imaginärteil. Die numerische Konvektion erfolgt mit negativer Richtung im Raum und die numerische Diffusion mit dem Vorzeichen der Geschwindigkeit in der Zeit.

> **Abbildung (Modifizierte Wellenzahl)**
>
> ![Modifizierte Wellenzahl](.Dateien/Bilder/Modifizierte_Wellenzahl_Rueckwaertsdifferenz_O1.svg)

</details>

---
<details>
<summary markdown="span"><b>Zeitliche Ableitung: explizites Euler-Verfahren</b></summary>
<br>

$$
\begin{align*}
\partial_t \Phi(0,t_l) &\approx \frac{1}{h_t}\left[ \Phi(0,t_{l+1})-\Phi(0,t_l) \right] \\
{\color{blue}-\mathrm{i}{\color{red}f} \mathrm{e}^{-\mathrm{i}f l h_t}} &\approx \frac{1}{h_t}\left[ \mathrm{e}^{-\mathrm{i}f (l+1) h_t}-\mathrm{e}^{-\mathrm{i}f l h_t} \right] \\
&\approx \frac{1}{h_t}\left[ \mathrm{e}^{-\mathrm{i}f h_t}-1 \right] \mathrm{e}^{-\mathrm{i}f l h_t} \\
&\approx -\mathrm{i}\left[ \frac{\mathrm{i}\mathrm{e}^{-\mathrm{i}f h_t}-\mathrm{i}}{h_t} \right] \mathrm{e}^{-\mathrm{i}f l h_t} \\
&\approx {\color{blue}-\mathrm{i}{\color{red}\underbrace{\left[ \frac{\sin(f h_t)}{h_t} + \mathrm{i}\frac{\cos(f h_t)-1}{h_t} \right]}_{\eqqcolon f^\prime}} \mathrm{e}^{-\mathrm{i}f l h_t}}
\end{align*}
$$

Die modifizierte Frequenz hat einen negativen Imaginärteil. Die numerische Diffusion erfolgt mit positiver Richtung in der Zeit und die numerische Konvektion entgegen dem Vorzeichen der Geschwindigkeit im Raum.

> **Abbildung (Modifizierte Frequenz)**
>
> ![Modifizierte Frequenz](.Dateien/Bilder/Modifizierte_Frequenz_Euler_explizit.svg)

</details>

---
<details>
<summary markdown="span"><b>Zeitliche Ableitung: implizites Euler-Verfahren</b></summary>
<br>

$$
\begin{align*}
\partial_t \Phi(0,t_l) &\approx \frac{1}{h_t}\left[ \Phi(0,t_l)-\Phi(0,t_{l-1}) \right] \\
{\color{blue}-\mathrm{i}{\color{red}f} \mathrm{e}^{-\mathrm{i}f l h_t}} &\approx \frac{1}{h_t}\left[ \mathrm{e}^{-\mathrm{i}f l h_t}-\mathrm{e}^{-\mathrm{i}f (l-1) h_t} \right] \\
&\approx \frac{1}{h_t}\left[ 1-\mathrm{e}^{\mathrm{i}f h_t} \right] \mathrm{e}^{-\mathrm{i}f l h_t} \\
&\approx -\mathrm{i}\left[ \frac{\mathrm{i}-\mathrm{i}\mathrm{e}^{\mathrm{i}f h_t}}{h_t} \right] \mathrm{e}^{-\mathrm{i}f l h_t} \\
&\approx {\color{blue}-\mathrm{i}{\color{red}\underbrace{\left[ \frac{\sin(f h_t)}{h_t} + \mathrm{i}\frac{1-\cos(f h_t)}{h_t} \right]}_{\eqqcolon f^\prime}} \mathrm{e}^{-\mathrm{i}f l h_t}}
\end{align*}
$$

Die modifizierte Frequenz hat einen positiven Imaginärteil. Die numerische Diffusion erfolgt mit negativer Richtung in der Zeit und die numerische Konvektion mit dem Vorzeichen der Geschwindigkeit im Raum.

> **Abbildung (Modifizierte Frequenz)**
>
> ![Modifizierte Frequenz](.Dateien/Bilder/Modifizierte_Frequenz_Euler_implizit.svg)

</details>


<!------------------------------------------------------------------------------
Stabilitätsanalyse nach John von Neumann
------------------------------------------------------------------------------->
## Stabilitätsanalyse nach John von Neumann

Numerische und physikalische Stabilität sind zweierlei. Erstere muss unabhängig von letzterer gewährleistet sein, sonst läuft die numerische Berechnung unter Umständen Gefahr ins Unermessliche zu entfachen, wobei das Ergebnis mit der Zeit divergiert. Als während des 2. Weltkrieges die erste Rechnertechnik zur Verfügung stand und numerische Berechnungen dadurch an Bedeutung gewannen, war die numerische Stabilität mit das erste zu bewerkstelligende Problem. Eine Stabilitätsanalyse ist von daher als Sicherheitsvorkehrung einer jeden numerischen Berechnung nach wie vor unabdingbar.

Eine numerisch stabile Simulation zeichnet sich durch Approximationsverfahren aus, welche die Konvektion (in positiver Strömungsrichtung) und die Diffusion (im positiven Zeitverlauf) begünstigen. Um diese Eigenschaft messbar zu machen, wird die eindimensionale Transportgleichung ohne Quellterm in ihre Bestandteile zerlegt und separat auf Stabilität untersucht. Dafür wird eine räumliche Fourier-Mode als Ansatz verwendet,

$$
\Phi(x_j,t_l) = \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x},
$$

sodass, nach der Stabilitätsbedingung,

$$
\left| \frac{\Phi(x_j,t_{l+1})}{\Phi(x_j,t_l)} \right| = \left| \xi \right| \overset{!}{\leq} 1,
$$

die absoluten Funktionswerte in der Zeit monoton fallen. Außerdem lässt sich nach Richard Courant, Kurt Friedrichs und Hans Lewy für jede partikuläre Gleichung eine dimensionslose Zahl definieren, welche die notwendigen Modellparameter substituiert und dadurch die Angabe einer Stabilitätsbedingung ermöglicht.

### Konvektionsgleichung

Die Konvektionsgleichung entspricht der diffusions- und quellfreien Transportgleichung. Im Eindimensionalen ist sie somit gegeben durch

$$
\partial_t \Phi(x,t) = -u\partial_x \Phi(x,t).
$$

Die CFL-Zahl für den Konvektionsterm mit der 1. räumlichen Ableitung nach _`x`_ wird wie folgt definiert:

$$
\mathrm{CFL}_x \coloneqq |u|\frac{h_t}{h_x}
$$

### Diffusionsgleichung

Die Diffusionsgleichung entspricht der konvektions- und quellfreien Transportgleichung. Im Eindimensionalen ist sie somit gegeben durch

$$
\partial_t \Phi(x,t) = c\partial_x^2 \Phi(x,t).
$$

Die CFL-Zahl für den Diffusionsterm mit der 2. räumlichen Ableitung nach _`x`_ wird wie folgt definiert:

$$
\mathrm{CFL}_{xx} \coloneqq |c|\frac{h_t}{h_x^2}
$$

### Beispiele

---
<details>
<summary markdown="span"><b>Konvektionsgleichung: explizites Euler-Verfahren mit Zentraldifferenz 2. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_l) \\
&\approx \Phi(x_j,t_l) + h_t \left[-u \partial_x \Phi(x_j,t_l) \right] \\
&\approx \Phi(x_j,t_l) - \frac{u h_t}{2 h_x} \left[ \Phi(x_{j+1},t_l)-\Phi(x_{j-1},t_l) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} - \frac{u h_t}{2 h_x} \left[ \xi^l \mathrm{e}^{\mathrm{i}\lambda (j+1) h_x}- \xi^l \mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
\xi &\approx 1 - \frac{u h_t}{2 h_x} \left[ \mathrm{e}^{\mathrm{i}\lambda h_x}-\mathrm{e}^{-\mathrm{i}\lambda h_x} \right] \\
&\approx 1 - \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot\mathrm{i}\sin(\lambda h_x) \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 - \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot\mathrm{i}\sin(\lambda h_x)\right| \\
\Rightarrow\quad \mathrm{CFL}_x &= 0
\end{align*}
$$

Die Kombination der Approximationsverfahren ist ausnahmslos instabil, da das explizite Euler-Verfahren die Konvektion entgegen der Strömungsrichtung entfacht, während die Zentraldifferenz keinen Einfluss hat.

</details>

---
<details>
<summary markdown="span"><b>Konvektionsgleichung: explizites Euler-Verfahren mit Rückwärtsdifferenz 1. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_l) \\
&\approx \Phi(x_j,t_l) + h_t \left[-u \partial_x \Phi(x_j,t_l) \right] \\
&\approx \Phi(x_j,t_l) - \frac{u h_t}{h_x} \left[ \Phi(x_{j},t_l)-\Phi(x_{j-1},t_l) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} - \frac{u h_t}{h_x} \left[ \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x}- \xi^l \mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
\xi &\approx 1 - \frac{u h_t}{h_x} \left[ 1-\mathrm{e}^{-\mathrm{i}\lambda h_x} \right] \\
&\approx 1 - \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot \left[ 1-\cos(\lambda h_x)+\mathrm{i}\sin(\lambda h_x) \right] \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 - \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot \left[ 1-\cos(\lambda h_x)+\mathrm{i}\sin(\lambda h_x) \right]\right| \\
&\geq \left[ 1 - \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot \left[ 1-\cos(\lambda h_x) \right] \right]^2 +\left[\mathrm{CFL}_x\cdot \sin(\lambda h_x) \right]^2 \\
&\geq \ldots \\
&\geq 1 + 2\left[ 1-\cos(\lambda h_x) \right] \left[ \mathrm{CFL}_x \left[ \mathrm{CFL}_x - \operatorname{sgn}(u) \right] \right] \\
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_x \begin{cases}\leq 1, &\operatorname{sgn}(u)=1 \\ =0, &\operatorname{sgn}(u)<1 \end{cases}
\end{align*}
$$

Die Kombination der Approximationsverfahren ist in Abhängigkeit von der Strömungsrichtung bedingt stabil, da die Rückwärtsdifferenz in positiver Strömungsrichtung dämpft und dabei die fehlgerichtete numerische Konvektion des expliziten Euler-Verfahrens gegenläufig kompensiert.

</details>

---
<details>
<summary markdown="span"><b>Konvektionsgleichung: implizites Euler-Verfahren mit Zentraldifferenz 2. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_{l+1}) \\
&\approx \Phi(x_j,t_l) + h_t \left[-u \partial_x \Phi(x_j,t_{l+1}) \right] \\
&\approx \Phi(x_j,t_l) - \frac{u h_t}{2 h_x} \left[ \Phi(x_{j+1},t_{l+1})-\Phi(x_{j-1},t_{l+1}) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} - \frac{u h_t}{2 h_x} \left[ \xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda (j+1) h_x}- \xi^{l+1}\mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
\xi &\approx \left[ 1 + \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot\mathrm{i}\sin(\lambda h_x) \right]^{-1} \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 + \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot\mathrm{i}\sin(\lambda h_x)\right|^{-1} \\
\Rightarrow\quad \mathrm{CFL}_x &\geq 0
\end{align*}
$$

Die Kombination der Approximationsverfahren ist unbedingt stabil, da das implizite Euler-Verfahren die Konvektion in Strömungsrichtung begünstigt, während die Zentraldifferenz keinen Einfluss hat.

</details>

---
<details>
<summary markdown="span"><b>Konvektionsgleichung: implizites Euler-Verfahren mit Rückwärtsdifferenz 1. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_{l+1}) \\
&\approx \Phi(x_j,t_l) + h_t \left[-u \partial_x \Phi(x_j,t_{l+1}) \right] \\
&\approx \Phi(x_j,t_l) - \frac{u h_t}{h_x} \left[ \Phi(x_{j},t_{l+1})-\Phi(x_{j-1},t_{l+1}) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} - \frac{u h_t}{h_x} \left[ \xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x}-\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
\xi &\approx \left[1 + \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot \left[ 1-\cos(\lambda h_x)+\mathrm{i}\sin(\lambda h_x) \right]\right]^{-1} \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 + \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot \left[ 1-\cos(\lambda h_x)+\mathrm{i}\sin(\lambda h_x) \right]\right|^{-1} \\
&\leq \left[ 1 + \operatorname{sgn}(u)\cdot \mathrm{CFL}_x \cdot \left[ 1-\cos(\lambda h_x) \right] \right]^2 +\left[\mathrm{CFL}_x\cdot \sin(\lambda h_x) \right]^2 \\
&\leq \ldots \\
&\leq 1 + 2\left[ 1-\cos(\lambda h_x) \right] \left[ \mathrm{CFL}_x \left[ \mathrm{CFL}_x + \operatorname{sgn}(u) \right] \right] \\
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_x \begin{cases}\geq 0, &\operatorname{sgn}(u)=1 \\ \leq 1, &\operatorname{sgn}(u)<1 \end{cases}
\end{align*}
$$

Die Kombination der Approximationsverfahren ist in Abhängigkeit von der Strömungsrichtung bedingt stabil, da das implizite Euler-Verfahren die Konvektion in Strömungsrichtung begünstigt, während die Rückwärtsdifferenz in positiver Strömungsrichtung dämpft und in negativer Strömungsrichtung entfacht.

</details>

---
<details>
<summary markdown="span"><b>Diffusionsgleichung: explizites Euler-Verfahren mit Zentraldifferenz 2. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_l) \\
&\approx \Phi(x_j,t_l) + h_t \left[c \partial_x^2 \Phi(x_j,t_l) \right] \\
&\approx \Phi(x_j,t_l) + \frac{c h_t}{h_x^2} \left[ \Phi(x_{j+1},t_l)-2\Phi(x_j,t_l)+\Phi(x_{j-1},t_l) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} + \frac{c h_t}{h_x^2} \left[ \xi^l \mathrm{e}^{\mathrm{i}\lambda (j+1) h_x} - 2 \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} + \xi^l \mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
\xi &\approx 1 + \frac{c h_t}{h_x^2} \left[ \mathrm{e}^{\mathrm{i}\lambda h_x} - 2 + \mathrm{e}^{-\mathrm{i}\lambda h_x} \right] \\
&\approx 1 + \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \cdot 2\left[\cos(\lambda h_x)-1\right] \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 + \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \cdot 2\left[\cos(\lambda h_x)-1\right]\right| \\
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_{xx} \begin{cases}\geq 0, &\operatorname{sgn}(c)=1 \\ = 0, &\operatorname{sgn}(c)<1 \end{cases}
\end{align*}
$$

Die Kombination der Approximationsverfahren ist für positive Diffusion stabil, da das explizite Euler-Verfahren die Diffusion in dieser Richtung begünstigt, und für negative Diffusion dementsprechend instabil, wobei die Zentraldifferenz keinen Einfluss hat.

</details>

---
<details>
<summary markdown="span"><b>Diffusionsgleichung: explizites Euler-Verfahren mit Rückwärtsdifferenz 1. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_l) \\
&\approx \Phi(x_j,t_l) + h_t \left[c \partial_x^2 \Phi(x_j,t_l) \right] \\
&\approx \Phi(x_j,t_l) + \frac{c h_t}{h_x^2} \left[ \Phi(x_j,t_l)-2\Phi(x_{j-1},t_l)+\Phi(x_{j-2},t_l) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} + \frac{c h_t}{h_x^2} \left[ \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} - 2 \xi^l \mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} + \xi^l \mathrm{e}^{\mathrm{i}\lambda (j-2) h_x} \right] \\
\xi &\approx 1 + \frac{c h_t}{h_x^2} \left[ 1 - 2 \mathrm{e}^{-\mathrm{i}\lambda h_x} + \mathrm{e}^{-2 \mathrm{i}\lambda h_x} \right] \\
&\approx 1 + \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \left[1-2\cos(\lambda h_x)+\cos(2\lambda h_x)+\mathrm{i}\left[2\sin(\lambda h_x)+\sin(2\lambda h_x)\right]\right] \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 + \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \left[1-2\cos(\lambda h_x)+\cos(2\lambda h_x)+\mathrm{i}\left[2\sin(\lambda h_x)+\sin(2\lambda h_x)\right]\right]\right| \\
&\geq \left[1 + \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \left[1-2\cos(\lambda h_x)+\cos(2\lambda h_x)\right]\right]^2+\left[\mathrm{CFL}_{xx}\left[2\sin(\lambda h_x)+\sin(2\lambda h_x)\right]\right]^2 \\
\Rightarrow\quad \mathrm{CFL}_{xx} &= 0
\end{align*}
$$

Die Kombination der Approximationsverfahren ist instabil, da die Rückwärtsdifferenz bei der zweiten Ableitung die Diffusion in entgegengesetzter Richtung entfacht und dadurch das Zeitschrittverfahren destabilisiert.

</details>

---
<details>
<summary markdown="span"><b>Diffusionsgleichung: implizites Euler-Verfahren mit Zentraldifferenz 2. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_{l+1}) \\
&\approx \Phi(x_j,t_l) + h_t \left[c \partial_x^2 \Phi(x_j,t_{l+1}) \right] \\
&\approx \Phi(x_j,t_l) + \frac{c h_t}{h_x^2} \left[ \Phi(x_{j+1},t_{l+1}) - 2\Phi(x_j,t_{l+1}) + \Phi(x_{j-1},t_{l+1}) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} + \frac{c h_t}{h_x^2} \left[ \xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda (j+1) h_x} -2 \xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} +\xi^{l+1}\mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} \right] \\
\xi &\approx \left[ 1 - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \cdot 2\left[\cos(\lambda h_x)-1\right] \right]^{-1} \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \cdot 2\left[\cos(\lambda h_x)-1\right]\right|^{-1} \\
&\leq 1 - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \cdot 4\left[\cos(\lambda h_x)-1\right] + \mathrm{CFL}_{xx}^2 \cdot 4\left[\cos(\lambda h_x)-1\right]^2 \\
0 &\leq - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} + \mathrm{CFL}_{xx}^2 \left[\cos(\lambda h_x)-1\right] \\
&\leq \mathrm{CFL}_{xx}\left[\mathrm{CFL}_{xx} \left[\cos(\lambda h_x)-1\right] - \operatorname{sgn}(c) \right] \\
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_{xx} \begin{cases}\leq 1, &\operatorname{sgn}(c)=-1 \\ = 0, &\operatorname{sgn}(c)>-1 \end{cases}
\end{align*}
$$

Die Kombination der Approximationsverfahren ist für positive Diffusion instabil, da das implizite Euler-Verfahren die Diffusion entgegen dieser Richtung entfacht, und für negative Diffusion dementsprechend bedingt stabil, wobei die Zentraldifferenz keinen Einfluss hat.

</details>

---
<details>
<summary markdown="span"><b>Diffusionsgleichung: implizites Euler-Verfahren mit Rückwärtsdifferenz 1. Ordnung</b></summary>
<br>

$$
\begin{align*}
\Phi(x_j,t_{l+1}) &\approx \Phi(x_j,t_l) + h_t \partial_t \Phi(x_j,t_{l+1}) \\
&\approx \Phi(x_j,t_l) + h_t \left[c \partial_x^2 \Phi(x_j,t_{l+1}) \right] \\
&\approx \Phi(x_j,t_l) + \frac{c h_t}{h_x^2} \left[ \Phi(x_j,t_{l+1}) - 2\Phi(x_{j-1},t_{l+1}) + \Phi(x_{j-2},t_{l+1}) \right] \\
\xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} &\approx \xi^l \mathrm{e}^{\mathrm{i}\lambda j h_x} + \frac{c h_t}{h_x^2} \left[ \xi^{l+1} \mathrm{e}^{\mathrm{i}\lambda j h_x} -2 \xi^{l+1}\mathrm{e}^{\mathrm{i}\lambda (j-1) h_x} + \xi^{l+1}\mathrm{e}^{\mathrm{i}\lambda (j-2) h_x} \right] \\
\xi &\approx \left[ 1 - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \left[1-2\cos(\lambda h_x)+\cos(2\lambda h_x)+\mathrm{i}\left[2\sin(\lambda h_x)+\sin(2\lambda h_x)\right]\right] \right]^{-1} \\\\
\forall(\lambda h_x)\in\mathbb{R}\colon\quad 1 &\geq \left|1 - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \left[1-2\cos(\lambda h_x)+\cos(2\lambda h_x)+\mathrm{i}\left[2\sin(\lambda h_x)+\sin(2\lambda h_x)\right]\right]\right|^{-1} \\
&\leq \left[1 - \operatorname{sgn}(c)\cdot \mathrm{CFL}_{xx} \left[1-2\cos(\lambda h_x)+\cos(2\lambda h_x)\right]\right]^2+\left[\mathrm{CFL}_{xx}\left[2\sin(\lambda h_x)+\sin(2\lambda h_x)\right]\right]^2 \\
\Rightarrow\quad \mathrm{CFL}_{xx} &= 0
\end{align*}
$$

Die Kombination der Approximationsverfahren ist instabil, da die Rückwärtsdifferenz bei der zweiten Ableitung die Diffusion in entgegengesetzter Richtung entfacht und dadurch das Zeitschrittverfahren destabilisiert.

</details>


<!------------------------------------------------------------------------------
Implementierung angepasster Verfahren
------------------------------------------------------------------------------->
## Implementierung angepasster Verfahren

Bei der Implementierung der Approximationsverfahren geht es neben dem Rechenaufwand vorrangig um Stabilität und Genauigkeit. Unter Anbetracht der Stabilitätsanalyse und des Übertragungsverhaltens lassen sich dafür gezielt Vorkehrungen treffen.

### Aufwind-Differenzenverfahren

Die Finite-Differenzenschemata entgegen der Strömungsrichtung anzulegen, hat nicht nur einen stabilisierenden Effekt, sondern geht auch mit einer höheren Genauigkeit einher, zumal die Informationen mit der Strömung transportiert und somit rechtzeitig abgegriffen werden. Dafür wird der Konvektionsterm überall nach dem Vorzeichen der Geschwindigkeit angepasst:

$$
\left[\operatorname{diag}(\boldsymbol{u}_-)\cdot\boldsymbol{D}_{x+}^{(1)}+\operatorname{diag}(\boldsymbol{u}_+)\cdot\boldsymbol{D}_{x-}^{(1)} + \operatorname{diag}(\boldsymbol{v}_-)\cdot\boldsymbol{D}_{y+}^{(1)}+\operatorname{diag}(\boldsymbol{v}_+)\cdot\boldsymbol{D}_{y-}^{(1)}\right] \cdot
$$

Wobei

$$
\begin{align*}
    \boldsymbol{u}_- &= \min(\boldsymbol{u},\boldsymbol{0}), \\
    \boldsymbol{u}_+ &= \max(\boldsymbol{u},\boldsymbol{0}), \\
    \boldsymbol{v}_- &= \min(\boldsymbol{v},\boldsymbol{0}), \\
    \boldsymbol{v}_+ &= \max(\boldsymbol{v},\boldsymbol{0}), \\
\end{align*}
$$

und

$$
\begin{align*}
\begin{rcases}
    \boldsymbol{D}_{x-}^{(1)} \\
    \boldsymbol{D}_{y-}^{(1)}
\end{rcases}&~
\text{Rückwärtsdifferenzen}, \\
\begin{rcases}
    \boldsymbol{D}_{x+}^{(1)} \\
    \boldsymbol{D}_{y+}^{(1)}
\end{rcases}&~
\text{Vorwärtsdifferenzen}.
\end{align*}
$$

### Kombiniertes Zeitschrittverfahren

Bei der Stabilitätsanalyse hat sich die zeitliche Entwicklung mittels expliziten Euler-Verfahrens für die Diffusionsgleichung und mittels impliziten Euler-Verfahrens für die Konvektionsgleichung als stabil erwiesen. Da sich die Wirbeltransportgleichung sowohl aus dem Konvektionsterm als auch dem Diffusionsterm zusammensetzt, ist es evident, beide Zeitschrittverfahren so zu kombinieren, dass für die gesamte Gleichung Stabilität gewährleistet ist.

$$
\omega_z(t+h_t) \approx \omega_z(t) + h_t \underbrace{\left[ \nu\nabla^2\omega_z(t)\right.}_\text{explizit} - \underbrace{\left.\boldsymbol{u}(\omega_z(t+h_t))\cdot\nabla\omega_z(t+h_t) \right]}_\text{implizit}
$$

Zusätzlich kann auch hier im Konvektionsterm das Aufwind-Differenzenverfahren eingebaut werden.

### Implizite Berechnung

Für die implizite Berechnung des nächsten Funktionswertes wird die Vorschrift erst nach null aufgelöst,

$$
f({\color{red}z},t) \coloneqq \omega_z(t) - {\color{red}z} + h_t \left[ \nu\nabla^2\omega_z(t) - \boldsymbol{u}({\color{red}z})\cdot\nabla{\color{red}z} \right]
$$

und dann das Newton-Verfahren zu jedem Zeitschritt sukzessiv angewendet:

$$
\forall t \in h_t\cdot\mathbb{N}\colon\quad f(z,t) \overset{!}{=} 0 \quad\Rightarrow\quad z \approx \omega_z(t+h_t)
$$

Der vorherige Funktionswert wird dabei jeweils als Startwert verwendet.
