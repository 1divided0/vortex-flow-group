<!------------------------------------------------------------------------------
Beispiele für das numerische Übertragungsverhalten
------------------------------------------------------------------------------->
# Beispiele für das numerische Übertragungsverhalten

Hier befinden sich Beispielrechnungen zum numerischen Übertragungsverhalten.

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
