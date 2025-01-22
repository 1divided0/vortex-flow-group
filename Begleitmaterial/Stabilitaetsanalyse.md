<!------------------------------------------------------------------------------
Beispiele für die Stabilitätsanalyse 
------------------------------------------------------------------------------->
# Beispiele für die Stabilitätsanalyse 

Hier befinden sich Beispielrechnungen zur Stabilitätsanalyse.

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
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_x \begin{cases}\leq 1, &\operatorname{sgn}(u)=1 \\ =0, &\operatorname{sgn}(u)<1 \end{cases}
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
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_x \begin{cases}\geq 0, &\operatorname{sgn}(u)=1 \\ \leq 1, &\operatorname{sgn}(u)<1 \end{cases}
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
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_{xx} \begin{cases}\geq 0, &\operatorname{sgn}(c)=1 \\ = 0, &\operatorname{sgn}(c)<1 \end{cases}
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
&\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\! \Rightarrow\quad \mathrm{CFL}_{xx} \begin{cases}\leq 1, &\operatorname{sgn}(c)=-1 \\ = 0, &\operatorname{sgn}(c)>-1 \end{cases}
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
