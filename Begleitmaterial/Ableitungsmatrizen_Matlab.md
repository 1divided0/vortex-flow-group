<!------------------------------------------------------------------------------
Erstellung von Ableitungsmatrizen in Matlab
------------------------------------------------------------------------------->
# Erstellung von Ableitungsmatrizen in Matlab

Hier wird ausführlich erläutert, wie die Ableitungsmatrizen erstellt werden können. Um die Poisson-Gleichung im mehrdimensionalen Fall als lineares Gleichungssystem lösen zu können, bietet es sich an, der Konvention der Programmiersprache entsprechend zu vektorisieren.

[TOC]

## Eindimensionaler Fall

Sind die Werte der Funktion und die Koeffizienten des Differenzenschemas an den _`N`_ Stützstellen bekannt, so kann die _`k`_-te Ableitung der Funktion an diesen Stützstellen approximiert werden. Dafür lässt sich ein Koeffizientenvektor definieren, welcher bereits durch die Schrittweitenpotenz dividiert ist.

$$
\forall i\in\{1,\ldots,N\}\colon\quad (\boldsymbol{d}_{x_i}^{(k)})_{1\times{N}} \coloneqq \frac{1}{h_x^k} \begin{bmatrix}\boldsymbol{0}_{1\times{i-l-1}} & \boldsymbol{\alpha}_{1\times{l+r+1}}^{(k)} & \boldsymbol{0}_{1\times{N-i-r}}\end{bmatrix}_{x_i}
$$

Im Eindimensionalen entspricht dann die Approximation der Ableitung ebenjenen Differenzenschema, angewendet auf die einzelnen Stützstellen.

$$
\partial_x^k \boldsymbol{\Phi}_{N\times{1}} \approx \begin{bmatrix}\rule[.5ex]{2.5ex}{0.5pt}&(\boldsymbol{d}_{x_1}^{(k)})_{1\times{N}}&\rule[.5ex]{2.5ex}{0.5pt}\\&\vdots&\\\rule[.5ex]{2.5ex}{0.5pt}&(\boldsymbol{d}_{x_N}^{(k)})_{1\times{N}}&\rule[.5ex]{2.5ex}{0.5pt}\end{bmatrix} \cdot \begin{bmatrix}\Phi_{x_1}\\\vdots\\\Phi_{x_N}\end{bmatrix} \eqqcolon (\boldsymbol{D}_x^{(k)})_{N\times{N}} \cdot \boldsymbol{\Phi}_{N\times{1}}
$$

## Zweidimensionaler Fall

Im Zweidimensionalen wird das Skalarfeld ebenso über ein Gitter diskretisiert, nur dass diesmal die Werte zunächst eine Matrix und keinen Vektor bilden. Für die Handhabung mittels Finiter-Differenzen-Methode ist das jedoch etwas unpraktisch, da so für die partiellen Ableitungen zwei unterschiedliche Dualräume entstehen.

Wenn die Zeilen den diskreten _`y`_-Werten und die Spalten den diskreten _`x`_-Werten entsprechen, dann wird das Skalarfeld für die _`k`_-te partielle Ableitung nach _`x`_ in den Zeilenraum der Differenzenmatrix abgebildet. Dabei ist zu beachten, dass diese Matrix im Vergleich zum Eindimensionalen transponiert ist.

$$
\partial_x^k \boldsymbol{\Phi}_{n\times{m}} \approx \begin{bmatrix}\rule[.5ex]{2.5ex}{0.5pt}&(\boldsymbol{\Phi}_{y_1})_{1\times{m}}&\rule[.5ex]{2.5ex}{0.5pt}\\&\vdots&\\\rule[.5ex]{2.5ex}{0.5pt}&(\boldsymbol{\Phi}_{y_n})_{1\times{m}}&\rule[.5ex]{2.5ex}{0.5pt}\end{bmatrix} \cdot \begin{bmatrix}\rule[-1ex]{0.5pt}{2.5ex}&&\rule[-1ex]{0.5pt}{2.5ex}\\(\boldsymbol{d}_{x_1}^{(k)})_{m\times{1}}&\cdots&(\boldsymbol{d}_{x_m}^{(k)})_{m\times{1}}\\\rule[-1ex]{0.5pt}{2.5ex}&&\rule[-1ex]{0.5pt}{2.5ex}\end{bmatrix} \eqqcolon \boldsymbol{\Phi}_{n\times{m}} \cdot (\boldsymbol{D}_x^{(k)})_{m\times{m}}^\top
$$

Für die _`k`_-te partielle Ableitung nach _`y`_ erfolgt die Abbildung dementsprechend in den Spaltenraum der Differenzenmatrix, so wie es auch im Eindimensionalen bewerkstelligt wurde.

$$
\partial_y^k \boldsymbol{\Phi}_{n\times{m}} \approx \begin{bmatrix}\rule[.5ex]{2.5ex}{0.5pt}&(\boldsymbol{d}_{y_1}^{(k)})_{1\times{n}}&\rule[.5ex]{2.5ex}{0.5pt}\\&\vdots&\\\rule[.5ex]{2.5ex}{0.5pt}&(\boldsymbol{d}_{y_n}^{(k)})_{1\times{n}}&\rule[.5ex]{2.5ex}{0.5pt}\end{bmatrix} \cdot \begin{bmatrix}\rule[-1ex]{0.5pt}{2.5ex}&&\rule[-1ex]{0.5pt}{2.5ex}\\(\boldsymbol{\Phi}_{x_1})_{n\times{1}}&\cdots&(\boldsymbol{\Phi}_{x_m})_{n\times{1}}\\\rule[-1ex]{0.5pt}{2.5ex}&&\rule[-1ex]{0.5pt}{2.5ex}\end{bmatrix} \eqqcolon (\boldsymbol{D}_y^{(k)})_{n\times{n}} \cdot \boldsymbol{\Phi}_{n\times{m}}
$$

Es ist anzumerken, dass bei dieser Approximation der partiellen Ableitungen immer die selbe Differenzenmatrix auf alle Stützstellen angewendet wird. Außerdem stellt sich heraus, dass die Rechnung sehr viel übersichtlicher wird, wenn das zweidimensionale Gitter vektorisiert ist. So können die beiden partiellen Ableitungen in ein und denselben Vektorraum abgebildet werden.

Werden nun also bei der Vektorisierung die Spalten verkettet – so wie in Matlab üblich – dann lassen sich die partiellen Ableitungen mittels Kronecker-Produkt ([`kron`](https://de.mathworks.com/help/matlab/ref/kron.html)) wie folgt schreiben:

$$
\begin{align*}
\partial_x^k \boldsymbol{\Phi}_{mn\times{1}} &\approx \underbrace{\left[(\boldsymbol{D}_x^{(k)})_{m\times{m}} \otimes \boldsymbol{I}_{n\times{n}}\right]}_{\eqqcolon(\boldsymbol{D}_x^{(k)})_{mn\times{mn}}} \cdot \boldsymbol{\Phi}_{mn\times{1}} \\\\
\partial_y^k \boldsymbol{\Phi}_{mn\times{1}} &\approx \underbrace{\left[\boldsymbol{I}_{m\times{m}} \otimes (\boldsymbol{D}_y^{(k)})_{n\times{n}}\right]}_{\eqqcolon(\boldsymbol{D}_y^{(k)})_{mn\times{mn}}} \cdot \boldsymbol{\Phi}_{mn\times{1}}
\end{align*}
$$

Die partiellen Ableitungen werden dabei, wie zuvor, jeweils über die Stützstellen berechnet, an denen sich die andere Koordinate nicht verändert, und es wird immer die selbe Differenzenmatrix angewendet. Wenn die Topologie des Strömungsgebiets jedoch komplizierter ist, dann reicht es womöglich nicht mehr aus alle Stützstellen gleich zu behandeln, sodass der Ausdruck mit dem Kronecker-Produkt individuell auf das Strömungsgebiet angepasst werden muss, damit das Differenzenschema auch auf den Rändern des Strömungsgebiets konsistent ist.

Unter Umständen ist das Strömungsgebiet nicht einfach zusammenhängend und enthält beispielsweise ein Hindernis. Die nachfolgende Abbildung soll diesen Sachverhalt veranschaulichen, wobei die Auflösung für den Demonstrationszweck reduziert ist und für eine praktikable Anwendung eigentlich erhöht werden müsste, damit zwischen den Rändern des Strömungsgebiets genügend Platz für konsistente Differenzenschemata vorhanden ist.

> **Abbildung (2d Ableitung mit Hindernis)**
>
> ![2d Ableitung mit Hindernis](Dateien/Bilder/2d_Ableitung_mit_Hindernis_Matlab.svg)

In diesem Fall müssen die Ableitungsmatrizen aus individuellen Blöcken zusammengesetzt werden.

$$
\forall{i,j}\in\{1,\ldots,5\}\colon\quad (\boldsymbol{D}_{x}^{(k)})_{ij} \coloneqq \begin{bmatrix}(\boldsymbol{D}_{x_1}^{(k)})_{ij}&&\\&(\boldsymbol{D}_{x_2}^{(k)})_{ij}&\\&&(\boldsymbol{D}_{x_3}^{(k)})_{ij}\end{bmatrix}
$$

$$
\boldsymbol{D}_x^{(k)} \coloneqq \begin{bmatrix}(\boldsymbol{D}_{x}^{(k)})_{11}&(\boldsymbol{D}_{x}^{(k)})_{12}&(\boldsymbol{D}_{x}^{(k)})_{13}&(\boldsymbol{D}_{x}^{(k)})_{14}&(\boldsymbol{D}_{x}^{(k)})_{15}\\(\boldsymbol{D}_{x}^{(k)})_{21}&(\boldsymbol{D}_{x}^{(k)})_{22}&(\boldsymbol{D}_{x}^{(k)})_{23}&(\boldsymbol{D}_{x}^{(k)})_{24}&(\boldsymbol{D}_{x}^{(k)})_{25}\\(\boldsymbol{D}_{x}^{(k)})_{31}&(\boldsymbol{D}_{x}^{(k)})_{32}&(\boldsymbol{D}_{x}^{(k)})_{33}&(\boldsymbol{D}_{x}^{(k)})_{34}&(\boldsymbol{D}_{x}^{(k)})_{35}\\(\boldsymbol{D}_{x}^{(k)})_{41}&(\boldsymbol{D}_{x}^{(k)})_{42}&(\boldsymbol{D}_{x}^{(k)})_{43}&(\boldsymbol{D}_{x}^{(k)})_{44}&(\boldsymbol{D}_{x}^{(k)})_{45}\\(\boldsymbol{D}_{x}^{(k)})_{51}&(\boldsymbol{D}_{x}^{(k)})_{52}&(\boldsymbol{D}_{x}^{(k)})_{53}&(\boldsymbol{D}_{x}^{(k)})_{54}&(\boldsymbol{D}_{x}^{(k)})_{55}\end{bmatrix},\quad
\boldsymbol{D}_y^{(k)} \coloneqq \begin{bmatrix}\boldsymbol{D}_{y_1}^{(k)}&&&&\\&\boldsymbol{D}_{y_2}^{(k)}&&&\\&&\boldsymbol{D}_{y_3}^{(k)}&&\\&&&\boldsymbol{D}_{y_4}^{(k)}&\\&&&&\boldsymbol{D}_{y_5}^{(k)}\end{bmatrix}
$$
