# Thema 3: Numerische Methodik

_... befindet sich in Arbeit ..._

[TOC]


<!------------------------------------------------------------------------------
Räumliche Diskretisierung
------------------------------------------------------------------------------->
## Räumliche Diskretisierung

Um ein zweidimensionales Gitter zu erstellen, werden zunächst die beiden Achsen einzeln diskretisiert und anschließend zu einer gemeinsamen Produktmenge verknüpft.

$$
\begin{align*}
    \Omega &\coloneqq X \times Y = \{(x,y) ~|~ x\in X, y\in Y\} \\
    &\Downarrow~\text{Diskretisierung} \\
    \boldsymbol{\Omega}_{ji} &~\widehat{=}~ (x_i,y_j)
\end{align*}
$$

Diese Produktmenge liegt dann als indizierbare Matrix vor. Die Reihenfolge der Achsen kann dabei prinzipiell beliebig festgelegt werden. Allerdings gibt es im Zweidimensionalen eine etwas sonderbare Konvention: Denn es wird zuerst die y-Achse und dann die x-Achse indiziert – also nicht chronologisch, wie es vielleicht zu vermuten wäre. Der Grund dafür ist, dass man sich die so entstehende Matrix in einem Koordinatensystem vorstellt.

---
> **Begleitmaterial (2d Gitter)**
>
> [![Matlab Beispiel](.Dateien/Plaketten/Matlab_Beispiel.svg)](Begleitmaterial/2d_Gitter_Matlab.ipynb) [![Python Beispiel](.Dateien/Plaketten/Python_Beispiel.svg)](Begleitmaterial/2d_Gitter_Python.ipynb)
>
> ![2d Gitter](.Dateien/Bilder/2d_Gitter.svg)


<!------------------------------------------------------------------------------
Vektorisierung des Rechengitters
------------------------------------------------------------------------------->
## Vektorisierung des Rechengitters

Auf dem erstellten Gitter werden die strömungsmechanischen Gleichungen numerisch gelöst. Dies geschieht bei der Finiten-Differenzen-Methode in Vektorform, damit sich die Rechenoperationen über Matrizen abwickeln lassen. Um die Gleichungen in dieser Form lösen zu können, müssen die auf dem Gitter abgespeicherten Werte von einer Matrix in einen Vektor überführt werden. Je nachdem werden dafür entweder die Spalten oder Zeilen verkettet.

---
> **Begleitmaterial (Vektorisierung)**
>
> [![Matlab Beispiel](.Dateien/Plaketten/Matlab_Beispiel.svg)](Begleitmaterial/Vektorisierung_Matlab.ipynb) [![Python Beispiel](.Dateien/Plaketten/Python_Beispiel.svg)](Begleitmaterial/Vektorisierung_Python.ipynb)
>
> ![2d Gitter](.Dateien/Bilder/Verkettung.svg)


<!------------------------------------------------------------------------------
Behandlung dünnbesetzter Matrizen
------------------------------------------------------------------------------->
## Behandlung dünnbesetzter Matrizen

Nach der Erstellung eines _`n×m`_-Gitters werden dessen Werte je nach Art der Vektorisierung über _`mn×mn`_- bzw. _`nm×nm`_-Matrizen abgebildet. Diese Matrizen sind bei der Finiten-Differenzen-Methode sehr groß und dünnbesetzt. Um die Laufzeit des Lösungsalgorithmus zu verbessern, macht es durchaus Sinn, diese dünnbesetzten Matrizen platzsparend zu speichern. Standardbibliotheken stellen dafür spezielle Datenstrukturen bereit.

In Matlab basiert die Implementierung für dünnbesetzte Matrizen auf dem CSC-Format:
- John Gilbert, Cleve Moler, Robert Schreiber. Sparse Matrices in MATLAB - Design and Implementation. <https://www.mathworks.com/help/pdf_doc/otherdocs/simax.pdf>. (1992)

In Python bietet [`scipy.sparse`](https://docs.scipy.org/doc/scipy/reference/sparse.html) ein Sammelsurium von unterschiedlichen Formaten (wie z. B. CSR, CSC und COO):
- SciPy Dokumentation. Sparse Arrays. <https://docs.scipy.org/doc/scipy/tutorial/sparse.html>

Weitere Informationen über die möglichen Formate dünnbesetzter Matrizen finden sich auch auf Wikipedia:
- Wikipedia Artikel. Sparse matrix. <https://en.wikipedia.org/wiki/Sparse_matrix>

---
<details>
<summary markdown="span"><b>Vergleich der Laufzeit</b></summary>
<br>
Hier lässt sich ein Vergleich der Laufzeit anstellen. Dafür können zwei Verfahren verwendet werden, die auch bei der algorithmischen Auslegung von Computerclustern eingesetzt werden:
- Bei dem sog. Strong-Scaling wird die Laufzeit für eine fixe Problemgröße mit zunehmender Anzahl von Prozessorkernen gemessen. Dieser Zusammenhang wird durch das Gesetz von Amdahl beschrieben.
- Bei dem sog. Weak-Scaling wird die Laufzeit für eine fixe Anzahl von Prozessorkernen mit zunehmender Problemgröße gemessen. Dieser Zusammenhang wird durch das Gesetz von Gustafson beschrieben.

In dieser Projektarbeit wird der Programmcode nicht parallelisiert, dennoch ist die Laufzeit gerade bei Skriptsprachen ein berechtigtes Bedenken. Denn im Gegensatz zur AOT-Kompilierung erfolgt die bei den Skriptsprachen verwendete JIT-Kompilierung während der Laufzeit, was die algorithmische Ausführung dementsprechend verlangsamt.

Als Beispiel für das Weak-Scaling wird hier die Laufzeit bei der Multiplikation einer Diagonalmatrix mit einem Vektor unter Verwendung unterschiedlicher Speicherformate gemessen. Das Ergebnis ist in jedem Fall das gleiche, nur die Laufzeit unterscheidet sich wesentlich.

$$
\underbrace{\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}}_\text{Rechenoperation mit Nullen} = \underbrace{\mathbf{1}_{N \times 1}\odot\mathbf{1}_{N \times 1}}_\text{Rechenoperation ohne Nullen} = \mathbf{1}_{N \times 1}
$$
</details>

---
> **Begleitmaterial (Vergleich der Laufzeit)**
>
> [![Matlab Beispiel](.Dateien/Plaketten/Matlab_Beispiel.svg)](Begleitmaterial/Sparse_Performance.m) [![Python Beispiel](.Dateien/Plaketten/Python_Beispiel.svg)](Begleitmaterial/Sparse_Performance.py)
>
> ![Vergleich der Laufzeit](.Dateien/Bilder/Sparse_Performance.svg)


<!------------------------------------------------------------------------------
Allgemeine Finite-Differenzen-Methode
------------------------------------------------------------------------------->
## Allgemeine Finite-Differenzen-Methode

Das Ziel einer jeden finiten Differenz ist es – basiert auf einem Schema, welches auf bestimmte Stützstellen vor- und zurückgreift – eine Ableitung zu approximieren. Dafür wird der Wert der abzuleitenden Funktion _`Φ`_ an der jeweiligen Stelle, sowie den _`l`_ Stellen links und den _`r`_ Stellen rechts davon durch Koeffizienten _`α`_ gewichtet und anschließend aufsummiert. Die Koeffizienten des Differenzenschemas leiten sich dabei aus der Taylor-Reihenentwicklung ab. Ist die Schrittweite _`h`_ äquidistant, dann lässt sich ihre entsprechende Potenz als Reziproke aus dem Differenzenschema faktorisieren.

$$
\partial_x^k \Phi_i = \frac{1}{h_x^k} \sum_{j=-l}^{r} \alpha_j \Phi_{i+j}
$$

> **Abbildung (Differenzenschema)**
>
> ![Differenzenschema](.Dateien/Bilder/Differenzenschema.svg)

Die Taylor-Reihenentwicklung ist somit gegeben durch

$$
\Phi_{i+j} = \left[ \sum_{p=0}^{l+r} (\partial_x^p \Phi_i) (j h_x)^p / p! \right] + \mathcal{O}(h_x^{l+r+1}),
$$

wobei bis zur Ordnung _`l+r+1`_ entwickelt wird, da dies der Anzahl von Stützstellen entspricht und damit auch der Größe des, für die Koeffizienten _`α`_, zu lösenden Gleichungssystems. Der Ausdruck für die Reihenentwicklung lässt sich sodann in das Differenzenschema einsetzen.

$$
\partial_x^k \Phi_i = \frac{1}{h_x^k} \sum_{j=-l}^{r} \alpha_j \left\{ \left[ \sum_{p=0}^{l+r} (\partial_x^p \Phi_i) (j h_x)^p / p! \right] + \mathcal{O}(h_x^{l+r+1}) \right\}
$$

Die Summen können vertauscht werden und der Ordnungsterm reduziert sich durch das Reziproke der Schrittweitenpotenz.

$$
\partial_x^k \Phi_i = \frac{1}{h_x^k} \left[ \sum_{p=0}^{l+r} \sum_{j=-l}^{r} \alpha_j (\partial_x^p \Phi_i) (j h_x)^p / p! \right] + \mathcal{O}(h_x^{l+r+1-k})
$$

Um aus der Vorschrift ein Gleichungssystem abzuleiten, wird die äußere Summe mit _`p=k`_ und _`p≠k`_ aufgespalten.

$$
\partial_x^k \Phi_i = \underbrace{\left[ \sum_{j=-l}^{r} j^k \alpha_j / k! \right]}_{\stackrel{!}{=}1} \partial_x^k \Phi_i
+ \sum_{\substack{p = 0 \\ p \ne k}}^{l+r} h_x^{p-k} \underbrace{\left[ \sum_{j=-l}^{r} j^p \alpha_j / p! \right]}_{\stackrel{!}{=}0} \partial_x^p \Phi_i
+ \mathcal{O}(h_x^{l+r+1-k})
$$

Damit die Vorschrift konsistent ist, muss der Ausdruck in der ersten Klammer 1 und in der zweiten Klammer 0 ergeben, sodass sich für die Koeffizienten folgendes Gleichungssystem ergibt:

$$
\forall p \in \{ 0,\ldots,l+r \} \ni k\colon\quad \sum_{j=-l}^{r} j^p \alpha_j = \delta_{kp} p!
$$

Daraus lässt sich ein beliebiges Differenzenschema für die _`k`_-te Ableitung mit der Fehlerordnung _`l+r+1-k`_ konstruieren. Die Ableitungsordnung muss dabei nur kleiner sein als die Anzahl von Stützstellen.

---
> **Aufgabe (Allgemeines Differenzenschema)**
>
> Schreibt ein Programm, welches mit _`l≥0`_, _`r≥0`_ und _`0<k≤l+r`_ die Koeffizienten des Differenzenschemas für die _`k`_-te Ableitung berechnet.


<!------------------------------------------------------------------------------
Erstellung eindimensionaler Ableitungsmatrizen
------------------------------------------------------------------------------->
## Erstellung eindimensionaler Ableitungsmatrizen

Sind die Werte der Funktion und die Koeffizienten des Differenzenschemas an den Stützstellen bekannt, so kann die Ableitung der Funktion approximiert werden. Im Eindimensionalen entspricht diese Ableitung ebenjenen Differenzenschema, angewendet auf die einzelnen Stützstellen.

$$
\partial_x^k \boldsymbol{\Phi}_{N\times{1}} \approx \begin{bmatrix}\rule[.5ex]{2.5ex}{0.5pt}\boldsymbol{d}_1^{(k)}\rule[.5ex]{2.5ex}{0.5pt}\\\vdots\\\rule[.5ex]{2.5ex}{0.5pt}\boldsymbol{d}_N^{(k)}\rule[.5ex]{2.5ex}{0.5pt}\end{bmatrix} \begin{bmatrix}\Phi_1\\\vdots\\\Phi_N\end{bmatrix} \eqqcolon \boldsymbol{D}_{N\times{N}}^{(k)} \boldsymbol{\Phi}_{N\times{1}}
$$

---
> **Aufgabe (Ableitung 1d)**
>
> Schreibt ein Programm, welches die erste Ableitung von _`sin(x)`_ mit der Zentraldifferenz _`l=r=1`_ berechnet und vergleicht das Ergebnis mit der analytischen Lösung. Welche Fehlerordnung hat dieses Ableitungsverfahren?


<!------------------------------------------------------------------------------
Erstellung zweidimensionaler Ableitungsmatrizen
------------------------------------------------------------------------------->
## Erstellung zweidimensionaler Ableitungsmatrizen

...
