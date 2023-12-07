# Teil 2: Numerische Methodik

...

[TOC]


<!------------------------------------------------------------------------------
Diskretisierung des zweidimensionalen Strömungsgebiets
------------------------------------------------------------------------------->
## Diskretisierung des zweidimensionalen Strömungsgebiets

Um ein zweidimensionales Gitter zu erstellen, werden zunächst die beiden Achsen einzeln diskretisiert und anschließend zu einer gemeinsamen Produktmenge verknüpft. Diese Produktmenge liegt dann als indizierbare Matrix vor. Da dieses Vorgehen einzig und allein der Strukturierung und Visualisierung der Daten dient, kann dabei die Reihenfolge der Achsen prinzipiell beliebig festgelegt werden. Allerdings gibt es im Zweidimensionalen eine etwas sonderbare Konvention: Denn es wird zuerst die y-Achse und dann die x-Achse indiziert – also nicht chronologisch, wie es vielleicht zu vermuten wäre. Der Grund dafür ist, dass man sich die so entstehende Matrix in einem Koordinatensystem visualisiert. Dabei werden die Zeilen vertikal entlang der y-Achse und die Spalten horizontal entlang der x-Achse dargestellt.

$$
    \Omega \coloneqq Y \times X = \{(y,x) | y\in Y, x\in X\}
$$

> **Begleitmaterial (2d Gitter)**
>
> [![Matlab Beispiel](../Dateien/Plaketten/Matlab_Beispiel.svg)](../Begleitmaterial/2d_Gitter_Matlab.ipynb) [![Python Beispiel](../Dateien/Plaketten/Python_Beispiel.svg)](../Begleitmaterial/2d_Gitter_Python.ipynb)
>
> ![2d Gitter](../Dateien/Bilder/2d_Gitter.svg)
