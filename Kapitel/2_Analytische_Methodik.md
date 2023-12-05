# Teil 2: Analytische Methodik

...

[TOC]


<!------------------------------------------------------------------------------
Wirbeltransportgleichung
------------------------------------------------------------------------------->
## Wirbeltransportgleichung

Um Wirbelströmungen adäquat untersuchen zu können, wird der Transport einer ganz bestimmten Strömungsgröße betrachtet – nämlich der Wirbelstärke $\boldsymbol{\omega}$. Sie ist wie folgt definiert:

$$
    \boldsymbol{\omega} \coloneqq \operatorname{rot}(\boldsymbol{u}) = \nabla\times\boldsymbol{u}
$$

Für zweidimensionale Strömungen in der x-y-Ebene ist die Geschwindigkeit in z-Richtung immer null,

$$ u_z=0 $$

und die Wirbelstärke hat demzufolge nur eine Komponente in z-Richtung.

$$
    \omega_z = \partial_x u_y - \partial_y u_x
$$

Für die Herleitung der inkompressiblen Wirbeltransportgleichung

$$
    \partial_t\boldsymbol{\omega} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{\omega} = (\boldsymbol{\omega}\cdot\nabla)\boldsymbol{u} + \nu\nabla^2\boldsymbol{\omega}
$$

wird die Rotation auf die inkompressible Navier-Stokes-Gleichung angewendet:

$$
    \nabla\times \left\{ \partial_t\boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u} \right\} = \nabla\times \left\{ -\nabla{p}/\rho + \nu\nabla^2\boldsymbol{u} + \boldsymbol{g} \right\}
$$

> **Aufgabe (Herleitung der Wirbeltransportgleichung)**
>
> Leitet die Wirbeltransportgleichung für inkompressible zweidimensionale Strömungen her. Verwendet dazu die folgenden Zusammenhänge.
> 
> $$ \begin{align*} (\boldsymbol{u}\cdot\nabla)\boldsymbol{u} &=\nabla\boldsymbol{u}^2/2 -\boldsymbol{u}\times\boldsymbol{\omega} \\ \nabla\times(\boldsymbol{u}\times\boldsymbol{\omega}) &= \boldsymbol{u}(\nabla\cdot\boldsymbol{\omega}) - \boldsymbol{\omega}(\nabla\cdot\boldsymbol{u}) + (\boldsymbol{\omega}\cdot\nabla)\boldsymbol{u} - (\boldsymbol{u}\cdot\nabla)\boldsymbol{\omega} \\ \nabla\cdot\boldsymbol{\omega} &=0 \quad\text{(Die Wirbelstärke erfüllt die Kontinuitätsgleichung)} \\ \nabla\times\nabla\phi &=0 \quad\text{(Gradientenfelder sind wirbelfrei)} \end{align*} $$
> 
> **Tipp:** Die Schwerkraft ist eine konservative Kraft und lässt sich somit als Gradient der potentiellen Energie schreiben.
