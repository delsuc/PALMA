# PALMA




*This document is associated to the manuscript*

**PALMA, an improved algorithm for DOSY signal processing**

*Afef Cherni Emilie Chouzenoux and Marc-André Delsuc¨*

[arxiv.org:1608.07055](https://arxiv.org/abs/1608.07055)

## Abstract
NMR is a tool of choice for the measure of diffusion coefficients of species in solution.
The DOSY experiment, a 2D implementation of this measure, has proven to be particularly useful for the study of complex mixtures, molecular interactions, polymers, etc.
However, DOSY data analysis requires to resort to inverse Laplace transform, in particular for polydisperse samples.
This is a known difficult numerical task, for which we present here a novel approach.
A new algorithm based on a splitting scheme and on the use of proximity operators is introduced.
Used in conjunction with a Maximum Entropy and $\ell_1$  hybrid regularisation, this algorithm converges rapidly and produces results robust against experimental noise.
This method has been called PALMA.
It is able to reproduce faithfully monodisperse as well as polydisperse systems, and numerous simulated and experimental examples are presented.
It has been implemented on the server http://palma.labo.igbmc.fr where users can have their datasets processed automatically.

## This repository
This repository contains the python and jpyter code used to generate the data and the figures in the manuscript


---
*This Program are provided under the Licence CeCILL 2.1*

*This program HAS NOT been tested intensively, it is believed to do what it is supposed to do, However, you are welcome to check it on your own data.*

    Author : M-A Delsuc (madelsuc@unistra.fr)
    Date : June 2016
    Version : 1.0
