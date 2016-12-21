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
This repository contains ipython notebook to create every figure and table of the ESI and the manuscript.


### Figures_Article
This repository contains ipython notebook to create every figure from the ESI from S1 to S16.
And you can find the code of other algorithms TRAIn, ITAMeD and ITAMeD with lp in the following links:

- TRAIn:  https://code.google.com/archive/p/diffusion-mri/source
- ITAMeD: http://nmr.cent.uw.edu.pl/index.php/download/category/1-diffusion-nmr 
- ITAMeD with lp: http://nmr.cent.uw.edu.pl/index.php/download/category/1-diffusion-nmr

### Figures_ESI
This repostory contains the ipython notebook code of table S2 and table S3 in the ESI.
All signals are available in this file.

### Code
The Code repository contains the code of PALMA algorithm and the code of all simulated signals presented in the Support Information.
These files will be used in all other tests through an alias, and a symbolic link to this folder should be created in each location where PALMA is to be used.
The `link.sh` script creates these links in all the required location for these tests.
So please run this script before running the different examples.


---
*This Program are provided under the Licence CeCILL 2.1*

*This program HAS NOT been tested intensively, it is believed to do what it is supposed to do, However, you are welcome to check it on your own data.*

    Author : Afef Cherni, Marc-André Delsuc (madelsuc@unistra.fr)
    Version : 1.0   Date : June 2016
    Version : 2.0   Date : December 2016
