# QValue $2\nu\beta\beta$ 
Author: Toby Dixon (toby.dixon.23@ucl.ac.uk)
A tool to fit the $2\nu\beta\beta$ decay spectrum treating $Q_{\beta\beta}$ as a free parameter.

### Mathematical formulation
The idea of this analysis is to extract the $Q_{\beta\beta}$ value for $2\nu\beta\beta$ decay directly from the $2\nu\beta\beta$ spectrum. This is described in https://arxiv.org/pdf/2210.08722.
The $2\nu\beta\beta$ spectrum can be described as:

$$
    \frac{d\Gamma}{dE}=A(E)\times (Q-E)^5
$$

Where $A(E)$ is the 'shape' factor, fairly complicated to compute depending on the atomic and nuclear physics of the decay.
The basic principle of our analysis is to perform a fit to extract an analytic approximation / interpolation of A(E), the $2\nu\beta\beta$ decay distribution can then be included in a background model fit with $Q_{\beta\beta}$ as a free parameter.
Within the improved description (https://arxiv.org/pdf/2307.14086) of $2\nu\beta\beta$ decay:

$$
    \frac{d\Gamma}{dE}=g_{A,eff}^4|M_{GT-1}|^2\Big[\frac{dG_0}{dE}+\xi_{3,1}\frac{dG_2}{dE}+...]
$$

Or alternatively to allow the parameterisation above:

$$
    \frac{d\Gamma}{dE}=g_{A,eff}^4|M_{GT-1}|^2(Q-E)^5\Big[A_0(E)+\xi_{3,1}A_2(E)+...]
$$

Thus we can parameterise independently each of the additive contributions to the spectrum, for inclusion in a fit considering both the spectral shape uncertainty and $Q_{\beta\beta}$ uncertainities.

### BetaSpectrumHandler
The first part of the analysis is contained in a script called `BetaSpectrumHandler.cxx`, which performs this parameterisation.
This code takes in a MC simulation of the $2\nu\beta\beta$ decay obtained with the $Q=3034.4$ keV.
It then fits this spectrum to $A(E)\times (Q-E)^5$, whether $A(E)$ is parameterised as a 6th order polynomial.
The code to run this can be compiled with:

    g++ -std=c++11 -g -o BetaSpectrumHandler BetaSpectrumHandler.cxx $(root-config --cflags --libs)  

It can then be run with:

    ./BetaSpectrumHandler -l [plot_path] -i [input_path] -p [polynomial order] -q [Qbb] -o [outpath] 

add the option -I to use the improved model and -h for this help

### QbbFit
The next part of the analysis chain consists of a BAT code called QbbFit which performs a background model fit on the spectrum using the analytic approximation of the $2\nu\beta\beta$ decay spectrum shape.

