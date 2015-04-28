Source code for a R package that fits and the seven models presented in Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at http://arxiv.org/abs/1406.4355

The main goal of the package is to compare the CCRW to the TLW. The package also allows to simulate some of these models. The files in the folder Simulations are R files to reproduce the simulations described in the manuscript.

To be able to use the package you need to build it from this source code. If unfamiliar with compiling code, I recommend using RStudio (http://www.rstudio.com/). See the information found on RStudio support website to get help on how to build a package with Rstudio: https://support.rstudio.com/hc/en-us/articles/200486488-Developing-Packages-with-RStudio. Depending on your operating system you will need compilers, for more information see: https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites

You can also use devtools to install the package:
> install_github(‘MarieAugerMethe/CCRWvsLW’)
