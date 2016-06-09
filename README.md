Source code for a R package that fits and the seven models presented in:
Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution 6:1179-1189. Preprint available at http://arxiv.org/abs/1406.4355
and the new model extensions presented in:
Auger-Methe, M., A.E. Derocher, C.A. DeMars, M.J. Plank, E.A. Codling, M.A. Lewis (2016-InPress) Evaluating random search strategies in three mammals from distinct feeding guilds. Journal of Animal Ecology


The main goal of the package is to compare the CCRW to the TLW. The package also allows to simulate some of these models.

The files in the folder Simulations are R files to reproduce the simulations described in the Auger-Methe et al. (2015) manuscript. Note that since the manuscript publication, I've added additional functions. If you want to reproduce the method in the Auger-Methe et al. (2015) manuscript see version 1: https://github.com/MarieAugerMethe/CCRWvsLW/tree/v1.0.

The files in the folder EmpiricalExamples demonstrate how to reproduce the analysis presented in Auger-Methe et al. (2016) manuscript. To reproduce exactly, see version 2: https://github.com/MarieAugerMethe/CCRWvsLW/tree/v2.0.

To be able to use the package you need to build it from this source code. If unfamiliar with compiling code, I recommend using RStudio (http://www.rstudio.com/). See the information found on RStudio support website to get help on how to build a package with Rstudio: https://support.rstudio.com/hc/en-us/articles/200486488-Developing-Packages-with-RStudio. Depending on your operating system you will need compilers, for more information see: https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites

You can also use devtools to install the package:
> install_github(‘MarieAugerMethe/CCRWvsLW’)
