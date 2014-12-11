library(CCRWvsLW)

mov <- simmCCRW(500,0.9,0.9,0.01,0.001,10,1,0.5)
plot(mov)

movLikelihoods(mov)
