################################################################################
#Cargar librerias
###############################################################################

#CTS
library(SymTS)
#CTS
library(MixedTS)



################################################################################
#Generar Datos
################################################################################


dat=rCTS(1000, 1.72, c = 0.099, ell = 0.084, mu = 0)
hist(dat,probability = T)
curve(dCTS(x, 1.72, c = 0.099, ell = 0.084, mu = 0), -0.5, 0.5,add=TRUE,col="red")
curve(dnorm(x, 0,sd(dat)), -0.5, 0.5,add=TRUE,col="Blue")






mean(dat)

sd(dat)





#https://towardsdatascience.com/why-tempered-stable-distribution-no-math-equations-5c91bb64e4e9

install.packages("SymTS")


#https://rdrr.io/rforge/MixedTS/


install.packages("MixedTS", repos="http://R-Forge.R-project.org")

#sudo apt install libgsl-dev

install.packages("GSL")

#https://CRAN.R-project.org/package=SymTS

#library(MixedTS)

library(SymTS)




