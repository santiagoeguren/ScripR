###############################################################################
#Library
################################################################################

library(stabledist)



################################################################################
#Generar muestra
###############################################################################

#------------------------------------------------------------------------------
#Cargar parametros

n=1000
alpha=1.72
beta=0
gamma=1
delta=0
pm=2

r=rstable(n = n, alpha = alpha, beta = beta, gamma =  gamma, delta = delta ,pm=pm)

#--------------------------------------------------------------------------------
#Graficar y calculos basicos
#------------------------------------------------------------------------------

plot(r, type = "l", main = "stable",
     col = "steelblue")
grid()

#-------------------------------------------------------------------------------
#Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white",
     col = "steelblue")

x=seq(-5, 5, 0.25)
lines(x, dstable(x, alpha = alpha, beta =  beta,  gamma =  gamma, delta = delta ,pm=pm,
                 tol= 1e-3), lwd = 2)

## Plot df and compare with true df:
plot(ecdf(r), do.points=TRUE, col = "steelblue",
     main = "Probabilities: ecdf(rstable(1000,..)) and true cdf F()")
rug(r)
lines(x, pstable(q = x,alpha = alpha, beta =  beta,  gamma =  gamma, delta = delta ,pm=pm),
      col="#0000FF88", lwd= 2.5)


## Compute quantiles:
qstable(0.025, alpha = alpha, beta =  beta,  gamma =  gamma, delta = delta ,pm=pm)
quantile(r,0.025)


#----------------------------------------------------------------------------
#Comparar dos distribución
#------------------------------------------------------------------------------

## Switching sign(beta) <==> Mirror the distribution around x == delta:
curve(dstable(x, alpha=1.2, beta = .8, gamma = 3, delta = 2), -10, 10)
curve(dstable(x, alpha=1.2, beta = -.8, gamma = 3, delta = 2),
      add=TRUE, col=2)


## or the same
curve(dstable(2*2-x, alpha=1.2, beta = +.8, gamma = 3, delta = 2),
      add=TRUE, col=adjustcolor("gray",0.2), lwd=5)
abline(v = 2, col = "gray", lty=2, lwd=2)
axis(1, at = 2, label = expression(delta == 2))


