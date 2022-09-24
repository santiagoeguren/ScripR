################################################################################
#Cargar librerias
################################################################################


library(MixedTS)


################################################################################
#Generar Datos
################################################################################


#-------------------------------------------------------------------------------
# Density of MixedTS with Gamma
#-------------------------------------------------------------------------------

ParamEx1=setMixedTS.param(mu0=0, mu=0, sigma=0.5, a=2,
                           alpha=1.72, lambda_p=1, lambda_m=1, 
                           Mixing="Gamma")







#-------------------------------------------------------------------------------
#Generar valores aleatorios

dat=rMixedTS(1000,object=ParamEx1)

dat@Data


#-------------------------------------------------------------------------------
#Determinar funcion densidad

x=seq(-2,2,length=1000)

dens=dMixedTS(x=x,object=ParamEx1,setSup=0.1,setInf=-0.1,N=2^7)



#-------------------------------------------------------------------------------
#Graficar

plot(dens)
hist(dat@Data, n = 30, probability = TRUE,border = "white",
     col = "lightgray",add=T)





################################################################################
#Estimacion
################################################################################





ParamEx1=setMixedTS.param(mu0=0.1, mu=1, sigma=0.1, a=0.1,
                           alpha=0.1, lambda_p=0.1, lambda_m=0.1, Mixing="Gamma")




Rand1=rMixedTS(x=5000,object=ParamEx1, setSup=10,setInf=-10,N=2^9)


Rand1@Data=rnorm(100,100,2)




est1=mle.MixedTS(object=Rand1 , setSup=10,setInf=-10,N=2^9)


summary(est1)







































#-------------------------------------------------------------------------------
# Density of MixedTS with IG
#-------------------------------------------------------------------------------





Mix<-"User"

logmgf<-("lamb/mu1*(1-sqrt(1-2*mu1^2/lamb*u))")

parMix<-list(lamb=1,mu1=1)

ParamEx2<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=logmgf,
                           alpha=0.8, lambda_p=4, lambda_m=1,
                           Mixing=Mix,paramMixing=parMix)

x<-seq(-3,1,length=100)

dens2<-dMixedTS(x=x,object=ParamEx2,setSup=10,setInf=-10,N=2^7)

plot(dens2)


















# First Example:
# We define the Mixed Tempered Stable using the function setMixedTS.param


ParamEx1<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=1.5,
                           alpha=0.8, lambda_p=4, lambda_m=1, Mixing="Gamma")

# We generate a sample using the rMixedTS method
set.seed(100)
Rand1 <- rMixedTS(x=5000,object=ParamEx1, setSup=10,setInf=-10,N=2^9)

# Estimate procedure
## Not run: 
est1<-mle.MixedTS(object=Rand1 , setSup=10,setInf=-10,N=2^9)
# Show results

summary(est1)






################################################################################
#Referencias
################################################################################


#https://www.rdocumentation.org/packages/MixedTS/versions/1.0.4/topics/dMixedTS-methods
#https://rdrr.io/rforge/MixedTS/man/mle.MixedTS.html
