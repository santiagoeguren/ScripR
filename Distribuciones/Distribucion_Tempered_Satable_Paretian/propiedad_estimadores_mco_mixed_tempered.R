################################################################################
#Cargar librerias
################################################################################


library(MixedTS)




#################################################################################
#Armar modelo lineal
################################################################################



ParamEx1=setMixedTS.param(mu0= 0 , mu=  0, sigma=   0.71361240  , a=  4.89189348  ,
                               alpha=    1.52875518, lambda_p=  0.15281382 , lambda_m=0.09713054, 
                                                         Mixing="Gamma")
t_values=NULL
y=NULL

ii=1


while (ii<=1000) {
  


Rand1=rMixedTS(x=201,object=ParamEx1, setSup=10,setInf=-10,N=2^9)


x=c(0:200)

i=1


while (i<=201) {
  
  
  
  y[i]=10+0*x[i]+Rand1@Data[i]
  
  i=i+1
  
}


plot(x,y)


model=lm(y~x)


summary(model)


t_values[ii]=summary(model)[["coefficients"]][, "t value"][2]

ii=ii+1

}



hist(t_values)


shapiro.test(t_values)
