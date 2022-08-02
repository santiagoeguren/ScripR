
#-------------------------------------------------------------------------------
#Librerias
#-------------------------------------------------------------------------------


library(purrr)  #Armar la matrix
library(Quandl) #Para descarcar los datos de la pagina
Quandl.api_key("d9EidiiDWoFESfdk5nPy") #Para descargar los datos de mi usuario
library(astsa) #para usar acf2
library(urca) #test de dickey fuller
library(tseries)#test de dickey fuller
library(vars) #calcular vars



#---------------------------------------------------------------------------



E=rnorm(240)
X=cumsum(E)
plot(X,type="l")

#----------------------------------------------------------------------------
#Augmented Dickey-Fuller Test
#----------------------------------------------------------------------------


#ver: https://freakonometrics.hypotheses.org/12729


#1 librerya urca

#H0: no estacionario Vs  H1: estacionario

#1.1 Dickey-Fuller - Standart

#yt=ro. yt-1 + et
#que se transforma
#deltayt=(1-ro). yt-1 + et


df=ur.df(X,type="none",lags=0)
df
summary(df)



#1.2 Dickey-Fuller - con un "y" lag

#deltayt=(1-ro). yt-1 +B.deltayt-1+ et




lags=1
z=diff(X)
n=length(z)
z.diff=embed(z, lags+1)[,1]
z.lag.1=X[(lags+1):n]
k=lags+1
z.diff.lag = embed(z, lags+1)[, 2:k]
plot(z.diff.lag,type="l")

summary(lm(z.diff~0+z.lag.1+z.diff.lag ))




df=ur.df(X,type="none",lags=1)
df
summary(df)




#1.3 Dickey-Fuller - drift

#deltayt=alpa+(1-ro). yt-1 +B.deltayt-1+ et




summary(lm(z.diff~1+z.lag.1+z.diff.lag ))


df=ur.df(X,type="drift",lags=1)
df
summary(df)



#1.4 Dickey-Fuller - trend

#deltayt=alpa+beta*t+(1-ro). yt-1 +B.deltayt-1+ et



temps=(lags+1):n
summary(lm(z.diff~1+temps+z.lag.1+z.diff.lag ))


df=ur.df(X,type="trend",lags=1)
summary(df)

#--------------------------------------------------------------------------------
#test  Philips-Perron
#-----------------------------------------------------------------------------

#H0: no estacionario Vs  H1: estacionario


#https://www.rdocumentation.org/packages/aTSA/versions/3.1.2/topics/pp.test


#pp.test(X, type = "Z(alpha)", lag.short = TRUE, output = TRUE)

#lag.short = TRUE, we use the default number of Newey-West lags
#"Z(alpha)" con drift
#"Z(t_alpha)" con tendencia con drift
#

#pp.test(X,type = "Z(alpha)")


#library urca

#ur.pp(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"), lags = c("short", "long"), use.lag = NULL)
# trunc(4*(n/100)^0.25), otherwise trunc(12*(n/100)^0.25) --->short an long

df=ur.pp(E, type ="Z-tau",model = c("trend"),lags = "short")
df
df@cval
df@teststat
df@lag

summary(df)
plot(df)



library(urca) #test de dickey fuller




#--------------------------------------------------------------------------------
#test  elliott rothenberg stock
#-----------------------------------------------------------------------------

#ur.ers(y, type = c("DF-GLS", "P-test"), model = c("constant", "trend"), lag.max = 4)


df=ur.ers(X, type =  "DF-GLS", model =  "trend", lag.max = 4)
df

summary(df)



#------------------------------------------------------------------------
#Zivot \& Andrews Unit Root Test
#------------------------------------------------------------------------


#ur.za(X, model = c("intercept", "trend", "both"), lag=NULL)

df=ur.za(X, model ="both", lag=NULL)
summary(df)
plot(df)







