################################################################################
#Cargar librerias
################################################################################


library("nvmix")


################################################################################
#Crear Modelo
################################################################################


#"""
#X=mu +sqrt(W)A_dxk Z_dx1
#"""








#-------------------------------------------------------------------------------
#Crear la matriz A


d=1
#scale=diag(d)
A=diag(d)
A

#-------------------------------------------------------------------------------
#Crear Mu

loc=rep(0, d)

#-------------------------------------------------------------------------------
#Crear W

#Inverse Gamma

#Definir los grados de libertad
df=2

#qmix=function(u, df){
   
 #  1 / qgamma(1-u, shape = df/2, rate = df/2)

#}


#-------------------------------------------------------------------------------
#1Forma: Crear una muestra





#rt=rnvmix(200, qmix = qmix, loc = loc, scale = A, df = df)
#hist(rt)



#-------------------------------------------------------------------------------
#2Forma: Crear una muestra


rt=rnvmix(200, qmix = "inverse.gamma", loc = loc, scale = A, df = df)
hist(rt)


#-------------------------------------------------------------------------------
#Estimation

#fit_t=fitnvmix(rt, qmix = qmix, mix.param.bounds = c(0.5, 10))
#summary(fit_t)

fit_t=fitnvmix(rt, qmix =  "inverse.gamma", mix.param.bounds = c(0.5, 10))
summary(fit_t)



#-------------------------------------------------------------------------------
#Calcular probabilidades


x=1:d

#funcion acumulada

pt=pnvmix(x, qmix = "inverse.gamma", loc = loc, scale = A, df = df)
pt

#funcion densidad

dt=dnvmix(x, qmix = "inverse.gamma", loc = loc, scale = A, df = df)
dt


#-------------------------------------------------------------------------------
#Sampling normal variance mixtures using rnvmix()





#Rate es lamnda
rate=1
n=500



r.exp.1=rnvmix(n, rmix = list("exp", rate = rate))
r.exp.2=rnvmix(n, rmix = function(n) rexp(n, rate = rate))

#Graficar
#stopifnot(all.equal(r.exp.1, r.exp.2))

plot(r.exp.1)
plot(r.exp.2)


#-------------------------------------------------------------------------------
#4.3. Estimating the (log-)density with dnvmix()


set.seed(271)
d=20
df=3.9
n=2000
x=rnvmix(n, qmix = "inverse.gamma", df = df/3, scale = diag(d))
dt.1=dnvmix(x, qmix = "inverse.gamma", df = df, log = TRUE)
dt.2=dnvmix(x, qmix = function(u, df) 1 / qgamma(1 - u, shape = df/2,
                                                rate = df/2), df = df, log = TRUE)

dt.3=dnvmix(x, qmix = function(u, df) 1 / qgamma(1 - u, shape = df/2,
         rate = df/2), df = df, control = list(dnvmix.doAdapt = FALSE),
         log = TRUE)   


#-------------------------------------------------------------------------------
#4.4. Fitting normal variance mixtures with fitnvmix()



set.seed(42)
d=4
n=100
nu.=1.5
scale=cov2cor(tcrossprod(matrix(runif(d * d), ncol = d)))
scale
x=rnvmix(n, qmix = "pareto", alpha = nu., scale = scale)
m.p.b=c(0.1, 50)
qmix.=function(u, nu) (1 - u)^(-1/nu)
system.time(fit.par1<- fitnvmix(x, qmix = "pareto", mix.param.bounds = m.p.b))
system.time(fit.par2<-fitnvmix(x, qmix = qmix., mix.param.bounds = m.p.b))    

fit.par1                                 
                                                                  


set.seed(1)
qq.par=qqplot_maha(x, qmix = "pareto", alpha = fit.par1$nu,
                   loc = fit.par1$loc, scale = fit.par1$scale, plot = FALSE)
plot(qq.par)
plot(qq.par, plot.pars = list(log = "xy"))                                     

qq.par

#-------------------------------------------------------------------------------
#5. Example application

library(qrmtools)
library(rugarch)


set.seed(123)
data("SP500_const", package = "qrmdata")




time=c("2010-01-01", "2012-12-31")

x=SP500_const[paste0(time, collapse = "/"),   SP500_const_info$Subsector == "REITs"]
 
X = -returns(x)
 
uspec=rep(list(ugarchspec(distribution.model = "std")), ncol(X))

fit.ARMA.GARCH=fit_ARMA_GARCH(X, ugarchspec.list = uspec, verbose = FALSE)
                      
fits=fit.ARMA.GARCH$fit                                    

resi=lapply(fits, residuals, standardize = TRUE)

X=as.matrix(do.call(merge, resi))

colnames(X)=colnames(x)
n = nrow(X)


 qmix_ =list(constant = "constant", inverse.gamma = "inverse.gamma",
                   inverse.burr = function(u, nu) (u^(-1/nu[2]) - 1)^(-1/nu[1]),
                   pareto = "pareto")
 
 
m.p.b_ = list(constant = c(0, 1e8), inverse.gamma = c(1, 8),
                    inverse.burr = matrix(c(0.1, 0.1, 8, 8), ncol = 2), pareto = c(1, 8))

#Ejecutar este comando se demora

fit.results = lapply(1:4, function(i) fitnvmix(X, qmix = qmix_[[i]], mix.param.bounds = m.p.b_[[i]]))
  
    

l.dens <- matrix(NA_real_, ncol = 4, nrow = n)
mahas <- matrix(NA_real_, ncol = 4, nrow = n)
 for (i in 1:4) {
  
    mahas[, i] <- sqrt(mahalanobis(X, center = fit.results[[i]]$loc,
                                   
                                     cov = fit.results[[i]]$scale))
    
      order.maha <- order(mahas[, i])
      
        mahas[, i] <- mahas[order.maha, i]
        
          l.dens[, i] <- dnvmix(X[order.maha, ], qmix = qmix_[[i]],
                                
                                  loc = fit.results[[i]]$loc, scale = fit.results[[i]]$scale,
                                
                                  nu = fit.results[[i]]$nu, log = TRUE)
 }




qq.results=lapply(1:4, function(i)
    qqplot_maha(fitnvmix_object = fit.results[[i]]))


n. <- 50
u. <- seq(0.95, to = 0.999, length.out = n.)
u.matrix <- matrix(u., nrow = n., ncol = ncol(X))
tailprobs <- sapply(1:4, function(i) pnvmixcopula(1 - u.matrix,
                                                     
                                                       qmix = qmix_[[i]], scale = cov2cor(fit.results[[i]]$scale),
                                                     
                                                       nu = fit.results[[i]]$nu, control = list(pnvmix.abstol = 1e-5)))




