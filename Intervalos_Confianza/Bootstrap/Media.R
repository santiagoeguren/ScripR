################################################################################
#Librerias
################################################################################

library(mosaic)



#-------------------------------------------------------------------------------
#Distribucion Normal

Distancia=rnorm(500,18,13)



media=mean(Distancia)
desvio=sd(Distancia)

#-------------------------------------------------------------------------------
#Bootstraping

media_boot=NULL
for(i in 1:5000){

media_boot[i]=mean(sample(Distancia,500,replace = T))

}


#-------------------------------------------------------------------------------
#Distribucion bootstrap
hist(media_boot)

shapiro.test(media_boot)

#"""
#La media va estar centrada en la media original
#"""

#-------------------------------------------------------------------------------
#Desviavion estándar bootstrap

#"""
#Al desvio estandar se llama error estandar sigma/sqrt(n)
#"""


std_boot=sd(media_boot)


#-------------------------------------------------------------------------------
#intervalo de confianza

mean(media_boot)-sd(media_boot)*2
mean(media_boot)+sd(media_boot)*2

#Esta es la forma

media-sd(media_boot)*2
media+sd(media_boot)*2

#Otra forma


quantile(media_boot,0.025)
quantile(media_boot,0.975)
