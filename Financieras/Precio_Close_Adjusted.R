########################################################################
#Cargar librerias
#########################################################################

#Descargar datos
library(quantmod)  
#Graficos
library(ggplot2) 
#Test de Dickey-Fuller
library(tseries)
#Zivot & Andrews Unit Root Test
library(urca) 
#Estimar correlograma  
library(astsa)
#Extraer stock sp500
library(rvest)
#Ordenar datos funcios selec
library(dplyr)
#Funcion tiempo
library(lubridate)
#Ordenar datos para graficos
library(tidyr)
#Garch
require(fGarch)

###############################################################################
#Comparar 
################################################################################





getSymbols("CTAS",
           from = '2012-10-01',
           to = '2022-11-22',
           periodicity = "daily")






precios_cierre=as.numeric(CTAS$CTAS.Close)
precios_adjusted=as.numeric(CTAS$CTAS.Adjusted)



diferencia=precios_cierre-precios_adjusted

plot(diferencia)

tail(diferencia)

#----------------------------------------------------------------------------------
#Transformar en rentabilidad inicio 1


valor_cierre=precios_cierre/precios_cierre[1]
valor_adjusted=precios_adjusted/precios_adjusted[1]

#valor_cierre=precios_cierre
#valor_adjusted=precios_adjusted

#Armar data frame  Close
data_rent_cierre=data.frame(fecha=index(CTAS), valor_cierre=valor_cierre)


#Armar data frame Adjusted
data_rent_Adjusted=data.frame(fecha=index(CTAS),valor_adjusted=valor_adjusted)





#----------------------------------------------------------------------------------
#Graficar

data_valores=cbind(data_rent_cierre,data_rent_Adjusted$valor_adjusted)

#Ordenar datos
colnames(data_valores)=c("fecha","valor_cierre","valor_adjusted")

data_valores_p=data_valores %>%
  select(fecha,valor_cierre, valor_adjusted) %>% gather(key = "variable", value = "value", -fecha)

#Grear grafico

g=ggplot(data_valores_p, aes(x = fecha, y = value))
g=g+geom_area(aes(color = variable, fill = variable), alpha = 0.05, position = position_dodge(0.8)) 
g=g+scale_color_manual(values = c("#00AFBB", "red")) 
g=g+scale_fill_manual(values = c("#00AFBB", "red"))
g=g+scale_x_date(date_labels = "20%y", date_breaks = "1 years")
g=g+xlab("Tiempo") + ylab("Rendimiento acumulado")
g=g+theme_minimal()
g






