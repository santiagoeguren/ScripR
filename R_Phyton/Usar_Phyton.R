########################################################################
#Cargar librerias
#########################################################################
#Conecta R y Phton
library(reticulate)



################################################################################
#Establecer direccion
###############################################################################

setwd("/home/santiago/ScripR/R_Phyton")
getwd()


################################################################################
#Hay dos formas de usa Phyton in R
#1 Formas: Usas las librerias de Phyton
#2 Forma: Exportas de un archivo de Phyton

#Tene en cuenta que la base es R
################################################################################


###############################################################################
#1 Formas: Usas las librerias de Phyton
################################################################################

#-------------------------------------------------------------------------------
#Crear un vector en Phyton, mendiante un vector en R
#-----------------------------------------------------------------------------

#Pandas
py_install(packages = c("pandas", "scikit-learn"))

py_install(packages = c("pyfolio")

#Crer el objeto de R, que usa las funciones en R
pd=import("pandas")


pd$array(c(1, 2, 3))


################################################################################
#2 Forma: Exportas de un archivo de Phyton
################################################################################


#-------------------------------------------------------------------------------
#Usar una funcion de phyton
#-------------------------------------------------------------------------------

source_python("Phyton.py")



linreg_python$fit(X = mtcars[,-1], y = mtcars$mpg)

data.frame(var = c("Intercept", names(mtcars)[-1]), 
           python_coef = c(linreg_python$intercept_, linreg_python$coef_))



#-------------------------------------------------------------------------------
#Usar un objeto de Phyton
#-------------------------------------------------------------------------------

source_python("Phyton.py")

df$`0`
df$`1`
df$`2`
