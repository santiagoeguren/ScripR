########################################################################
#Cargar librerias
#########################################################################
#Conecta R y Phton
library(reticulate)


#"scikit-plot
py_install(packages = c("scikit-plot"))

#pycaret[full]
py_install(packages = c("pycaret[full]"))

#Pandas y scikit-learn
py_install(packages = c("pandas", "scikit-learn"))

