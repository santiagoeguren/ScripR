

stocks=stocks[!stocks %in% "LHX"]
stocks=stocks[!stocks %in% "FE"]
stocks=stocks[!stocks %in% "DG"]

#Saber si existe un stocks en la lista

stocks[stocks=="KO"]
