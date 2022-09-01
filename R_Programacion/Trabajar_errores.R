
#Descargar datos
library(quantmod)  


eq="KO"

i=1

while (i<=10) {
  
  
  class(try(getSymbols(eq, src = 'yahoo', from = f_init, to=f_final,
                 env = data_eq, auto.assign = T, periodicity = "d")))=="try-error"
  
  suppressWarnings(try(for(i in ls(data_eq)) data_eq[[i]] =
                         adjustOHLC(data_eq[[i]], use.Adjusted=T),silent = TRUE)) 


    i=i+1
}





#https://rpubs.com/Mentors_Ubiqum/Try_Function


2+"2"

class(try(2+"2")) == "try-error"


if(class(try(z < 5, silent = TRUE)) == "try-error"){
  print(try(z < 5, silent = TRUE)[1])






if(class(try(z < 5, silent = TRUE)) == "try-error"){
  print(try(z < 5, silent = TRUE)[1])
}else{
  print(paste("z < 5 is ", try(z < 5)))
}
