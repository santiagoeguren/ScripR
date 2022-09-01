#Extraer stock sp500
library(tidyquant)


stocks=tq_index("SP500")
stocks=stocks$symbol
