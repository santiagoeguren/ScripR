
install.packages("nvmix")


library(nvmix)
library(RColorBrewer)
library(lattice)


#https://cran.r-project.org/web/packages/nvmix/vignettes/nvmix_functionality.html

alfa=0
beta=0
sdt=1

y=NULL


i=1

while (i<10000) {
  


v=rgamma(1,1,1)

y[i]=alfa+beta*v+sdt*sqrt(v)*rnorm(1,0,1)
  
help("rgamma")

i=i+1

}


hist(y,breaks = 60)





