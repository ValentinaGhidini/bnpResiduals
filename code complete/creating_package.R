# Create and install package
library(devtools)
library(roxygen2)
setwd("C:\\Users\\valen\\Desktop") # directory with bnpResidual
create("bnpResiduals")

setwd("./bnpResiduals")

document()

setwd("..")

install("bnpResiduals")


