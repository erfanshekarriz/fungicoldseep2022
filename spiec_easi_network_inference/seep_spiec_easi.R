library(SpiecEasi)
library(phyloseq)
library(igraph)
library(tidyr)
library(rstudioapi)
library(fantaxtic)
library(tidyverse)



setwd(dirname(getActiveDocumentContext()$path))       # Set working directory to source file location
getwd()                                               # Check updated working directory



### SPIEC-EASI

### Bacteria + Archaea + Fungi ###
se.bacteria.archaea.fungi <- spiec.easi(list(physeq16S.seep, physeqITS.seep),
                                 method='slr',
                                 nlambda=40,
                                 lambda.min.ratio=0.005,
                                 pulsar.params = list(thresh = 0.05, rep.num=20, ncores=4))
getStability(se.bacteria.archaea.fungi) # must be at 0.05, to bump up, increase nlambda which will give you a denser network


