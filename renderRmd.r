# script to run all Rdm files to generate figure for the paper by 
# Danilova et al. "Sequencing analysis reveals evidence of immune activation in advanced HER2 negative breast cancer responders treated with entinostat + nivolumab + ipilimumab"

library(knitr)
# list all rmd files
files = list.files(path = ".", pattern = "*.Rmd$", full.names = T)
# render all rmds
for(i in files) rmarkdown::render(i, output_dir = ".")
