
# Title: dispersal-specialisation trade-offs in marine demersal fish

# load relevant libraries

library(tidyverse)
library(here)
library(broom)
library(corrplot)
library(gsheet)
library(viridis)

# make a folder to export figures or tables
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# load the data from google sheets






