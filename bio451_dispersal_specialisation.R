
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
url1 <- "https://docs.google.com/spreadsheets/d/1NyupV8XW7bWqi_Pxw476bbNF-sovFM8JVXmH2PgdqRw/edit#gid=0"
fish_dat <- read_csv(construct_download_url(url1))

# load the

# convert the decimals into points
fish_dat %>%
  mutate_all(all_vars(), gsub(",", ".", .))

fish_dat[1:10, 5:10] %>%
  mutate_all(~gsub(",", ".", .))

# Get variable names that are numbers
var_nums <- lapply(fish_dat, 
       function(r) {if_else(grepl(pattern = paste(c(0:9), collapse = "|"), x = r) == TRUE, 1 , 0) %>% sum()} ) %>%
  Filter(function(x) sum(x) > 0, .) %>%
  names()
var_nums

# Get variables where the commas need to be replaced
lapply(fish_dat, 
       function(r) {if_else(grepl(pattern = ",", x = r) == TRUE, 1 , 0) %>% sum()} ) %>%
  Filter(function(x) sum(x) > 0, .) %>%
  names()


# replace all the commas with points

match()



