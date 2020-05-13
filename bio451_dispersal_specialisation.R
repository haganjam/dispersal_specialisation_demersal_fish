
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
url1 <- "https://docs.google.com/spreadsheets/d/1NyupV8XW7bWqi_Pxw476bbNF-sovFM8JVXmH2PgdqRw/edit#gid=489914387"
fish_dat_raw <- read_csv(construct_download_url(url1))

# filter the Nanninga and Manica data for now
fish_dat <- fish_dat_raw %>%
  filter(Reference == "nanninga_manica_2018")

# remove some useless columns
names(fish_dat)

fish_dat <- select(fish_dat, -"Comment 1", -"Comment 2", -"binomial_2")

# get variable names that are numbers
var_nums <- lapply(fish_dat, 
       function(r) {if_else(grepl(pattern = paste(c(0:9), collapse = "|"), x = r) == TRUE, 1 , 0) %>% sum()} ) %>%
  Filter(function(x) sum(x) > 0, .) %>%
  names()
var_nums

# remove reference from this because we don't want this as a numeric variable
var_nums <- var_nums[var_nums != "Reference"]

# get variables where the commas need to be replaced
var_coms <- lapply(fish_dat, 
       function(r) {if_else(grepl(pattern = ",", x = r) == TRUE, 1 , 0) %>% sum()} ) %>%
  Filter(function(x) sum(x) > 0, .) %>%
  names()

# get unique variable numbers without commas
var_nums <- var_nums[-match(var_coms, var_nums)]

# check these variables
select(fish_dat, var_coms) %>%
  lapply(function(x) { if_else(unique(x) == "NA", 1, 0) })

var_coms

fish_dat$min_ucrit %>% as.numeric()

gsub(",", ".", fish_dat$latitude) %>% as.numeric()

bind_cols(x1 = fish_dat$latitude, x2 = gsub(",", ".", fish_dat$latitude) %>% as.numeric()) %>%
  View()

# convert the decimals into points and then numbers
fish_dat <- fish_dat %>%
  mutate_at(vars(var_coms), ~gsub(",", ".", .) %>% as.numeric())

# convert the numbers without decimals to numeric
fish_dat <- fish_dat %>%
  mutate_at(vars(var_nums), ~as.numeric(.))

c(1, 2, NA) %>% as.numeric()

# replace all the commas with points

match()



