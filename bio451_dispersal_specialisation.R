
# Title: dispersal-specialisation trade-offs in marine demersal fish

# load relevant libraries

library(tidyverse)
library(here)
library(broom)
library(corrplot)
library(gsheet)
library(viridis)
library(MuMIn)

# make a folder to export figures or tables
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# load the data from google sheets
# url1 <- "https://docs.google.com/spreadsheets/d/1NyupV8XW7bWqi_Pxw476bbNF-sovFM8JVXmH2PgdqRw/edit#gid=489914387"
# fish_dat_raw <- read_csv(construct_download_url(url1))

# for now, just read in the data to examine the predictors of ibd slope

ibd_dat <- read_csv(here("data/prel_dat.csv"))

# subset out the Nanninga and Manica (2018) data
#ibd_dat <- fish_dat_raw %>%
  #filter(Reference == "nanninga_manica_2018")

# filter out the rows without ibd values
ibd_dat <- filter(ibd_dat, !is.na(ibd))

# filter out rows with pld greater than zero
# this excludes the only direct developing species
ibd_dat <- filter(ibd_dat, mean_pld > 0)

# select certain variables that might predict the ibd slope
ibd_vars <- c("mean_ucrit", "ibd", "mean_pld", "max_ucrit", "mean_larval_size")

# log10 transform variables that could predict the ibd slope
ibd_dat <- mutate_at(ibd_dat, vars(ibd_vars),
                     ~log(.))

# count total number of data points
nrow(ibd_dat)

# check the correlation among these variables 
ibd_dat %>% 
  select(ibd_vars) %>%
  cor() %>%
  corrplot(method = "number")

ibd_dat %>% 
  select(ibd_vars) %>%
  pairs()

# mean_ucrit and max_ucrit are extremely correlated so we can ignore max_ucrit
ibd_vars <- ibd_vars[ibd_vars != c("max_ucrit")]

# check the distributions of these variables
ggplot(data = ibd_dat %>% 
         gather(ibd_vars, key = "group", value = "val"),
       mapping = aes(x = val, fill = group)) +
  geom_histogram() +
  facet_wrap(~ group, scales = "free") +
  theme_classic()

# how many predictor variables are appropriate for this dataset
nrow(ibd_dat)

# maximum four predictor variables given the 23 data points

# does ibd differ between egg types
ggplot(data = ibd_dat,
       mapping = aes(x = mean_ucrit, y = ibd, colour = egg_type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

ggplot(data = ibd_dat,
       mapping = aes(x = mean_pld, y = ibd, colour = egg_type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

# fit all possible models to compare their r2 values
lm_glob <- lm(ibd ~ mean_ucrit + mean_pld + egg_type + mean_ucrit:mean_pld +
                mean_ucrit:egg_type + mean_pld:egg_type, 
              data = ibd_dat,
              na.action = "na.fail")

all_mods <- dredge(lm_glob, extra = "r.squaredLR", fixed = "mean_ucrit")

# view these models
all_mods %>% 
  as_tibble() %>% 
  arrange(desc(r.squaredLR))  %>%
  View()

# based on this model selection, mean_pld, mean_ucrit and egg_type are all important variables

# use a principle components analysis to make a single dispersal variable
ibd_dat %>% names()

pca_glob <- princomp(~ mean_ucrit + mean_pld, data = ibd_dat,
               cor = TRUE)

pca_glob %>% 
  summary()

# extract pca component 1 scores from the global pca object
pca_scores <- pca_glob$scores %>% 
  as_tibble() %>%
  pull(Comp.1)

# add this as a variable to ibd_dat data
ibd_dat <- ibd_dat %>%
  mutate(pca_comp.1 = pca_scores)


# use a loop to fit different models to ibd data

# set up the different models explanatory variables
exp_vars <- list(c("mean_ucrit"),
                 c("mean_pld"), 
                 c("pca_comp.1"),
                 c("pca_comp.1", "egg_type"))

# set the response variable
resp_var <- c("ibd")

# set the dataset
lm_dat <- ibd_dat

# set up output lists
cof_out <- vector("list", length = length(exp_vars) )
diag_out <- vector("list", length = length(exp_vars) )

# run the loop to fit the models
for (i in seq_along(1:length(exp_vars)) ) {
  
  lm_out <- lm( reformulate(exp_vars[[i]], resp_var), data = lm_dat)
  
  cof_out[[i]] <- tidy(lm_out) %>%
    mutate(vars = paste(exp_vars[[i]], collapse = "_") )
  
  diag_out[[i]] <- glance(lm_out) %>% 
    mutate(aicc = AICc(lm_out),
           vars = paste(exp_vars[[i]], collapse = "_"))

  }

# check the diagnostic models
diag_out %>% 
  bind_rows(.id = "model") %>%
  arrange(aicc) %>%
  View()

# check the coefficients of the best model
cof_out

# based on this, we accept that pca_comp.1 is the best predictor of the ibd slope

# how is the pca_com.1 axis related to mean_ucrit and mean_pld?
# what does this axis mean?

ggplot(data = ibd_dat %>%
         gather(mean_ucrit, mean_pld, key = "larval_trait", value = "value"),
       mapping = aes(x = pca_comp.1, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~larval_trait) +
  theme_classic()




### getting rid of commas etc.

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



