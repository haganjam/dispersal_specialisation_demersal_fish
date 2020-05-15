
# Title: dispersal-specialisation trade-offs in marine demersal fish

# load relevant libraries

library(tidyverse)
library(here)
library(broom)
library(corrplot)
library(gsheet)
library(viridis)
library(MuMIn)
library(vegan)

# make a folder to export figures or tables
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# load the raw main database
fish_dat <- read_csv( here("data/1_alldata.csv") )
str(fish_dat)

# load the diet taxon matrix
diet_mat <- read_csv( here("data/diet taxon.csv") )
str(diet_mat)

# load the diet data for each species
diet_spp <- read_csv( here("data/2_diet_specialisation.csv") )
str(diet_spp)

# there are two incorrect column names
diet_spp <- 
  diet_spp %>%
  rename(plants_other_plants_terrestrial_plants = `plants_other_plants_terrestrial plants`,
         zooplankton_fish_early_stages_fish_eggs_larvae = `zooplankton_fish(early_stages)_fish_eggs_larvae`)

# make sure the diet names are harmonised
# if names are different give it a zero and then subset the zeros
bind_cols(mat = diet_mat$unique_id,
          spp = diet_spp %>% 
            select(-species) %>% 
            names()
          ) %>%
  mutate(harm = if_else(mat == spp, 1, 0)) %>%
  filter(harm == 0) %>%
  View()

# they are not harmonised

# harmonise the names
diet_mat$unique_id <- 
  diet_spp %>%
  select(-species) %>% 
  names()
  
# check that the names were harmonised properly
bind_cols(mat = diet_mat$unique_id,
          spp = diet_spp %>% 
            select(-species) %>% 
            names()
          ) %>%
  mutate(harm = if_else(mat == spp, 1, 0)) %>%
  filter(harm == 0)

# now the names are harmonised


# load the habitat taxon matrix
hab_mat <- read_csv( here("data/habitat taxon matrice.csv") )

# load the habitat data for each species
hab_spp <- read_csv( here("data/3_habitat_specialisation.csv") )

# rename X1 column as species
hab_spp <- 
  hab_spp %>%
  rename(species = "X1")


### justifying u-crit and pld values for as dispersal traits

# for this, we only use the Nanninga and Manica data
ibd_dat <- 
  fish_dat %>% 
  filter(reference == "nanninga_manica_2018")

# filter out the rows without ibd values
ibd_dat <- filter(ibd_dat, !is.na(ibd))



# filter out rows with pld greater than zero
# this excludes the only direct developing species in the dataset
ibd_dat <- filter(ibd_dat, mean_pld > 0)

# select certain variables that might predict the ibd slope
ibd_vars <- c("mean_ucrit", "ibd", "mean_pld", "max_ucrit", "mean_larval_size")

# plot histograms of these variables
ibd_dat %>%
  gather(ibd_vars, key = "var", value = "val") %>%
  ggplot(data = .,
       mapping = aes(x = val)) +
  geom_histogram() +
  facet_wrap(~var, scales = "free") +
  theme_classic()

# transform and plot again
mutate_at(ibd_dat, vars(ibd_vars), ~log(.)) %>%
  gather(ibd_vars, key = "var", value = "val") %>%
  ggplot(data = .,
         mapping = aes(x = val)) +
  geom_histogram() +
  facet_wrap(~var, scales = "free") +
  theme_classic()


# log10 transform variables that could predict the ibd slope to improve their distributions
ibd_dat <- mutate_at(ibd_dat, vars(ibd_vars), ~log(.))

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


### calculate dietary diversity for each species

# for this, we need the taxon matrix and species matrix for the diets
diet_mat

diet_spp

# need to convert diet_mat into a regular dataframe (due to the vegan package function)
diet_mat <- as.data.frame(diet_mat)

# set row names for the diet_mat
row.names(diet_mat) <- diet_mat$unique_id

# remove the unique_id column
diet_mat <- select(diet_mat, -unique_id)

# set the food_i, food_ii and food_iii variables names as Family, Genus, Species
names(diet_mat) <- c("Family", "Genus", "Species")

# use the taxa2dist function to create a dissimilarity matrix
diet_dis <- taxa2dist(diet_mat, varstep = TRUE, check = TRUE)


# use the species data and the dissimilarity matrix to calculate taxonomic distinctness

# check for NAs
lapply(diet_spp, function(x) {sum(if_else(is.na(x), 1, 0)) })

lapply(diet_spp, function(x) { 
  bind_cols(id = diet_spp$species, na_row = is.na(x) ) %>%
    filter(na_row == TRUE) %>%
    pull(id) 
  } )

# why is the "zoobenthos_sponges_tunicates_sponges" column have missing data?

# first, remove the rows without any diet data (i.e. the nas)
# second, remove the species column
# third conver to presence absence
diet_div_dat <- 
  diet_spp %>%
  filter_all( all_vars( (!is.na(.)) ) ) %>%
  select(-species) %>%
  as.data.frame()

diet_div <- 
  taxondive(diet_div_dat, diet_dis)

# attach these diet specialisation indices to a dataframe with the species names
spp_diet_div <- 
  diet_spp %>%
  filter_all( all_vars( (!is.na(.)) ) ) %>%
  select(species) %>%
  mutate(food_groups = diet_div$Species,
         diet_dplus = diet_div$Dplus)

pairs(select(spp_diet_div, -species))






# An example

# Make the taxon matrix using Food_I, Food_II and Food_II
df <- data.frame(Family = c("Detritus", "Detritus", "plants"),
                 Genus = c("detritus", "detritus", "phytoplankton"),
                 Species = c("debris", "carcasses", "blue_green_algae"))
df

# Give it rown names as a unique name (binomial equivalent)
row.names(df) <- c("Detritus_detritus_debris", "Detritus_detritus_carcasses", "plants_phytoplankton_blue_green_algae")

df

# Make a dataframe for three different species (rows) and whether they eat different foods
df_pres <- data.frame(Detritus_detritus_debris = c(1, 1, 0),
                      Detritus_detritus_carcasses = c(1, 1, 1),
                      plants_phytoplankton_blue_green_algae = c(0, 1, 1))
df_pres

# Calculate a dissimilarity matrix from the taxon matrix
df_dis <- taxa2dist(df, varstep = TRUE)
df_dis

# Use the species data and the dissimilarity matrix to calculate taxonomic distinctness
taxondive(df_pres, df_dis)







