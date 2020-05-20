
# Title: dispersal-specialisation trade-offs in marine demersal fish

# load relevant libraries

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)
library(corrplot)
library(MuMIn)
library(vegan)


# make a folder to export figures or tables
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}


### justifying u-crit and pld values for as dispersal traits

# load the raw main database
ibd_dat <- read_csv( here("data/1_alldata.csv") )
str(ibd_dat)

# for this, we only use the Nanninga and Manica data
ibd_dat <- 
  ibd_dat %>% 
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

# ibd slope is heavily skewed so we log10 transform it

# log10 transform variables that could predict the ibd slope to improve their distributions
ibd_dat <- 
  ibd_dat %>%
  mutate_at(vars(ibd_vars), ~log(.))

# check the distribution of the variables again
ibd_dat %>%
  gather(ibd_vars, key = "var", value = "val") %>%
  ggplot(data = .,
         mapping = aes(x = val)) +
  geom_histogram() +
  facet_wrap(~var, scales = "free") +
  theme_classic()

# this improves the distribution considerably

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

# how many predictor variables are appropriate for this dataset
nrow(ibd_dat)

# maximum four predictor variables given the 23 data points

# does ibd differ between egg types
ggplot(data = ibd_dat,
       mapping = aes(x = mean_ucrit, y = ibd, colour = egg_type)) +
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
  arrange(AICc)  %>%
  View()

# based on this model selection, mean_pld, mean_ucrit and egg_type are all important variables

# use a principle components analysis to make a single dispersal variable
ibd_dat %>% 
  names()

pca_glob <- princomp(~ mean_ucrit + mean_pld, data = ibd_dat,
               cor = TRUE)

pca_glob %>% 
  summary()

# extract pca component 1 scores from the global pca object
pca_scores <- 
  pca_glob$scores %>% 
  as_tibble() %>%
  pull(Comp.1)

# add this as a variable to ibd_dat data
ibd_dat <- 
  ibd_dat %>%
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
# pca is based on ln-transformed mean_ucrit and ln-transformed mean_pld

# how is the pca_com.1 axis related to mean_ucrit and mean_pld?
# what does this axis mean?

ggplot(data = ibd_dat %>%
         gather(mean_ucrit, mean_pld, key = "larval_trait", value = "value"),
       mapping = aes(x = pca_comp.1, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~larval_trait) +
  theme_classic()

# high pca values indicate high dispersal



### load the raw data and correct the problems with it

all_raw <- read_csv( here("data/1_alldata.csv") )
str(all_raw)

# remove the direct developing species
all_raw <- 
  all_raw %>% filter(mean_pld > 0 | is.na(mean_pld))

# fix the incorrect names (from Niklas)

# amphirion clarkii has a latitude range of 60 but min_lat = max_lat=30 
# this must be a mistake, so we need to change min_lat on that row to -30
all_raw$min_lat[which(all_raw$binomial == "amphirion_clarki")] <- (-30)

# change demersal egg type (chromis_chromis) to benthic 
all_raw$egg_type[which(all_raw$egg_type == "dem")] <- c("ben")

# also change "Sparidae" to "sparidae" in the family column
all_raw$family[which(all_raw$family == "Sparidae")] <- c("sparidae")

# amphirion should be amphiprionin binomial and genus names 
all_raw$genus[which(all_raw$genus == "amphirion")] <- c("amphiprion")

# binomial for certain amphiprion species are incorrect: (amphirion_clarki, amphirion_melanopus, amphirion_percula)
all_raw$binomial[which(all_raw$binomial == "amphirion_clarki")] <- c("amphiprion_clarki")
all_raw$binomial[which(all_raw$binomial == "amphirion_melanopus")] <- c("amphiprion_melanopus")
all_raw$binomial[which(all_raw$binomial == "amphirion_percula")] <- c("amphiprion_percula")

# Amphiprioninae AND Pomacentrinae should be Pomacentridae in family names
all_raw$family[which(all_raw$family == "amphiprioninae")] <- c("pomacentridae")
all_raw$family[which(all_raw$family == "pomacentrinae")] <- c("pomacentridae")

# additional mistakes with family
all_raw$family[which(all_raw$family == "pomacanthidae")] <- c("pomacentridae")
all_raw$family[which(all_raw$family == "lutjanidae")] <- c("caesionidae")
all_raw$family[which(all_raw$family == "pomacanthidae")] <- c("pomacentridae")
all_raw$family[which(all_raw$family == "monacanthidae")] <- c("nemipteridae")

# additional mistake with genus: abudefdud should be abudefduf
all_raw$genus[which(all_raw$genus == "abudefdud")] <- c("abudefduf")

# the problem here is that we now have extra data points that we need to correct

### amphiprion_melanopus
all_raw %>%
  filter(binomial == "amphiprion_melanopus") %>%
  View()

# calculate the mean_ucrit for amphiprion_melanopus
a_m <- 
  all_raw %>%
  filter(binomial == "amphiprion_melanopus") %>%
  pull(mean_ucrit) %>%
  mean()

# replace the ucrit value with this mean
all_raw$mean_ucrit[which(all_raw$binomial == "amphiprion_melanopus")] <- a_m

# amphiprion_melanopus also has two different latitude values: convert to -11
all_raw$latitude_species[which(all_raw$binomial == "amphiprion_melanopus")] <- c(-11)

### amphiprion_percula
all_raw %>%
  filter(binomial == "amphiprion_percula") %>%
  View()

# calculate the mean_ucrit for amphiprion_percula
a_p <- 
  all_raw %>%
  filter(binomial == "amphiprion_percula") %>%
  pull(mean_ucrit) %>%
  mean()

# replace the ucrit value with this mean
all_raw$mean_ucrit[which(all_raw$binomial == "amphiprion_percula")] <- a_p

### atherina_presbyter, there are several inconsistencies:
all_raw %>%
  filter(binomial == "atherina_presbyter") %>%
  View()

# mean_larval_size is missing for one row of the data and must be filled in
all_raw$mean_larval_size[which(all_raw$binomial == "atherina_presbyter")] <- 16

# food_ i should be one
all_raw$food_i[which(all_raw$binomial == "atherina_presbyter")] <- 1

# iucn_hab_sub is unclear but I will make it 0.35
all_raw$habitat_iucn_sub[which(all_raw$binomial == "atherina_presbyter")] <- 0.35


# take the distinct rows out (i.e. remove any duplicates)
all_raw <- 
  all_raw %>%
  select(-reference, -mean_ucrit_subpopulation, -min_ucrit, -max_ucrit, -latitude_subpopulation,
         -comment_1, -comment_2, -binomial_2) %>%
  distinct()

nrow(all_raw)
all_raw$binomial %>%
  unique() %>%
  length()

# now the data are correct I think... but we will probably find more mistakes


### create a dispersal trait variable from mean_ucrit and mean_pld

# create a copy of the all raw to be modified
disp_axis <- all_raw

# check for NAs
lapply(disp_axis, function(x) {  sum( if_else(is.na(x), 1, 0) )  })

# subset out only columns that are useful for this analysis
names(disp_axis)

disp_axis <- 
  disp_axis %>%
  select("family", "genus", "species", "binomial", "mean_ucrit", "mean_pld")

# filter rows without both mean_ucrit and mean_pld values
disp_axis <- 
  disp_axis %>%
  filter_at(vars(c("mean_ucrit", "mean_pld")), all_vars(!is.na(.)))

lapply(disp_axis, function(x) {range(x, na.rm = TRUE)})

# ln-transform the mean_ucrit and mean_pld data
disp_axis <- 
  disp_axis %>%
  mutate_at(vars(c("mean_ucrit", "mean_pld")), ~log(.))

# dispersal pca axis for all species
pca_disp <- princomp(~ mean_ucrit + mean_pld, data = disp_axis, cor = TRUE)

# how much variation does this axis explain?
summary(pca_disp)

comp.1_scores <- 
  pca_disp$scores %>% 
  as_tibble() %>%
  pull(Comp.1)

# add this pca axis score to the disp_axis data
disp_axis <- 
  disp_axis %>%
  mutate(dispersal_trait_axis = comp.1_scores)

# how is this new trait axis related to mean_ucrit and mean_pld
disp_axis %>%
  gather(mean_ucrit, mean_pld, key = "key", value = "value") %>%
  ggplot(mapping = aes(x = dispersal_trait_axis, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~key) +
  theme_classic()

# dispersal trait axis is positively related to both mean_pld and mean_ucrit values
# we know that this axis predicts ibd well in a subset of these species
# we assume that this axis predicts dispersal for other species as well

disp_axis



### calculate dietary diversity for each species

# load the diet taxon matrix
diet_mat <- read_csv( here("data/diet taxon.csv") )
str(diet_mat)

# load the diet data for each species
diet_spp <- read_csv( here("data/2_diet_specialisation2.csv") )
str(diet_spp)

# check species names of diet_spp
diet_spp$species %>% 
  unique() %>%
  sort()

# correct the incorrectly spelt names
diet_spp$species[which(diet_spp$species == "amphirion_clarki")] <- c("amphiprion_clarki")
diet_spp$species[which(diet_spp$species == "amphirion_melanopus")] <- c("amphiprion_melanopus")
diet_spp$species[which(diet_spp$species == "amphirion_percula")] <- c("amphiprion_percula")

# remove the acanthochromis_polyacanthus (direct developer)
diet_spp <- 
  diet_spp %>%
  filter(species != c("acanthochromis_polyacanthus"))

# there are two incorrect column names
diet_spp <- 
  diet_spp %>%
  dplyr::rename(plants_other_plants_terrestrial_plants = `plants_other_plants_terrestrial plants`,
         zooplankton_fish_early_stages_fish_eggs_larvae = `zooplankton_fish(early_stages)_fish_eggs_larvae`)

# make sure the diet names are harmonised
# if names are different give it a zero and then subset the zeros
bind_cols(mat = diet_mat$unique_id,
          spp = diet_spp %>% 
            select(-species) %>% 
            names()
) %>%
  mutate(harm = if_else(mat == spp, 1, 0)) %>%
  filter(harm == 0)

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


# convert diet_mat into a regular dataframe (due to the vegan package function)
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

# first, remove the rows without any diet data (i.e. the nas)
# second, remove the species column
# third conver to presence absence
diet_div_dat <- 
  diet_spp %>%
  filter_all( all_vars( (!is.na(.)) ) ) %>%
  select(-species) %>%
  as.data.frame() %>%
  decostand(., method = "pa")

diet_div <- 
  taxondive(diet_div_dat, diet_dis)

# attach these diet specialisation indices to a dataframe with the species names
spp_diet_div <- 
  diet_spp %>%
  filter_all( all_vars( (!is.na(.)) ) ) %>%
  select(species) %>%
  mutate(food_groups = diet_div$Species,
         diet_dplus = diet_div$Dplus)

# check the correlations among these variables
pairs(select(spp_diet_div, -species))

# this dataframe can now be joined to the giant main dataframe
spp_diet_div


### calculate habitat diversity for each species

# load the habitat taxon matrix
hab_mat <- read_csv( here("data/habitat taxon matrice.csv") )
str(hab_mat)

# load the habitat data for each species
hab_spp <- read_csv( here("data/3_habitat_specialisation.csv") )

# rename X1 column as species
hab_spp <- 
  hab_spp %>%
  dplyr::rename(species = "X1")
str(hab_spp)

# check the species names
hab_spp$species %>%
  unique() %>%
  sort()

# correct the incorrectly spelt names
hab_spp$species[which(hab_spp$species == "amphirion_clarki")] <- c("amphiprion_clarki")
hab_spp$species[which(hab_spp$species == "amphirion_melanopus")] <- c("amphiprion_melanopus")
hab_spp$species[which(hab_spp$species == "amphirion_percula")] <- c("amphiprion_percula")

# remove the acanthochromis_polyacanthus species as it is the direct developer
hab_spp <- 
  hab_spp %>%
  filter(species != c("acanthochromis_polyacanthus"))

# check that the names are harmonised
bind_cols(mat = hab_mat$id_number,
          spp = hab_spp %>% 
            select(-species) %>% 
            names()
) %>%
  mutate(harm = if_else(mat == spp, 1, 0)) %>%
  filter(harm == 0)

# these names are harmonised

# fill the NAs in the hab_mat so that it is unique
hab_mat <- hab_mat %>%
  group_by(habitat_i, habitat_ii) %>%
  mutate(habitat_iii = if_else(is.na(habitat_iii), "a", habitat_iii))

# as previously, we convert hab_mat into a regular dataframe (due to the vegan package function)
hab_mat <- as.data.frame(hab_mat)

# set row names for the hab_mat
row.names(hab_mat) <- hab_mat$id_number

# remove the unique_id column
hab_mat <- select(hab_mat, -id_number)

# set the food_i, food_ii and food_iii variables names as Family, Genus, Species
names(hab_mat) <- c("Family", "Genus", "Species")

# use the taxa2dist function to create a dissimilarity matrix
hab_dis <- taxa2dist(hab_mat, varstep = TRUE, check = TRUE)
hab_dis


# use the species data and the dissimilarity matrix to calculate taxonomic distinctness

# check for NAs
lapply(hab_spp, function(x) {sum(if_else(is.na(x), 1, 0)) })

lapply(hab_spp, function(x) { 
  bind_cols(id = hab_spp$species, na_row = is.na(x) ) %>%
    filter(na_row == TRUE) %>%
    pull(id) 
} )


# first, remove the rows without any habitat data (i.e. the nas)
# second, remove the species column
# third convert to presence absence
hab_div_dat <- 
  hab_spp %>%
  filter_all( all_vars( (!is.na(.)) ) ) %>%
  select(-species) %>%
  as.data.frame() %>%
  decostand(., method = "pa")

hab_div <- 
  taxondive(hab_div_dat, hab_dis)

# attach these habitat specialisation indices to a dataframe with the species names
spp_hab_div <- 
  hab_spp %>%
  filter_all( all_vars( (!is.na(.)) ) ) %>%
  select(species) %>%
  mutate(hab_groups = hab_div$Species,
         hab_dplus = hab_div$Dplus)

# check the correlations among these variables
pairs(select(spp_hab_div, -species))

# this dataframe can now be joined to the giant main dataframe
spp_hab_div



### join new dispersal and specialisation variables into the main dataset

# copy the main database without duplicates
fish_dat <- all_raw
str(fish_dat)

# let's see what data we have
fish_dat

spp_diet_div

spp_hab_div

disp_axis

# join the diet data and habitat data
match(spp_diet_div$species, spp_hab_div$species)

special <- 
  full_join(spp_diet_div, spp_hab_div, by = "species") %>%
  rename(binomial = "species")

# join the dispersal axis to the full data
disp <- 
  full_join(fish_dat, 
            select(disp_axis, -mean_ucrit, -mean_pld),
            by = c("family", "genus", "species", "binomial"))


# join the full data and the specialisation data
ds_dat <- 
  full_join(disp, special, by = "binomial")


### make a single specialisation axis for the environmental conditions

names(ds_dat)

# depth_range
# tpref_range
# salinity_range

spec_vars <- c("depth_range", "tpref_range", "salinity_range")

spec_dat <- 
  ds_dat %>%
  select(binomial, spec_vars)
  
# remove rows with NA's
spec_dat <- 
  spec_dat %>%
  filter_at(vars(spec_vars), all_vars(!is.na(.)))

pca_spec <- princomp(reformulate(spec_vars), data = spec_dat, cor = TRUE)
summary(pca_spec)

biplot(pca_spec)

# add the Comp.1 scores as a variable
spec_dat$env_special <- 
  pca_spec$scores %>%
  as_tibble() %>%
  pull(Comp.1)
  
# join this variable onto the ds_dat data
ds_dat <- 
  full_join(ds_dat, select(spec_dat, binomial, env_special), by = c("binomial"))


# you will work with this ds_dat dataset
# note that the important variables are:
# dispersal traits: dispersal_trait_axis, mean_ucrit, mean_pld
# dispersal_trait_axis (high = high dispersal, low = low dispersal)
# diet_dplus and hab_dplus (high = generalist, low = specialist)
# env_special (high = generalist, low = specialist)

ggplot(data = ds_dat %>% 
         group_by(family) %>%
         mutate(latitude_species = abs(latitude_species)) %>%
         summarise_at(vars(c("latitude_species", "hab_dplus")), ~mean(., na.rm = TRUE)),
       mapping = aes(x = latitude_species, y = hab_dplus) ) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

ggplot(data = ds_dat %>%
         filter(family == "pomacentridae"),
       mapping = aes(x = env_special, y = dispersal_trait_axis)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()


# now you can explore patterns in the data as recommended in the last meeting





