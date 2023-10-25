### This script is part of the available material associated with the article 
### "Photoacclimation in the kleptoplastidic ciliate Mesodinium rubrum and its 
### cryptophyte prey Teleaulax amphioxeia: phenotypic variability and implications
### for red tide remote sensing", submitted to the Journal of Plankton Research on 2023/09/29

### This script and all the associated data files are made available publicly to follow the
### FAIR (Findable, Accessible, Interoperable, Reusable) principles for scientific data.

### This script was written entirely by the authors of the article. Please give credit if reusing all or parts
### of the script

####### Script Photoacclim 2.0: Absorption_processing ###
### Date of last modification: 2023/10/25
# Author: Victor Pochic, ISOMer, Nantes Universite

# This scripts processes and plots absorption data of the photoacclimation experiment presented
# in the article.

### Packages needed:

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)

#### Importing data ####

### TECAN - PE extracts absorption ###
# For data obtained with the TECAN spectrophotometer, the data need to be first transformed into a clean
# .csv file

ODPE <- read.csv('OD/ODPE_photoacclim2.csv', sep = ';', header = TRUE, dec = ',')

### Perkin-Elmer - Live cells absorption ###

ODCuvette <- read.csv2('OD/OD_TeleMeso_cuvette_2022112829.csv', header = TRUE)

### Perkin-Elmer - Absorption on filters ###

ODFilters <- read.csv('OD/OD_TeleMeso_filter_20221215.csv', sep = ';', header = TRUE)
ODFilters2 <- read.csv('OD/OD_Meso_filter_20220429.csv', sep = ';', header = TRUE)

### HPLC data - pigment concentrations ###
# HPLC data needs to be adapted into mg.m-3 beforehand

HPLCdata <- read.csv('HPLC_photoacclim/HPLCdata_photoacclim2.csv', sep = ';', dec = ',', header = TRUE)

### Pigments from Annick (pigments spectra from Bricaud et al., 2003) ###

PAnnick <- read.csv('HPLC_photoacclim/Pigments_from_annick.csv', sep = ';', dec = ',', header = TRUE)

### Importing a datframe with a list of the samples ###

Sample_list <- read.csv('Experiment_samples/Samples_photoacclim2.csv', sep = ';', header = TRUE)

### Simple cell counts extracted from the file HPLCdata
Cell_count <- select(HPLCdata, c('Name', 'Cell_count', 'Species'))
# As soon as the cell counts are extracted, we exclude the Meso20N lines in the dataframe to avoid conflicts
HPLCdata <- filter(HPLCdata, Sample != 'Na')


##### Processing absorption from OD #####
#### For TECAN OD data ####

OD_to_abs.PE <- function(OD, l, Vs, Vb, p_surn) # l le chemin optique en m ; 
  # Vs et Vb les volumes d echantillon et de solvant, respetivement ; 
  # p_surn le pourcentage de cellules restant dans le surnageant apres centrifugation
{
  log(10)*((OD+OD*p_surn)/(Vs/Vb))/l
}
## Avec toutes les valeurs qui vont bien pour la fonction
# Chemin optique du tecan en m pour un volume = 200 microL --> calcule avec la methode expliquee dans le doc word
l_tecan <- 0.006
# Volume de l echantillon pour la PE
Vs <- 40
# Volume de solvant pour l extraction de PE
Vb <- 5
# Pourcentage de cellules restant dans le surnageant
p_surnTele = 0.148
p_surnMeso = 0

###
#### Computing absorption from OD
AbsPE <- ODPE %>%
  mutate_at(vars(-'Wavel.'), list(~ OD_to_abs.PE(., l_tecan, Vs, Vb, p_surnTele))) %>%
  mutate_at(vars(-c('Wavel.', 'C1', 'C2', 'C3')), 
            list(~ .-((C1 + C2 + C3)/3))) %>% # Subtracting the blank from the other spectra
  filter(Wavel. >= 410) %>%
  filter(Wavel. <= 700) # Adjusting lower and upper boundaries of the spectrum to match Pigments from Annick


#### For in vivo measurements with the Perkin-Elmer ####

OD_to_abs.cuvette <- function(OD, l) # l le chemin optique en m
{
  log(10)*(OD/l)
}

l_cuvette = 0.01 # The optical path in the cuvette is 1 cm

###

#### Computing absorption from OD

# First, we subtract the mean of OD between 750 and 850 nm -> absorption should be equal to 0 in this area
# We create ODCuvette_ghost, that will serve as a reference
ODCuvette_ghost <- ODCuvette
ODCuvette_cor <- ODCuvette

for (i in 2:length(colnames(ODCuvette_cor))) { # we start at 2 to avoid the wavelength column
  # This for loop subtracts the mean values btw 750 and 850 nm -> removing any offset
  ODCuvette_cor <- mutate_at(ODCuvette_cor, vars(colnames(ODCuvette_cor)[i]), 
                              list(~ as.numeric(.)-mean(as.numeric(ODCuvette_ghost[,i][1:101]))))
}

# To obtain a usable version of ODCuvette_cor, you need to apply the "renaming" and "removing blank" pipes
# below, but without the "conversion to absorption" part

AbsCuvette <- mutate_at(ODCuvette_cor, vars(-'nm'), list(~ OD_to_abs.cuvette(as.numeric(.), l_cuvette))) %>%
  mutate('Wavel.' = nm) %>%
  # Subtracting the mean of the blanks for Tele
  mutate_at(vars(-c('nm', 'Wavel.', 'Blank01', 'Blank02', 'Blank03',
                    'Blank04', 'Blank05', 'Blank06', 'Blank07',
                    'Blank08', 'Blank09', 'Blank10', 'Blank11',
                    'Autozero_2', 'Air_2', 'CDOM_MES20N_01', 'CDOM_MES80_02', 
                    'CDOM_MES200_03', 'Blank01_2', 'MES20N_R1', 'MES20N_R2', 
                    'MES20N_R3', 'MES80_R1', 'MES80_R2', 'MES80_R3', 'MES200_R1', 
                    'MES200_R2', 'MES200_R3')),
            list(~ .-((Blank01 + Blank02 + Blank03 +
                         Blank04 + Blank05 + Blank06 + Blank07 +
                         Blank08 + Blank09 + Blank10 + Blank11)/11))) %>%
  # And then for Meso
  mutate_at(vars(c('CDOM_MES20N_01', 'CDOM_MES80_02','CDOM_MES200_03','MES20N_R1', 'MES20N_R2', 'MES20N_R3', 'MES80_R1', 
                   'MES80_R2', 'MES80_R3', 'MES200_R1', 'MES200_R2', 'MES200_R3')),
            list(~ .-(Blank01_2))) %>%
  
# Subtracting the CDOM for Tele
mutate_at(vars(c('TEL20_R1_01', 'TEL20_R1_02', 'TEL20_R2_01', 'TEL20_R3_01', 'TEL20_R3_02')),
          list(~ .-((CDOM_TEL20_01 + CDOM_TEL20_02)/2))) %>%
mutate_at(vars(c('TEL80_R1_01', 'TEL80_R1_02', 'TEL80_R2_01', 'TEL80_R2_02', 'TEL80_R3_01', 'TEL80_R3_02')),
          list(~ .-((CDOM_TEL80_01 + CDOM_TEL80_02 + CDOM_TEL80_03)/3))) %>%
mutate_at(vars(c('TEL200_R1_01', 'TEL200_R2_01', 'TEL200_R3_01')),
          list(~ .-(CDOM_TEL200_02))) %>%
# Subtracting the CDOM for Meso
mutate_at(vars(c('MES20N_R1', 'MES20N_R2', 'MES20N_R3')),
          list(~ .-(CDOM_MES20N_01))) %>%
mutate_at(vars(c('MES80_R1', 'MES80_R2', 'MES80_R3')),
          list(~ .-(CDOM_MES80_02))) %>%
mutate_at(vars(c('MES200_R1', 'MES200_R2', 'MES200_R3')),
          list(~ .-(CDOM_MES200_03))) %>%

# Matching names with the other dataframes used
rename(Tele20_1 = TEL20_R1_01) %>%
  rename(Tele20_2 = TEL20_R2_01) %>%
  rename(Tele20_3 = TEL20_R3_01) %>%
  rename(Tele80_1 = TEL80_R1_01) %>%
  rename(Tele80_2 = TEL80_R2_01) %>%
  rename(Tele80_3 = TEL80_R3_01) %>%
  rename(Tele200_1 = TEL200_R1_01) %>%
  rename(Tele200_2 = TEL200_R2_01) %>%
  rename(Tele200_3 = TEL200_R3_01) %>%
  rename(Meso20_1 = MES20N_R1) %>%
  rename(Meso20_2 = MES20N_R2) %>%
  rename(Meso20_3 = MES20N_R3) %>%
  rename(Meso80_1 = MES80_R1) %>%
  rename(Meso80_2 = MES80_R2) %>%
  rename(Meso80_3 = MES80_R3) %>%
  rename(Meso200_1 = MES200_R1) %>%
  rename(Meso200_2 = MES200_R2) %>%
  rename(Meso200_3 = MES200_R3) %>%
  mutate_all(list(~ as.numeric(.))) %>%
  select(all_of(c(Sample_list$Name, 'Wavel.'))) # And select the relevant values

### Plotting the absorption spectra of some replicates (just to visualize)
ggplot(data = AbsCuvette) +
  geom_line(aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = 2)

# Standardizing the absorption value by the cell count
AbsCells <- AbsCuvette %>%
  mutate(Tele20_1 = Tele20_1/Cell_count$Cell_count[1]) %>%
  mutate(Tele20_2 = Tele20_2/Cell_count$Cell_count[2]) %>%
  mutate(Tele20_3 = Tele20_3/Cell_count$Cell_count[3]) %>%
  mutate(Tele80_1 = Tele80_1/Cell_count$Cell_count[4]) %>%
  mutate(Tele80_2 = Tele80_2/Cell_count$Cell_count[5]) %>%
  mutate(Tele80_3 = Tele80_3/Cell_count$Cell_count[6]) %>%
  mutate(Tele200_1 = Tele200_1/Cell_count$Cell_count[7]) %>%
  mutate(Tele200_2 = Tele200_2/Cell_count$Cell_count[8]) %>%
  mutate(Tele200_3 = Tele200_3/Cell_count$Cell_count[9]) %>%
  # Here we get the cell counts for MesoN, though we renamed accordingly to match the HPLC pigment data
  mutate(Meso20_1 = Meso20_1/Cell_count$Cell_count[10]) %>%
  mutate(Meso20_2 = Meso20_2/Cell_count$Cell_count[11]) %>%
  mutate(Meso20_3 = Meso20_3/Cell_count$Cell_count[12]) %>%
  mutate(Meso80_1 = Meso80_1/Cell_count$Cell_count[13]) %>%
  mutate(Meso80_2 = Meso80_2/Cell_count$Cell_count[14]) %>%
  mutate(Meso80_3 = Meso80_3/Cell_count$Cell_count[15]) %>%
  mutate(Meso200_1 = Meso200_1/Cell_count$Cell_count[16]) %>%
  mutate(Meso200_2 = Meso200_2/Cell_count$Cell_count[17]) %>%
  mutate(Meso200_3 = Meso200_3/Cell_count$Cell_count[18])

# Compute the means
AbsCells.means <- mutate(AbsCells, Tele20 = (Tele20_1 + Tele20_2 + Tele20_3)/3) %>%
  mutate(Tele80 = (Tele80_1 + Tele80_2 + Tele80_3)/3) %>%
  mutate(Tele200 = (Tele200_1 + Tele200_2 + Tele200_3)/3) %>%
  mutate(Meso20 = (Meso20_1 + Meso20_2 + Meso20_3)/3) %>%
  mutate(Meso80 = (Meso80_1 + Meso80_2 + Meso80_3)/3) %>%
  mutate(Meso200 = (Meso200_1 + Meso200_2 + Meso200_3)/3) %>%
  select(c('Wavel.', 'Tele20', 'Tele80', 'Tele200', 'Meso20', 'Meso80', 'Meso200')) %>%
  mutate_at(vars(-c('Wavel.')), list(~ ./(10^6))) %>% # Dividing by 10^6 because the cell
  # counts are in cells/mL -> so that we have absorption cross section in m2/cell
  filter(Wavel. >= 400 & Wavel. <= 751)

#### Supp plot : absorption spectra of T. amphioxeia and M. rubrum ####

AbsCells.std <- AbsCells.means %>%
  mutate(Meso80.std = Meso80/(Meso80[77])) %>%
  mutate(Tele80.std = Tele80/(Tele80[77]))
  
### Plotting
plot_abs_2sp <- ggplot(data = AbsCells.std) +
  geom_line(aes(x=Wavel., y = Meso80.std, color = 'M. rubrum'), linewidth = .75, linetype = 1) +
  geom_line(aes(x=Wavel., y = Tele80.std, color = 'T. amphioxeia'), linewidth = .75, linetype = 1) +
  labs(y = NULL, x = NULL, color = c(expression(paste("Species:"))))+
  # Format of x axis
  scale_x_continuous(breaks = seq(400, 750, by = 50)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,1.5, by = 0.5), limits = c(0,1.55)) +
  # Title and axis labels
  # ggtitle('Mesodinium absorption spectra') +
  labs(y = c(expression(paste(italic('a')['p'],'(',lambda,')/',italic('a')['p'],'(675) (-)'))), 
       x = (expression(paste('Wavelength (nm)')))) +
  scale_color_manual(name = "Species:", breaks = c('M. rubrum', 'T. amphioxeia'),
                     values = c('#440154','#21918c'),
                     labels = c((expression(paste(italic('M. rubrum')))),
                                (expression(paste(italic('T. amphioxeia'))))))+
  # Theme
  theme(plot.title = element_text(size = 9, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=9, color = 'black'),
        axis.text.y = element_text(size=9, color = 'black'),
        # Legend
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), 
        legend.position = 'bottom', legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray25', linewidth = .2), panel.grid.minor = element_line(color = 'gray70', linewidth = .12),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(nrow=1, byrow = TRUE))

plot_abs_2sp

# ggsave('abs_2sp.plot_supp.tiff', dpi = 600, width = 164, height = 100, units = 'mm')


#### For OD measurements on filters with the Perkin-Elmer ####

OD_to_abs.filter <- function(ODf, V, D)
{
  ODs  <- 0.323*ODf^1.0867 ## standard pathlength correction 
  V_m3 <- V*1e-6 ## Volume from ml to m3
  A_m2 <- (pi/4)*(D*1e-3)^2 ## Surface area in m2 with D in mm
  log(10)*ODs*A_m2/V_m3
}

# Parametres
# Parameters of filtered volume and colored area diameter are reported in csv files
Parameters_filters <- read.csv2('Optics_photoacclim/Parameters_filters_20221215.csv', header = TRUE)
Parameters_filters2 <- read.csv2('Optics_photoacclim/Parameters_filters_20220429.csv', header = TRUE)
# Volume filtre en mL
# Diametre de la tache de pigments en mm

#### For filters -> experiment at 21 degrees Celsius ####

# First, we subtract the mean of OD between 750 and 850 nm -> absorption should be equal to 0 in this area
# We create ODFilters_ghost, that will serve as a reference
ODFilters_ghost <- ODFilters
ODFilters_cor <- ODFilters

for (i in 2:length(colnames(ODFilters_cor))) { # we start at 2 to avoid the wavelength column
  ODFilters_cor <- mutate_at(ODFilters_cor, vars(colnames(ODFilters_cor)[i]), 
                         list(~ as.numeric(.)-mean(as.numeric(ODFilters_ghost[,i][1:101]))))
}

# To obtain a usable version of ODFilters_cor, you need to apply the "renaming" and "removing blank" pipes
# below, but not the "conversion to absorption" pipe

AbsFilters <- ODFilters_cor %>%
  # Subtracting the mean of the blanks
  mutate_at(vars(all_of(Parameters_filters$Sample)),
            list(~ .-((Blank00 + Blank01 + Blank02 + Blank03 +
                         Blank04 + Blank05 + Blank06)/7))) %>%
  select(all_of(c(Parameters_filters$Sample, 'nm'))) %>% # only mutating the samples spectra (not autozero, air or blank)
  rename(Wavel. = nm) # Renaming the wavelength vector to match other dataframes

# This for loop applies the function that computes absorption to all OD spectra, using the specific parameters recorded in
# the Parameters_filters file
for (i in 1:length(Parameters_filters$Sample)) {
  AbsFilters <- mutate_at(AbsFilters, vars(colnames(AbsFilters)[i]), 
                          list(~ OD_to_abs.filter(., Parameters_filters$Volume[i], Parameters_filters$Diameter[i])))
}


AbsFilters <- AbsFilters %>%
  # Rename all spectra with names that match those in other dataframes,
  # and computing the mean for the filters that were measured in 2 replicates because heterogenous.
  mutate(Meso20_1 = (MES20_R1_01+MES20_R1_02)/2) %>%
  rename(Meso20_2 = MES20_R2_01) %>%
  rename(Meso20_3 = MES20_R3_01) %>%
  rename(Meso80_1 = MES80_R1_01) %>%
  rename(Meso80_2 = MES80_R2_01) %>%
  rename(Meso80_3 = MES80_R3_01) %>%
  rename(Meso200_1 = MES200_R1_01) %>%
  rename(Meso200_2 = MES200_R2_01) %>%
  rename(Meso200_3 = MES200_R3_01) %>%
  rename(Tele20_1 = TEL20_R1_01) %>%
  rename(Tele20_2 = TEL20_R2_01) %>%
  mutate(Tele20_3 = (TEL20_R3_01+TEL20_R3_01)/2) %>%
  rename(Tele80_1 = TEL80_R1_01) %>%
  rename(Tele80_2 = TEL80_R2_01) %>%
  rename(Tele80_3 = TEL80_R3_01) %>%
  rename(Tele200_1 = TEL200_R1_01) %>%
  rename(Tele200_2 = TEL200_R2_01) %>%
  mutate(Tele200_3 = (TEL200_R3_01+TEL200_R3_01)/2) %>%
  select(all_of(c(Sample_list$Name, 'Wavel.'))) # And select the relevant values

### Plotting the absorption spectra of some replicates (just to visualize)
ggplot(data = AbsFilters) +
  geom_line(aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = 2)

### Comparing the absorption spectra of suspended cells and filters
# For Tele
ggplot() +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = .5) +
  geom_line(data = AbsCuvette, aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette, aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette, aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = .5, linetype = 2)

# For Meso
ggplot() +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = .5) +
  geom_line(data = AbsCuvette, aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette, aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette, aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = .5, linetype = 2)

# Standardizing absorption in cuvettes to match absorption on filters (by the absorption value at 675 nm)
AbsCuvette_cor <- select(AbsCuvette, all_of(c(Sample_list$Name, 'Wavel.')))

for (i in 1:length(Sample_list$Name)) {
  AbsCuvette_cor <- mutate_at(AbsCuvette_cor, vars(colnames(AbsCuvette_cor)[i]), 
                              list(~ .*(AbsFilters[,i][176]/AbsCuvette_cor[,i][176])))
}

# What does it give us?
# For Tele
ggplot() +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = .5) +
  geom_line(data = AbsCuvette_cor, aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette_cor, aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette_cor, aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = .5, linetype = 2)

# For Meso
ggplot() +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = .5) +
  geom_line(data = AbsFilters, aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = .5) +
  geom_line(data = AbsCuvette_cor, aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette_cor, aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCuvette_cor, aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = .5, linetype = 2)

# Computing the absorption per cell for filter samples, as we did for cuvette samples
# Standardizing the absorption value by the cell count
AbsCells.F <- AbsFilters %>%
  mutate(Tele20_1 = Tele20_1/Cell_count$Cell_count[1]) %>%
  mutate(Tele20_2 = Tele20_2/Cell_count$Cell_count[2]) %>%
  mutate(Tele20_3 = Tele20_3/Cell_count$Cell_count[3]) %>%
  mutate(Tele80_1 = Tele80_1/Cell_count$Cell_count[4]) %>%
  mutate(Tele80_2 = Tele80_2/Cell_count$Cell_count[5]) %>%
  mutate(Tele80_3 = Tele80_3/Cell_count$Cell_count[6]) %>%
  mutate(Tele200_1 = Tele200_1/Cell_count$Cell_count[7]) %>%
  mutate(Tele200_2 = Tele200_2/Cell_count$Cell_count[8]) %>%
  mutate(Tele200_3 = Tele200_3/Cell_count$Cell_count[9]) %>%
  # Here, in contrary to cuvette samples, we get the cell counts for Meso and not MesoN
  mutate(Meso20_1 = Meso20_1/Cell_count$Cell_count[13]) %>%
  mutate(Meso20_2 = Meso20_2/Cell_count$Cell_count[14]) %>%
  mutate(Meso20_3 = Meso20_3/Cell_count$Cell_count[15]) %>%
  mutate(Meso80_1 = Meso80_1/Cell_count$Cell_count[16]) %>%
  mutate(Meso80_2 = Meso80_2/Cell_count$Cell_count[17]) %>%
  mutate(Meso80_3 = Meso80_3/Cell_count$Cell_count[18]) %>%
  mutate(Meso200_1 = Meso200_1/Cell_count$Cell_count[19]) %>%
  mutate(Meso200_2 = Meso200_2/Cell_count$Cell_count[20]) %>%
  mutate(Meso200_3 = Meso200_3/Cell_count$Cell_count[21])

# Let's compare again
# For Tele
ggplot() +
  geom_line(data = AbsCells.F, aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = .5) +
  geom_line(data = AbsCells.F, aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = .5) +
  geom_line(data = AbsCells.F, aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = .5) +
  geom_line(data = AbsCells, aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCells, aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCells, aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = .5, linetype = 2)

# For Meso
ggplot() +
  geom_line(data = AbsCells.F, aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = .5) +
  geom_line(data = AbsCells.F, aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = .5) +
  geom_line(data = AbsCells.F, aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = .5) +
  geom_line(data = AbsCells, aes(x=Wavel., y = Meso20_1, color = 'Meso20_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCells, aes(x=Wavel., y = Meso80_1, color = 'Meso80_1'), linewidth = .5, linetype = 2) +
  geom_line(data = AbsCells, aes(x=Wavel., y = Meso200_1, color = 'Meso200_1'), linewidth = .5, linetype = 2)

#### For filters -> experiment at 17.5 degrees Celsius ####

# First, we subtract the mean of OD between 750 and 850 nm -> absorption should be equal to 0 in this area
# We create ODFilters_ghost, that will serve as a reference
ODFilters2_ghost <- ODFilters2

for (i in 2:length(colnames(ODFilters2))) { # we start at 2 to avoid the wavelength column
  ODFilters2 <- mutate_at(ODFilters2, vars(colnames(ODFilters2)[i]), 
                         list(~ as.numeric(.)-mean(as.numeric(ODFilters2_ghost[,i][1:101]))))
}

AbsFilters2 <- ODFilters2 %>%
  # Subtracting the mean of the blanks
  mutate_at(vars(all_of(Parameters_filters2$Sample)),
            list(~ .-((Blank01_01 + Blank01_02 + Blank02_01 + Blank02_02 + Blank03_01 + Blank03_02 
                       + Blank04_01 + Blank04_02 + Blank05_01 + Blank05_02 + Blank06_01 + Blank06_02)/12))) %>%
  select(all_of(c(Parameters_filters2$Sample, 'Wavel.'))) # only mutating the samples spectra (not autozero, air or blank)

# This for loop applies the function that computes absorption to all OD spectra, using the specific parameters recorded in
# the Parameters_filters file
for (i in 1:length(Parameters_filters2$Sample)) {
  AbsFilters2 <- mutate_at(AbsFilters2, vars(colnames(AbsFilters2)[i]), 
                           list(~ OD_to_abs.filter(., Parameters_filters2$Volume[i], Parameters_filters2$Diameter[i])))
}

AbsFilters2 <- AbsFilters2 %>%
  # Rename all spectra with names that match those in other dataframes,
  # and computing the mean for all the filters that were measured in 2 replicates
  mutate(Meso20_1a = (MES20_R1_01+MES20_R1_02)/2) %>%
  mutate(Meso20_2a = (MES20_R2_01+MES20_R2_02)/2) %>%
  mutate(Meso20_3a = (MES20_R3_01+MES20_R3_02)/2) %>%
  mutate(Meso65_1 = (MES65_R1_01+MES65_R1_02)/2) %>%
  mutate(Meso65_2 = (MES65_R2_01+MES65_R2_02)/2) %>%
  mutate(Meso65_3 = (MES65_R3_01+MES65_R3_02)/2) %>%
  mutate(Meso120_1 = (MES120_R1_01+MES120_R1_02)/2) %>%
  mutate(Meso120_2 = (MES120_R2_01+MES120_R2_02)/2) %>%
  mutate(Meso120_3 = (MES120_R3_01+MES120_R3_02)/2) %>%
  mutate(Meso220_1 = (MES220_R1_01+MES220_R1_02)/2) %>%
  mutate(Meso220_2 = (MES220_R2_01+MES220_R2_02)/2) %>%
  mutate(Meso220_3 = (MES220_R3_01+MES220_R3_02)/2) %>%
  select(all_of(c(Parameters_filters2$Name, 'Wavel.'))) # And select the relevant values

### Plotting the absorption spectra of some replicates (just to visualize)
ggplot(data = AbsFilters2) +
  geom_line(aes(x=Wavel., y = Meso20_1a, color = 'Meso20_1a'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso65_1, color = 'Meso65_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso120_1, color = 'Meso120_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso220_1, color = 'Meso220_1'), linewidth = 2)

##### Absorption of pigments in the sample #####

# Computing the spectra of pigments based on the HPLC concentration
Pigment_spectra <- data.frame(Wavel. = PAnnick$Wavelength, Tele20_1 = rep(NA, 146), 
                              Tele20_2 = rep(NA, 146), Tele20_3 = rep(NA, 146),
                              Tele80_1 = rep(NA, 146), 
                              Tele80_2 = rep(NA, 146), Tele80_3 = rep(NA, 146),
                              Tele200_1 = rep(NA, 146), 
                              Tele200_2 = rep(NA, 146), Tele200_3 = rep(NA, 146),
                              Meso20_1 = rep(NA, 146), 
                              Meso20_2 = rep(NA, 146), Meso20_3 = rep(NA, 146),
                              Meso80_1 = rep(NA, 146), 
                              Meso80_2 = rep(NA, 146), Meso80_3 = rep(NA, 146),
                              Meso200_1 = rep(NA, 146), 
                              Meso200_2 = rep(NA, 146), Meso200_3 = rep(NA, 146))

# Modifying the cell count for pigment data by removing the Meso20N rows
Cell_count_pigs <- filter(Cell_count, grepl('Meso20N', Name, fixed = TRUE) == FALSE)


for (i in 1:length(Sample_list$Name)) {
  Work <- PAnnick %>%
    select(Wavelength, Chl.a, Chlc12, Allox, alpha.car) %>%
    mutate(Chl.a = Chl.a*HPLCdata$Chl.a[i]*1000) %>% # We multiply by 1000 to have pigment concentrations
    mutate(Chl.c2 = Chlc12*HPLCdata$Chl.c2[i]*1000) %>% # in mg.m-3 (they were in mg.L-1)
    mutate(Allox = Allox*HPLCdata$Allox[i]*1000) %>%
    mutate(alpha.car = alpha.car*HPLCdata$alpha.car[i]*1000) %>%
    select(Wavelength, Chl.a, Chl.c2, Allox, alpha.car) %>%
    mutate(PE = AbsPE[,i+1]) %>% # Adding the absorption spectrum of PE for the sample
    mutate(full_spectrum = (Chl.a+Chl.c2+Allox+alpha.car+PE)/(Cell_count_pigs$Cell_count[i]*(10^6))) # We divide by the cell count 
  # in cells/m3 so we have absorption cross-section (in m2/cell)
  Pigment_spectra[,i+1] <- Work$full_spectrum
}

# Computing means for each condition
Pigment_spectra.means <- Pigment_spectra %>%
  mutate(Tele20 = (Tele20_1 + Tele20_2 + Tele20_3)/3) %>%
  mutate(Tele80 = (Tele80_1 + Tele80_2 + Tele80_3)/3) %>%
  mutate(Tele200 = (Tele200_1 + Tele200_2 + Tele200_3)/3) %>%
  mutate(Meso20 = (Meso20_1 + Meso20_2 + Meso20_3)/3) %>%
  mutate(Meso80 = (Meso80_1 + Meso80_2 + Meso80_3)/3) %>%
  mutate(Meso200 = (Meso200_1 + Meso200_2 + Meso200_3)/3) %>%
  select(c('Wavel.', 'Tele20', 'Tele80', 'Tele200', 'Meso20', 'Meso80', 'Meso200'))


### Plotting this sh*t
ggplot(data = Pigment_spectra) +
  geom_line(aes(x=Wavel., y = Tele20_1, color = 'Tele20_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele80_1, color = 'Tele80_1'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele200_1, color = 'Tele200_1'), linewidth = 2)


# Note : we don't include pigments such as crocoxanthin, phaeophytin etc... because they are not included
# in the dataset of Bricaud et al. (2003)
# We're not using the data from Clementson and Wojtasiewicz (2019) because these are absorption spectra
# measured in organic solvents
#### Computing the package effect ####

# Combining the absorption spectra of cells and pigments
Absorption_props <- AbsCells.means %>%
  mutate(Wavel. = as.integer(Wavel.)) %>% # For some reason Wavel. is still a character in the absorption dataframe
  left_join(Pigment_spectra.means, suffix = c('.cel', '.pig'), by = 'Wavel.') %>%
  filter(is.na(Tele20.pig) == FALSE) %>%
  mutate(Tele20.Qstar_a = Tele20.cel/Tele20.pig) %>% # Calculate the package effect for Tele
  mutate(Tele80.Qstar_a = Tele80.cel/Tele80.pig) %>%
  mutate(Tele200.Qstar_a = Tele200.cel/Tele200.pig) %>%
  mutate(Meso20.Qstar_a = Meso20.cel/Meso20.pig) %>% # Calculate the package effect for Meso
  mutate(Meso80.Qstar_a = Meso80.cel/Meso80.pig) %>%
  mutate(Meso200.Qstar_a = Meso200.cel/Meso200.pig)

ggplot(data = Absorption_props) +
  geom_line(aes(x=Wavel., y = Tele20.Qstar_a, color = 'Tele20.Qstar_a'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele80.Qstar_a, color = 'Tele80.Qstar_a'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Tele200.Qstar_a, color = 'Tele200.Qstar_a'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso20.Qstar_a, color = 'Meso20.Qstar_a'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso80.Qstar_a, color = 'Meso80.Qstar_a'), linewidth = 2) +
  geom_line(aes(x=Wavel., y = Meso200.Qstar_a, color = 'Meso200.Qstar_a'), linewidth = 2)

##### Calculating the specific absorption coefficients of chl a and PE 545 #####

AbsCoeffs <- AbsCuvette %>%
  melt(id.vars = 'Wavel.', variable.name = 'Name', value.name = 'Abs') %>%
  filter(Wavel. == 675 | Wavel. == 552 | Wavel. == 443 | Wavel. == 665) %>%
  filter(Name %in% Sample_list$Name)
# We created a dataframe with the absorption values at the wavelengths of interest

# Simplified dataframe of pigment concentrations in samples, computed from HPLC data
Pigment_data <- read.csv2('HPLC_photoacclim/Pigments_conc_photoacclim2.csv')

AbsCoeffs2 <- left_join(AbsCoeffs, Pigment_data, by = 'Name') %>%
  mutate(sigma_a = (Abs/Cell_count/(1000000))) %>% # Calculating absorption cross section
  # by dividing the absorption by the cell count /!\ IN CELLS/M3 /!\
  mutate(aphy_star = ifelse(Wavel. == 552, Abs/(PE), Abs/(Chl.a))) %>%
  # Calculating specific absorption coefficient of chl.a (443 and 675 nm) and
  # PE 545 (552 nm) using the pigment concentration /!\ in g/m3 /!\ (=mg/L)
  mutate(Temp = 21) # Adding the temperature

# Separating wavelengths for plots
AbsCoeffs.443 <- filter(AbsCoeffs2, Wavel. == 443)
AbsCoeffs.552 <- filter(AbsCoeffs2, Wavel. == 552)
AbsCoeffs.665 <- filter(AbsCoeffs2, Wavel. == 665)
AbsCoeffs.675 <- filter(AbsCoeffs2, Wavel. == 675)

# Computing stats
AbsCoeffs_means <- AbsCoeffs2 %>%
  group_by(Temp, Species, Irradiance, Wavel.) %>%
  summarise_at(vars('sigma_a', 'aphy_star'), list( ~ mean(.)), .groups = 'keep')

AbsCoeffs_sds <- AbsCoeffs2 %>%
  group_by(Temp, Species, Irradiance, Wavel.) %>%
  summarise_at(vars('sigma_a', 'aphy_star'), list( ~ sd(.)), .groups = 'keep')

AbsCoeffs_stats <- left_join(AbsCoeffs_means, AbsCoeffs_sds, by = c('Temp', 'Species', 'Irradiance', 'Wavel.'), 
                             suffix = c('_mean', '_sd'))

AbsCoeffs_stats.443 <- filter(AbsCoeffs_stats, Wavel. == 443)

# Formatting data for plotting - 552 nm
AbsCoeffs_stats.552 <- filter(AbsCoeffs_stats, Wavel. == 552) %>%
  mutate_at(vars(c('aphy_star_mean', 'aphy_star_sd')), list(~ ./1000)) %>% # Matching the units
  # The units of the first experiment for aphy_star is m2.mg-1, whereas it is in m2.g-1 for the second experiment
  mutate_at(vars('sigma_a_mean', 'sigma_a_sd'), list(~ ifelse(Species == 'Meso', ./50, .))) # Dividing the sigma_a for Meso by 50 for plotting

AbsCoeffs.552 <- mutate(AbsCoeffs.552, aphy_star = aphy_star/1000) %>% # Transforming the unit in m2.mg-1
  mutate_at(vars('sigma_a'), list(~ ifelse(Species == 'Meso', ./50, .))) 
# Dividing the sigma_a for Meso by 50 for plotting

# Formatting data for plotting - 665 nm
AbsCoeffs_stats.665 <- filter(AbsCoeffs_stats, Wavel. == 665) %>%
  mutate_at(vars(c('aphy_star_mean', 'aphy_star_sd')), list(~ ./1000)) %>% # Matching the units
  # The units of the first experiment for aphy_star is m2.mg-1, whereas it is in m2.g-1 for the second experiment
  mutate_at(vars('sigma_a_mean', 'sigma_a_sd'), list(~ ifelse(Species == 'Meso', ./50, .))) # Dividing the sigma_a for Meso by 50 for plotting

AbsCoeffs.665 <- mutate(AbsCoeffs.665, aphy_star = aphy_star/1000) %>% # Transforming the unit in m2.mg-1
  mutate_at(vars('sigma_a'), list(~ ifelse(Species == 'Meso', ./50, .))) 
# Dividing the sigma_a for Meso by 50 for plotting

# Formatting data for plotting - 675 nm
AbsCoeffs_stats.675 <- filter(AbsCoeffs_stats, Wavel. == 675) %>%
  mutate_at(vars(c('aphy_star_mean', 'aphy_star_sd')), list(~ ./1000)) %>% # Matching the units
  # The units of the first experiment for aphy_star is m2.mg-1, whereas it is in m2.g-1 for the second experiment
  mutate_at(vars('sigma_a_mean', 'sigma_a_sd'), list(~ ifelse(Species == 'Meso', ./50, .))) # Dividing the sigma_a for Meso by 50 for plotting

AbsCoeffs.675 <- mutate(AbsCoeffs.675, aphy_star = aphy_star/1000) %>% # Transforming the unit in m2.mg-1
  mutate_at(vars('sigma_a'), list(~ ifelse(Species == 'Meso', ./50, .))) 
# Dividing the sigma_a for Meso by 50 for plotting

# Saving the data as csv files
# write.csv2(AbsCoeffs_stats.675, 'AbsCoeffs_stats.675.csv', row.names = FALSE)
# write.csv2(AbsCoeffs_stats.552, 'AbsCoeffs_stats.552.csv', row.names = FALSE)

#### Plotting absorption coefficients ####

# Plotting the specific absorption coefficient for PE 545
Aphy_star.552.plot <- ggplot(AbsCoeffs_stats.552, aes(x = Irradiance, y = aphy_star_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = AbsCoeffs.552, aes(x = Irradiance, y = aphy_star, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,200, by = 50), limits = c(0,205)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.006, by = 0.002), limits = c(0,0.0061)) +
  # Title and axis labels
  ggtitle('D') +
  labs(y = c(expression(paste(~italic('a'['p']),'*(552 nm) (m'^'2','.mg'['PE 545'],''^'-1',')'))), x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')')))) +
  # Errorbars
  geom_errorbar(aes(ymin=aphy_star_mean-aphy_star_sd, ymax=aphy_star_mean+aphy_star_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Aphy_star.552.plot

sigma_a.552.plot <- ggplot(AbsCoeffs_stats.552, aes(x = Irradiance, y = sigma_a_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = AbsCoeffs.552, aes(x = Irradiance, y = sigma_a, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,200, by = 50), limits = c(0,205)) +
  scale_y_continuous(sec.axis = sec_axis(trans=~.*50,
                                         name = c(expression(paste(~italic('M. rubrum'),' ', ~italic(sigma['a']),'(552 nm) (m'^'2','.cell'^'-1',')'))))) +
  # Format of y axis
  # scale_y_continuous(breaks = seq(0,75, by = 25), limits = c(0,95)) +
  # Title and axis labels
  ggtitle('C') +
  labs(y = c(expression(paste(~italic('T. amphioxeia'),' ', ~italic(sigma['a']),'(552 nm) (m'^'2','.cell'^'-1',')'))), x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')')))) +
  # Errorbars
  geom_errorbar(aes(ymin=sigma_a_mean-sigma_a_sd, ymax=sigma_a_mean+sigma_a_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.title.y.right = element_text(angle = 90), # putting the secondary axis title in the same direction
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

sigma_a.552.plot

# Plotting the specific absorption coefficient for chl a - 675 nm
Aphy_star.675.plot <- ggplot(AbsCoeffs_stats.675, aes(x = Irradiance, y = aphy_star_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = AbsCoeffs.675, aes(x = Irradiance, y = aphy_star, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,200, by = 50), limits = c(0,205)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0.005,0.020, by = 0.005), limits = c(0.005,0.022)) +
  # Title and axis labels
  ggtitle('B') +
  labs(y = c(expression(paste(~italic('a'['p']),'*(675 nm) (m'^'2','.mg'['chl '~italic('a')],''^'-1',')'))),
       x = NULL) +
  # Errorbars
  geom_errorbar(aes(ymin=aphy_star_mean-aphy_star_sd, ymax=aphy_star_mean+aphy_star_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Aphy_star.675.plot

# Save the plot
# ggsave('Aphy_star.675.plot.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

sigma_a.675.plot <- ggplot(AbsCoeffs_stats.675, aes(x = Irradiance, y = sigma_a_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = AbsCoeffs.675, aes(x = Irradiance, y = sigma_a, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,200, by = 50), limits = c(0,205)) +
  scale_y_continuous(sec.axis = sec_axis(trans=~.*50,
                     name = c(expression(paste(~italic('M. rubrum'),' ', ~italic(sigma['a']),'(675 nm) (m'^'2','.cell'^'-1',')'))))) +
  # Format of y axis
  # scale_y_continuous(breaks = seq(0,75, by = 25), limits = c(0,95)) +
  # Title and axis labels
  ggtitle('A') +
  labs(y = c(expression(paste(~italic('T. amphioxeia'),' ', ~italic(sigma['a']),'(675 nm) (m'^'2','.cell'^'-1',')'))),
       x = NULL) +
       # Errorbars
  geom_errorbar(aes(ymin=sigma_a_mean-sigma_a_sd, ymax=sigma_a_mean+sigma_a_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.title.y.right = element_text(angle = 90), # putting the secondary axis title in the same direction
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

sigma_a.675.plot

# Save the plot
# ggsave('sigma_a.675.plot.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

# Arrange the plots and save
ggarrange(sigma_a.675.plot, Aphy_star.675.plot, sigma_a.552.plot, Aphy_star.552.plot, common.legend = TRUE, legend ='bottom')
# ggsave('AbsCoeffs.2wl.plot4.tiff', dpi = 600, width = 164, height = 164, units = 'mm')

##### For processing satellite data - visual checking of abs coeffs at 665 nm #####
# Plotting the specific absorption coefficient for chl a - 665 nm
Aphy_star.665.plot <- ggplot(AbsCoeffs_stats.665, aes(x = Irradiance, y = aphy_star_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = AbsCoeffs.665, aes(x = Irradiance, y = aphy_star, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,250, by = 50), limits = c(0,255)) +
  # Format of y axis
  # scale_y_continuous(breaks = seq(0,75, by = 25), limits = c(0,95)) +
  # Title and axis labels
  ggtitle('B') +
  labs(y = c(expression(paste(~italic('a'['ph']),'*(665 nm) (m'^'2','.mg'['chl '~italic('a')],''^'-1',')'))), x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')')))) +
  # Errorbars
  geom_errorbar(aes(ymin=aphy_star_mean-aphy_star_sd, ymax=aphy_star_mean+aphy_star_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Aphy_star.665.plot

# Save the plot
# ggsave('Aphy_star.675.plot.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

sigma_a.665.plot <- ggplot(AbsCoeffs_stats.665, aes(x = Irradiance, y = sigma_a_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = AbsCoeffs.665, aes(x = Irradiance, y = sigma_a, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,250, by = 50), limits = c(0,255)) +
  scale_y_continuous(sec.axis = sec_axis(trans=~.*50)) +
  # Format of y axis
  # scale_y_continuous(breaks = seq(0,75, by = 25), limits = c(0,95)) +
  # Title and axis labels
  ggtitle('A') +
  labs(y = c(expression(paste(~italic(sigma['a']),'(665 nm) (m'^'2','.cell'^'-1',')'))), x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')')))) +
  # Errorbars
  geom_errorbar(aes(ymin=sigma_a_mean-sigma_a_sd, ymax=sigma_a_mean+sigma_a_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

sigma_a.665.plot

#### Supp plot : irradiance spectra of culture light ####

# Import data
Meso_spectra <- read.csv('Optics_photoacclim/Mesodinium_irradiance_spectra2.csv', sep = ';', dec = '.', header = TRUE)
Tele_spectra <- read.csv('Optics_photoacclim/Teleaulax_irradiance_spectra2.csv', sep = ';', dec = '.', header = TRUE)

# Modifying the format of the data
# For Meso
Meso_spectra20 <- select(Meso_spectra, wl, Ed_20) %>%
  mutate(Species = 'Meso') %>%
  mutate(Ed_level = 20) %>%
  rename(Ed=Ed_20)

Meso_spectra80 <- select(Meso_spectra, wl, Ed_80) %>%
  mutate(Species = 'Meso') %>%
  mutate(Ed_level = 80) %>%
  rename(Ed=Ed_80)

Meso_spectra200 <- select(Meso_spectra, wl, Ed_200) %>%
  mutate(Species = 'Meso') %>%
  mutate(Ed_level = 200) %>%
  rename(Ed=Ed_200)

# And Tele
Tele_spectra20 <- select(Tele_spectra, wl, Ed_20) %>%
  mutate(Species = 'Tele') %>%
  mutate(Ed_level = 20) %>%
  rename(Ed=Ed_20)

Tele_spectra80 <- select(Tele_spectra, wl, Ed_80) %>%
  mutate(Species = 'Tele') %>%
  mutate(Ed_level = 80) %>%
  rename(Ed=Ed_80)

Tele_spectra200 <- select(Tele_spectra, wl, Ed_200) %>%
  mutate(Species = 'Tele') %>%
  mutate(Ed_level = 200) %>%
  rename(Ed=Ed_200)

### Joining the data together
Ed_spectra <- bind_rows(Meso_spectra20, Meso_spectra80, Meso_spectra200,
                        Tele_spectra20, Tele_spectra80, Tele_spectra200) %>%
  group_by(Species, Ed_level)

### Plotting
plot_Ed_spectra <- ggplot(data = Ed_spectra, aes(x = wl, y = Ed, color = factor(Ed_level), 
                                                 linetype = Species)) +
  geom_line(linewidth = .75) +
  labs(y = NULL, x = NULL, color = c(expression(paste("Irradiance level:"))),
       linetype = c(expression(paste("Species:"))))+
  # Format of x axis
  scale_x_continuous(breaks = seq(300, 1100, by = 100)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.1, by = 0.02), limits = c(0,0.1)) +
  # Title and axis labels
  ggtitle('B') +
  labs(y = c(expression(paste('Irradiance (W.m'^'-2','.nm'^'-1',')'))), 
       x = (expression(paste('Wavelength (nm)')))) +
  scale_color_manual(name = 'Irradiance level:', breaks = c('20', '80', '200'),
                     values = c('#440154','#21918c','#fde725'),
                     labels = c('20', '80', '200'))+
  scale_linetype_manual(name = 'Species:', breaks = c('Meso', 'Tele'),
                        values = c(1,2),
                        labels = c((expression(paste(italic('M. rubrum')))),
                                   (expression(paste(italic('T. amphioxeia'))))))+
  # Theme
  theme(plot.title = element_text(size = 9, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=9, color = 'black'),
        axis.text.y = element_text(size=9, color = 'black'),
        # Legend
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), 
        legend.position = 'bottom', legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major.y = element_line(colour = 'gray10', linewidth = .5), panel.grid.minor = NULL,
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
        # Artificial axes
        geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
        geom_vline(aes(xintercept = 300), color = 'black', linewidth = .5) +
        guides(color = guide_legend(nrow=1, byrow = TRUE), linetype = guide_legend(nrow=1, byrow = TRUE))

plot_Ed_spectra

# Saving the plot
# ggsave('Ed_spectra.plot_supp.tiff', dpi = 600, width = 164, height = 93, units = 'mm')
