### This script is part of the available material associated with the article 
### "Photoacclimation in the kleptoplastidic ciliate Mesodinium rubrum and its 
### cryptophyte prey Teleaulax amphioxeia: phenotypic variability and implications
### for red tide remote sensing", submitted to the Journal of Plankton Research on 2023/09/29

### This script and all the associated data files are made available publicly to follow the
### FAIR (Findable, Accessible, Interoperable, Reusable) principles for scientific data.

### This script was written entirely by the authors of the article. Please give credit if reusing all or parts
### of the script

####### Script Photoacclim 2.0: Optics and remote sensing ###
### Date of last modification: 2023/10/25
# Author: Victor Pochic, ISOMer, Nantes Universite

# This scripts processes and plots data of optics and remote sensing presented in the article.

### Packages needed:

library(dplyr)
library(viridis)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(cmocean)

#### Field data ####
## Gons algorithm
gons_inversion_aphy <- function(Rrs_B4, Rrs_B5, Rrs_B7){
  Rrs_B4    <- as.numeric(paste(Rrs_B4))
  Rrs_B5    <- as.numeric(paste(Rrs_B5))
  Rrs_B7    <- as.numeric(paste(Rrs_B7))
  bb_gons   <- 1.56*pi*Rrs_B7/(0.082 - 0.6*pi*Rrs_B7)
  aphy_gons <- (0.7+ bb_gons)*Rrs_B5/Rrs_B4 - 0.4 - bb_gons^1.02
  aphy_gons
}

## Import data
field_HPLC <- read.csv2("HPLC_photoacclim/20210329_pigments_HPLC_all_export.csv", header = TRUE)
field_Rrs <- read.csv2("Optics_photoacclim/20210329_radiometry_Laura_averaged_convolution_S2_transposed.csv", header = TRUE)

field_data <- left_join(field_HPLC, field_Rrs, by = 'Station') %>%
  mutate(aphy_665 = gons_inversion_aphy(Rrs_665, Rrs_705, Rrs_783)) %>% # compute aphy(665nm) from reflectance
  # with gons inversion algorithm
  mutate_at(vars(-c('Station')), list(~ as.numeric(.))) %>%
  filter(Station != 'PT09') # remove point 9

#### AbsCuvette vs AbsFilter ####

AbsMeso_Cuvette_Filter <- read.csv2('Optics_photoacclim/AbsMeso_Cuvette_Filter.csv', header = TRUE) %>%
  group_by(Irradiance)

# Defining the function of pathlength amplification correction by Stramski et al. (2015)
# to use in the plot
Stramski_plcor <- function(x){
  0.323*(x^1.0867)
}

# And the custom function of pathlength amplification correction fitted on 
# our data
Homemade_plcor <- function(x){
  0.2927*(x^0.8036)
}

# And the custom functions of pathlength amplification correction fitted on 
# our data
# The correction from the mean of all 9 samples
Homemade_plcor_mean <- function(x){
  0.3203*(x^0.7743)
}
# The lowest correction (mean-std error)
Homemade_plcor_min <- function(x){
  0.2388*(x^0.7194)
}
# The highest correction (mean+std error)
Homemade_plcor_max <- function(x){
  0.4018*(x^0.8354)
}

# Plotting

ggplot(data = AbsMeso_Cuvette_Filter, aes(x = Abs_Filter, y = Abs_Cuvette, 
                                          shape = factor(Irradiance), color = factor(Irradiance),
                                          fill = factor(Irradiance))) +
  geom_point(size = 1, alpha = .4) +
  scale_shape_manual(name = 'Irradiance:', labels = c('20', '80', '200'), 
                     values = c(21, 24, 23)) +
  scale_color_manual(name = 'Irradiance:', labels = c('20', '80', '200'), 
                     values = c('#440154', '#21918c', 'black')) +
  scale_fill_manual(name = 'Irradiance:', labels = c('20', '80', '200'), 
                    values = c('#440154', '#21918c', '#fde725')) +
    # Adding the standard pathlength correction of Stramski et al. (2015)
  geom_function(fun = Stramski_plcor, color = 'red', linewidth = .5) +
  # Adding the homemade pathlength corrections
  geom_function(fun = Homemade_plcor_min, color = 'darkred', linewidth = .5, linetype = 2) +
  geom_function(fun = Homemade_plcor_max, color = 'darkred', linewidth = .5, linetype = 2) +
  geom_function(fun = Homemade_plcor_mean, color = 'darkred', linewidth = .5) +
  # geom_function(fun = Homemade_plcor, color = 'black', linewidth = .5, linetype = 2) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,60, by = 10), limits = c(0,60)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,8, by = 1), limits = c(0,8)) +
  # Title and labels
  ggtitle('A') +
  labs(y = c(expression(paste('a'['p']^'uncorr',' (Cuvette) (m'^'-1',')'))), 
       x = (expression(paste('a'['ph'],' (Filter) (m'^'-1',')')))) +
  # Artificial axes
  geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
  geom_vline(aes(xintercept = 0), color = 'black', linewidth = .5) +
  # Theme
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.85, 0.32), legend.box.margin = margin(0, 0, 0, 0, 'cm'), 
        legend.background = element_rect(color = 'gray10'), legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray70', linewidth = .5),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides()


# ggsave('aph_cuvette_vs_filter.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

### What if we do the same thing but only with means for each irradiance ?

AbsMeso_Cuvette_Filter_stats <- read.csv2('Optics_photoacclim/AbsMeso_Cuvette_Filter_stats.csv', header = TRUE) %>%
  group_by(Irradiance)

# Plotting

plot_aph_cuvette_vs_filter <- ggplot(data = AbsMeso_Cuvette_Filter_stats, aes(x = Abs_Filter_mean, y = Abs_Cuvette_mean, 
                                          shape = factor(Irradiance),
                                          fill = Wavel.)) +
  geom_point(size = 3, alpha = .55) +
  scale_shape_manual(name = 'Irradiance:', labels = c('20', '80', '200'), values = c(21, 24, 23)) +
  scale_fill_viridis(name = 'Wavelength:') +
  # Adding the standard pathlength correction of Stramski et al. (2015)
  geom_function(fun = Stramski_plcor, color = 'red', linewidth = .5) +
  # Adding the homemade pathlength corrections
  geom_function(fun = Homemade_plcor_min, color = 'darkred', linewidth = .5, linetype = 2) +
  geom_function(fun = Homemade_plcor_max, color = 'darkred', linewidth = .5, linetype = 2) +
  geom_function(fun = Homemade_plcor_mean, color = 'darkred', linewidth = .5) +
  # geom_function(fun = Homemade_plcor, color = 'black', linewidth = .8, linetype = 2) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,50, by = 10), limits = c(0,50)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,7, by = 1), limits = c(0,7)) +
  # Title and labels
  ggtitle('A') +
  labs(y = c(expression(paste(italic('a')['p'],' (Cuvette) (m'^'-1',')'))), 
       x = c(expression(paste(italic('a')['p']^'uncorr',' (Filter) (m'^'-1',')')))) +
  # Artificial axes
  geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
  geom_vline(aes(xintercept = 0), color = 'black', linewidth = .5) +
  # Theme
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = 'right', legend.box.margin = margin(0, 0, 0, 0, 'cm'), 
        legend.background = element_rect(color = 'gray10'), legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray70', linewidth = .5),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides()

plot_aph_cuvette_vs_filter

#ggsave('aph_cuvette_vs_filter_stats3.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

#### AbsMeasured vs AbsModelled ####

ap_measured <- read.csv2('Optics_photoacclim/20210329_absorption_ap_avg_spectra.csv', header = TRUE)

ap_measured.665 <- filter(ap_measured, wl == 665)
# write.csv2(ap_measured.665, 'ap_measured.665.csv', row.names = FALSE)
# A little processing with excel and we have something nice
ap_measured.665_ordered <- read.csv2('Optics_photoacclim/ap_measured.665_ordered.csv', header = TRUE) %>%
  filter(Station != 'PT09') %>% # Removing point 9 because it is dubious data
  select(c('Station', 'aphy_665', 'aphy_665_homemade')) %>% # Selecting only relevant values for our plot
  mutate(aphy_665 = as.numeric(aphy_665)) %>%
  mutate(aphy_665_homemade = as.numeric(aphy_665_homemade))
  

aphy_665_comparison <- select(field_data, c('Station', 'aphy_665')) %>%
  left_join(by = 'Station', ap_measured.665_ordered, suffix = c('.modelled', '.measured'))

#### Evaluation of the influence of aphy_star ####

# Function for calculating chl a from aphy and aphy_star
modelled_chla <- function(aphy_665, aphy_star_665){
  chla <- aphy_665/aphy_star_665
  chla
}

# Creating aphy_star vector for the plot
aphy_star_vector <- as.data.frame(seq(0.005,0.025, by = 0.0001))
chla_vector <- as.data.frame(seq(0,349, by = 350/201))
vectors_plot <- bind_cols(aphy_star_vector, chla_vector)
colnames(vectors_plot) <- c('aphy_star', 'chla')
vectors_plot <- vectors_plot %>%
  mutate(chla_aphy0.2 = 0.2/aphy_star) %>%
  mutate(chla_aphy1 = 1/aphy_star) %>%
  mutate(chla_aphy3 = 3/aphy_star)

vectors_plot_table <- vectors_plot %>%
  filter(((aphy_star == 0.0071) | (aphy_star == 0.0111)) |
                       ((aphy_star == 0.0133) | (as.character(aphy_star) == '0.0146'))
                        | ((as.character(aphy_star)) == '0.0217'))
# write.csv2(vectors_plot_table, file = 'eval_aphy_star_table.csv', row.names = FALSE)

# And creating labels for the 3 different curves
lab_aphy0.2 <- as.character(c(expression(paste(italic('a')['phy'],'(665) = 0.2'))))
lab_aphy1 <- as.character(c(expression(paste(italic('a')['phy'],'(665) = 1'))))
lab_aphy3 <- as.character(c(expression(paste(italic('a')['phy'],'(665) = 3'))))


plot_eval_aphy_star <- ggplot(data = vectors_plot, aes(x = aphy_star, y = chla)) +
  # 3 curves for 3 aphy_star with the function of chla modelling
  geom_line(aes(x=aphy_star, y = chla_aphy0.2), color = '#65B769', linewidth = .75) +
  geom_line(aes(x=aphy_star, y = chla_aphy1), color = '#0B7748', linewidth = .75) +
  geom_line(aes(x=aphy_star, y = chla_aphy3), color = '#122414', linewidth = .75) +
  # Vertical lines to illustrate different aphy_star values
  geom_vline(aes(xintercept = 0.0073), color = 'dodgerblue3', linewidth = .5, linetype = 2) +
  geom_vline(aes(xintercept = 0.0111), color = 'dodgerblue3', linewidth = .5, linetype = 2) +
  geom_vline(aes(xintercept = 0.0144), color = 'midnightblue', linewidth = .5, linetype = 2) +
  geom_vline(aes(xintercept = 0.0146), color = 'gray25', linewidth = .5, linetype = 1) +
  #geom_vline(aes(xintercept = 0.021656), color = 'red', linewidth = .5, linetype = 1) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0.005,0.025, by = 0.005), limits = c(0.005,0.025)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,600, by = 50), limits = c(0,450)) +
  # Title and labels
  ggtitle('B') +
  labs(y = c(expression(paste('[Chl ', italic('a'), '] (mg.m'^'-3',')'))), 
       x = (expression(paste(italic('a')['phy']^'*','(665) (m'^'2','.mg'['chl '~italic('a')]^'-1',')')))) +
  # Artificial axes
  geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
  geom_vline(aes(xintercept = 0.005), color = 'black', linewidth = .5) +
  # Text
  geom_text(aes(x = 0.0182, y = 28, label = lab_aphy0.2), 
                color = '#65B769', parse = TRUE, size = 3) +
  geom_text(aes(x = 0.018, y = 90, label = lab_aphy1), 
            color = '#0B7748', parse = TRUE, size = 3) +
  geom_text(aes(x = 0.018, y = 215, label = lab_aphy3), 
            color = '#122414', parse = TRUE, size = 3) +
  # Theme
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.85, 0.32), legend.box.margin = margin(0, 0, 0, 0, 'cm'), 
        legend.background = element_rect(color = 'gray10'), legend.key = element_blank(),
        # Panel
        panel.grid.major.y = element_line(color = 'gray70', linewidth = .5),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides()

plot_eval_aphy_star

# ggsave('evaluation_aphystar.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

#### Abs Cuvette vs Filter : spectral comparison ####

### Importing OD data ###
ODCuvette <- read.csv2('OD/OD_TeleMeso_cuvette_2022112829.csv', header = TRUE)
ODFilters <- read.csv('OD/OD_TeleMeso_filter_20221215.csv', sep = ';', header = TRUE)

###Computing absorption from the OD data ---> For cuvette samples

# Function that computes absorption from OD ofr cuvette samples
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

AbsCuvette_comp <- mutate(AbsCuvette, Meso20_mean = (Meso20_1+Meso20_2+Meso20_3)/3) %>%
  mutate(Meso80_mean = (Meso80_1+Meso80_2+Meso80_3)/3) %>%
  mutate(Meso200_mean = (Meso200_1+Meso200_2+Meso200_3)/3) %>%
  select(all_of(c('Meso20_mean', 'Meso80_mean', 'Meso200_mean', 'Wavel.'))) %>%
  filter((Wavel. <= 750) & (Wavel. >= 400))

###Computing absorption from the OD data ---> For cuvette samples

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

# Function that computes absorption from OD for filter samples, modified to remove the pathlength correction
OD_to_abs.filter_nopl <- function(ODf, V, D)
{
  ODs  <- ODf ## Applying no pathlength correction 
  V_m3 <- V*1e-6 ## Volume from ml to m3
  A_m2 <- (pi/4)*(D*1e-3)^2 ## Surface area in m2 with D in mm
  log(10)*ODs*A_m2/V_m3
}

### Computing absorption from OD (filters) -> experiment at 21 degrees Celsius
AbsFilters_nopl <- ODFilters_cor %>%
  # Subtracting the mean of the blanks
  mutate_at(vars(all_of(Parameters_filters$Sample)),
            list(~ .-((Blank00 + Blank01 + Blank02 + Blank03 +
                         Blank04 + Blank05 + Blank06)/7))) %>%
  select(all_of(c(Parameters_filters$Sample, 'nm'))) %>% # only mutating the samples spectra (not autozero, air or blank)
  rename(Wavel. = nm) # Renaming the wavelength vector to match other dataframes

# This for loop applies the function that computes absorption to all OD spectra, using the specific parameters recorded in
# the Parameters_filters file
for (i in 1:length(Parameters_filters$Sample)) {
  AbsFilters_nopl <- mutate_at(AbsFilters_nopl, vars(colnames(AbsFilters_nopl)[i]), 
                               list(~ OD_to_abs.filter_nopl(., Parameters_filters$Volume[i], Parameters_filters$Diameter[i])))
}


AbsFilters_nopl <- AbsFilters_nopl %>%
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

AbsFilter_comp <- mutate(AbsFilters_nopl, Meso20_mean = (Meso20_1+Meso20_2+Meso20_3)/3) %>%
  mutate(Meso80_mean = (Meso80_1+Meso80_2+Meso80_3)/3) %>%
  mutate(Meso200_mean = (Meso200_1+Meso200_2+Meso200_3)/3) %>%
  select(all_of(c('Meso20_mean', 'Meso80_mean', 'Meso200_mean', 'Wavel.'))) %>%
  mutate_at(vars(c('Meso20_mean', 'Meso80_mean', 'Meso200_mean')), list(~ Homemade_plcor(.))) %>%
  filter((Wavel. <= 750) & (Wavel. >= 400))

Abs_comp <- left_join(AbsCuvette_comp, AbsFilter_comp, by = 'Wavel.', 
                      suffix = c('.cuvette', '.filter'))

### End of data processing

### Plotting
plot_Meso_abs <- ggplot(data = Abs_comp) +
  geom_line(aes(x=Wavel., y = Meso80_mean.cuvette, color = 'Cuvette'), linewidth = .75, linetype = 1) +
  geom_line(aes(x=Wavel., y = Meso80_mean.filter, color = 'Filter'), linewidth = .75, linetype = 1) +
  labs(y = NULL, x = NULL, color = c(expression(paste("Measurement technique:"))))+
  # Format of x axis
  scale_x_continuous(breaks = seq(400, 750, by = 50)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,7, by = 2), limits = c(0,7)) +
  # Title and axis labels
  # ggtitle('Mesodinium absorption spectra') +
  labs(y = c(expression(paste(a[p],' (m'^'-1',')'))), 
       x = (expression(paste('Wavelength (nm)')))) +
  scale_color_manual(name = "Measurement technique:", breaks = c('Cuvette', 'Filter'),
                     values = c('gray65','black'),
                     labels = c('Cuvette', 'Filter'))+
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
  guides(color = guide_legend(nrow=2, byrow = TRUE))

plot_Meso_abs

# ggsave('abscvf_Meso80_mean.plot.tiff', dpi = 600, width = 80, height = 100, units = 'mm')

#### Comparing different aphy_star values ####

field_data <- left_join(field_data, aphy_665_comparison, by = 'Station')

# Performing linear regression on modelled aphy vs chl a
lm_aphy_modelled_vs_chla <- lm(data = field_data, 
                      formula = field_data$aphy_665.modelled ~ field_data$Chlorophyll_a)
lm_aphy_modelled_vs_chla # There is an intercept different from 0
# We now do a modified lm with intercept fixed to 0
lm_aphy_modelled_vs_chla_0 <- lm(data = field_data, 
                        formula = field_data$aphy_665.modelled ~ 
                          0 + field_data$Chlorophyll_a)
summary(lm_aphy_modelled_vs_chla_0)

# Same linear regression with filter-pad aphy vs chl a (standard pathlength correction, Stramski et al. 2015)
lm_aphy_measured_vs_chla <- lm(data = field_data, 
                               formula = field_data$aphy_665.measured ~ field_data$Chlorophyll_a)
lm_aphy_measured_vs_chla # There is an intercept different from 0
# We now do a modified lm with intercept fixed to 0
lm_aphy_measured_vs_chla_0 <- lm(data = field_data, 
                                 formula = field_data$aphy_665.measured ~ 
                                   0 + field_data$Chlorophyll_a)
summary(lm_aphy_measured_vs_chla_0)

# Same linear regression with filter-pad aphy vs chl a (specific pathlength correction)
lm_aphy_homemade_vs_chla <- lm(data = field_data, 
                               formula = field_data$aphy_665_homemade ~ field_data$Chlorophyll_a)
lm_aphy_homemade_vs_chla # There is an intercept different from 0
# We now do a modified lm with intercept fixed to 0
lm_aphy_homemade_vs_chla_0 <- lm(data = field_data, 
                                 formula = field_data$aphy_665_homemade ~ 
                                   0 + field_data$Chlorophyll_a)
summary(lm_aphy_homemade_vs_chla_0)

#### Plots of aphy_vs chl a : several versions

## graph inversed aphy vs Chl full
plot_aph_vs_chla_field_data_full <- ggplot(data = field_data, aes(x = Chlorophyll_a, y = aphy_665)) +
  # points of 8 stations with aphy_modelled
  geom_point(color = 'midnightblue', shape = 1, size = 3, stroke = .75) +
  # 2 plcors for measured aphy represented with points
  geom_point(aes(x = Chlorophyll_a, y = aphy_665.measured), color = 'red', shape = 0, size = 3, stroke = .75) +
  geom_point(aes(x = Chlorophyll_a, y = aphy_665_homemade), color = 'darkred', shape = 2, size = 3, stroke = .75) +
  # 2 aphy_star from measured aphy
  geom_abline(aes(intercept=0, slope=lm_aphy_measured_vs_chla_0$coefficients), color ='red', linewidth = .5, linetype = 2) +
  geom_abline(aes(intercept=0, slope=lm_aphy_homemade_vs_chla_0$coefficients), color ='darkred', linewidth = .5, linetype = 2) +
  # aphy_star from Gons et al. (2002) and a_sol_chla from Bricaud et al. (2004)
  geom_abline(aes(intercept=0, slope=0.0146), color ='gray25', linewidth = .6, linetype = 1) +
  geom_abline(aes(intercept=0, slope=0.021656), color ='black', linewidth = .6, linetype = 1) +
  geom_abline(aes(intercept=0, slope=lm_aphy_modelled_vs_chla_0$coefficients), color ='midnightblue', linewidth = .5, linetype = 2) +
  geom_abline(aes(intercept=0, slope=0.0073), color ='dodgerblue3', linewidth = .5, linetype = 2) +
  geom_abline(aes(intercept=0, slope=0.0111), color ='dodgerblue3', linewidth = .5, linetype = 2) +
  
  # Format of x axis
  scale_x_continuous(breaks = seq(0,25, by = 5), limits = c(0,25)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.6, by = 0.05), limits = c(0,0.6)) +
  # Artificial axes
  geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
  geom_vline(aes(xintercept = 0), color = 'black', linewidth = .5) +
  # Title and labels
  ggtitle('A') +
  labs(y = c(expression(paste(italic(' a')['phy'],' (665 nm) (m'^'-1',')'))), 
       x = (expression(paste('[Chl',~italic('a'),'] (mg.m'^'-3',')')))) +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

plot_aph_vs_chla_field_data_full

#ggsave('aph_vs_chla_field_data3_full.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

## graph inversed aphy vs Chl alternative: figure
plot_aph_vs_chla_field_data_alt_fig <- ggplot(data = field_data, aes(x = Chlorophyll_a, y = aphy_665)) +
  # points of 8 stations with aphy_modelled
  geom_point(color = 'midnightblue', shape = 1, size = 3, stroke = .75) +
  # 2 plcors for measured aphy represented with points
  # geom_point(aes(x = Chlorophyll_a, y = aphy_665.measured), color = 'red', shape = 0, size = 3, stroke = .75) +
  # geom_point(aes(x = Chlorophyll_a, y = aphy_665_homemade), color = 'darkred', shape = 2, size = 3, stroke = .75) +
  # 2 aphy_star from measured aphy
  # geom_abline(aes(intercept=0, slope=lm_aphy_measured_vs_chla_0$coefficients), color ='red', linewidth = .5, linetype = 2) +
  # geom_abline(aes(intercept=0, slope=lm_aphy_homemade_vs_chla_0$coefficients), color ='darkred', linewidth = .5, linetype = 2) +
  # aphy_star from Gons et al. (2002) and a_sol_chla from Bricaud et al. (2004)
  geom_abline(aes(intercept=0, slope=0.0146), color ='gray25', linewidth = .6, linetype = 1) +
  geom_abline(aes(intercept=0, slope=0.021656), color ='black', linewidth = .6, linetype = 1) +
  geom_abline(aes(intercept=0, slope=lm_aphy_modelled_vs_chla_0$coefficients), color ='midnightblue', linewidth = .5, linetype = 2) +
  geom_abline(aes(intercept=0, slope=0.0073), color ='dodgerblue3', linewidth = .5, linetype = 2) +
  geom_abline(aes(intercept=0, slope=0.0111), color ='dodgerblue3', linewidth = .5, linetype = 2) +
  
  # Format of x axis
  scale_x_continuous(breaks = seq(0,25, by = 5), limits = c(0,25)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.4, by = 0.05), limits = c(0,0.4)) +
  # Artificial axes
  geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
  geom_vline(aes(xintercept = 0), color = 'black', linewidth = .5) +
  # Title and labels
  ggtitle('A') +
  labs(y = c(expression(paste(italic(' a')['phy'],' (665 nm) (m'^'-1',')'))), 
       x = (expression(paste('[Chl',~italic('a'),'] (mg.m'^'-3',')')))) +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

plot_aph_vs_chla_field_data_alt_fig

# ggsave('aph_vs_chla_field_data3_alt_supp.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

## graph inversed aphy vs Chl alternative: supplementary figure
plot_aph_vs_chla_field_data_alt_supp <- ggplot(data = field_data, aes(x = Chlorophyll_a, y = aphy_665)) +
  # points of 8 stations with aphy_modelled
  geom_point(color = 'midnightblue', shape = 1, size = 3, stroke = .75) +
  # 2 plcors for measured aphy represented with points
  geom_point(aes(x = Chlorophyll_a, y = aphy_665.measured), color = 'red', shape = 0, size = 3, stroke = .75) +
  geom_point(aes(x = Chlorophyll_a, y = aphy_665_homemade), color = 'darkred', shape = 2, size = 3, stroke = .75) +
  # 2 aphy_star from measured aphy
  geom_abline(aes(intercept=0, slope=lm_aphy_measured_vs_chla_0$coefficients), color ='red', linewidth = .5, linetype = 2) +
  geom_abline(aes(intercept=0, slope=lm_aphy_homemade_vs_chla_0$coefficients), color ='darkred', linewidth = .5, linetype = 2) +
  # aphy_star from Gons et al. (2002) and a_sol_chla from Bricaud et al. (2004)
  geom_abline(aes(intercept=0, slope=0.0146), color ='gray25', linewidth = .6, linetype = 1) +
  geom_abline(aes(intercept=0, slope=0.021656), color ='black', linewidth = .6, linetype = 1) +
  geom_abline(aes(intercept=0, slope=lm_aphy_modelled_vs_chla_0$coefficients), color ='midnightblue', linewidth = .5, linetype = 2) +
  # geom_abline(aes(intercept=0, slope=0.0073), color ='dodgerblue3', linewidth = .5, linetype = 2) +
  # geom_abline(aes(intercept=0, slope=0.0111), color ='dodgerblue3', linewidth = .5, linetype = 2) +
  
  # Format of x axis
  scale_x_continuous(breaks = seq(0,25, by = 5), limits = c(0,25)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.6, by = 0.05), limits = c(0,0.6)) +
  # Artificial axes
  geom_hline(aes(yintercept = 0), color = 'black', linewidth = .5) +
  geom_vline(aes(xintercept = 0), color = 'black', linewidth = .5) +
  # Title and labels
  ggtitle('B') +
  labs(y = c(expression(paste(italic(' a')['phy'],' (665 nm) (m'^'-1',')'))), 
       x = (expression(paste('[Chl',~italic('a'),'] (mg.m'^'-3',')')))) +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=10), axis.title.y =element_text(size=10), axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

plot_aph_vs_chla_field_data_alt_supp

# ggsave('aph_vs_chla_field_data3_alt_supp.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

#### Patchwork plots for optics ####

### Patchwork of the 4 plots - full
ggarrange(plot_aph_vs_chla_field_data_full, plot_eval_aphy_star, ncol = 1, common.legend = FALSE)

# Saving the plot
# ggsave('plot_optics_remote_sensing3_full.tiff', dpi = 600, width = 80, height = 170, units = 'mm')

### Patchwork of the 4 plots - alternative
ggarrange(plot_aph_vs_chla_field_data_alt_fig, plot_eval_aphy_star, ncol = 1, common.legend = FALSE)

# Saving the plot
# ggsave('plot_optics_remote_sensing3_alt_fig.tiff', dpi = 600, width = 80, height = 170, units = 'mm')

### Supplementary figure for plcor
ggarrange(plot_aph_cuvette_vs_filter, plot_aph_vs_chla_field_data_alt_supp, ncol = 2, common.legend = FALSE)

# Saving the plot
# ggsave('plot_optics_plcor_alt_supp.tiff', dpi = 600, width = 180, height = 95, units = 'mm')
 
#### Creating a legend for chl a maps ####

value <- as.data.frame(seq(0,1, by = 0.01))
chla <- as.data.frame(seq(0,100, by = 1))
legend <- bind_cols(chla, value)
colnames(legend) <- c('chla','value')
legend <- mutate(legend,var = 'chla' )

ggplot(data = legend, aes(chla, var)) +
  geom_tile(aes(x = chla, fill = value), color = 'black', linewidth = .2) +
  scale_fill_gradientn(colors = cmocean(name = 'algae')(100)) +
  # Title and axis labels
  # ggtitle('Legend:') +
  labs(y = NULL, 
       x = (expression(paste('Chl ', italic('a'),' concentration (mg.m'^'-3',')')))) +
  # Theme
  theme(plot.title = element_text(size = 9, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=9, color = 'black'),
        axis.text.x = element_text(size=7, color = 'black'),
        axis.title.y = element_text(color = 'white'),
        axis.text.y = element_text(color = 'white'),
        # Legend
        legend.text = NULL, 
        legend.position = NULL,
        legend.key = NULL,
        # Panel
        panel.grid.major.x = element_line(color = 'black', linewidth = .2), panel.grid.minor = NULL,
        panel.grid.major.y = NULL,
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(fill = 'none')

# ggsave('legend_chla_maps.tiff', dpi = 600, width = 178, height = 25, units = 'mm')
