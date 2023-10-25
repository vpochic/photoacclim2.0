### This script is part of the available material associated with the article 
### "Photoacclimation in the kleptoplastidic ciliate Mesodinium rubrum and its 
### cryptophyte prey Teleaulax amphioxeia: phenotypic variability and implications
### for red tide remote sensing", submitted to the Journal of Plankton Research on 2023/09/29

### This script and all the associated data files are made available publicly to follow the
### FAIR (Findable, Accessible, Interoperable, Reusable) principles for scientific data.

### This script was written entirely by the authors of the article. Please give credit if reusing all or parts
### of the script

####### Script Photoacclim 2.0: Pigment concentrations ###
### Date of last modification: 2023/10/25
# Author: Victor Pochic, ISOMer, Nantes Universite

# This scripts processes and plots data of pigment concentrations (HPLC + hydrosoluble pigments) 
# for the photoacclimation experiment presented in the article.

### Packages needed:

library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)

#### Import data ####

# For data obtained with the TECAN spectrophotometer, the data need to be first transformed into a clean
# .csv file

ODPE <- read.csv('OD/ODPE_photoacclim2.csv', sep = ';', header = TRUE, dec = ',') %>%
  mutate_at(vars(-'Wavel.'), list(~ .-min(.))) %>% # Putting all spectra at the same base level, with min equal to 0
  mutate_at(vars(-c('Wavel.')), list(~ .-((C1+C2+C3)/3))) # Subtracting the blank

# For the data from the experiment at lower temperature (Meso only)
ODPE2 <- read.csv('OD/ODPE_temp17.5.csv', sep = ';', header = TRUE, dec = ',') %>%
  mutate_at(vars(-'Wavel.'), list(~ .-min(.))) %>% # Putting all spectra at the same base level, with min equal to 0
  mutate_at(vars(-c('Wavel.')), list(~ .-((F7+F8+F9)/3))) # Subtracting the blank

# Dataframe with HPLC pigments concentrations and cell counts
HPLCdata <- read.csv('HPLC_photoacclim/HPLCdata_photoacclim2.csv', sep = ';', header = TRUE, dec = ',')
HPLCdata2 <- read.csv('HPLC_photoacclim/HPLCdata_temp17.5.csv', sep = ';', header = TRUE, dec = ',')

# Dataframe with info about cell_counts, biovolumes etc...
Samples_photoacclim2 <- read.csv('Experiment_samples/Samples_photoacclim2.csv', sep = ';', dec = ',', header = TRUE)

# For the exp temp = 17.5, another dataframe
Samples_temp_17.5 <- read.csv2('Experiment_samples/Samples_photoacclim2_temp17.5.csv', header = TRUE)


#### Calculations of pigment concentrations ####

# Parameters for computing the PE concentration
p_surn.Tele = 0.148 # fraction of cells lost during centrifugation /!\ For Teleaulax only /!\
# For Meso, p_surn = 0
# p_surn was calculated using flow cytometer data, in the script 'Cyto_growth_rates'
epsilon = 12.6 # molar extinction coefficient of PE 545 in L.g-1.cm-1 (MacColl et al. 1976)
d = 0.6 # the optical path in the spectrophotometer in cm
Vb = 5 # the volume of buffer used for extraction, in mL
Vs = 40 # the volume of sample, in mL
MMPE = 57030 # the molar mass of PE 545, in g.mol-1

# Parameters that changed for the experiment at temp = 17.5
Vb2 = 4 # the volume of buffer used for extraction, in mL
Vs2 = 10 # the volume of sample, in mL

# Exp at temp = 21
PE_conc <- ODPE %>%
  filter(Wavel. == 552) %>%
  mutate_at(vars(-'Wavel.'), list(~ (.)/(epsilon*d)*(Vb/Vs)*1000)) %>% # Calculating the PE concentration in mg.L,
  # and only applying the correction for percentage of lost cells for Teleaulax later
  select(-c('Wavel.', 'C1', 'C2', 'C3')) %>% # Removing the wavelength and the blanks
  # We now have the PE concentration in mg.L
  melt(variable.name = 'Sample', value.name = 'PE') %>%
  mutate(PE = ifelse(Sample %in% c('B1', 'B2', 'B3'), PE*(8/7.5), PE)) %>% # Applying a correction
  # factor for Meso_20, because Vb and Vs were different
  left_join(HPLCdata, by = 'Sample') %>%
  select(c('Sample', 'Name', 'Species', 'Irradiance', 'Cell_count', 'Unit_cell_count', 
           'Unit_conc_pigments', 'PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 
           'Fuco', 'Diadino', 'Diato', 'Croco', 'Monado', 'Phaeophy', 'Phaeophor')) %>% # Selecting only necessary variables
  mutate(PE = ifelse(Species == 'Tele', PE+PE*p_surn.Tele, PE)) %>% # Correcting for cells lost in centri for Tele only
  mutate_at(vars(-c('Sample', 'Name', 'Species', 'Irradiance', 'Cell_count','Unit_cell_count', 'Unit_conc_pigments')),
            list(~ (./(Cell_count*1000))*10^12)) %>% # Dividing by the cell count to obtain cellular concentrations of pigments. 
  # *1000 because cell count is in cells.mL, *10^12 to obtain the result in fg per cell
  mutate(Unit_conc_pigments = 'fg.cell')

# Exp at temp = 17.5
PE_conc2 <- ODPE2 %>%
  filter(Wavel. == 552) %>%
  mutate_at(vars(-'Wavel.'), list(~ (.)/(epsilon*d)*(Vb2/Vs2)*1000)) %>% # Calculating the PE concentration in mg.L
  select(-c('Wavel.', 'F7', 'F8', 'F9')) %>% # Removing the wavelength and the blanks
  # We now have the PE concentration in mg.L
  melt(variable.name = 'Sample', value.name = 'PE') %>%
  left_join(HPLCdata2, by = 'Sample') %>%
  select(c('Sample', 'Name', 'Species', 'Irradiance', 'Cell_count', 'Unit_cell_count', 
           'Unit_conc_pigments', 'PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 
           'Fuco', 'Diadino', 'Diato', 'Croco', 'Monado', 'Phaeophy', 'Phaeophor')) %>% # Selecting only necessary variables
  mutate_at(vars(-c('Sample', 'Name', 'Species', 'Irradiance', 'Cell_count','Unit_cell_count', 'Unit_conc_pigments')),
            list(~ (./(Cell_count*1000))*10^12)) %>% # Dividing by the cell count to obtain cellular concentrations of pigments. 
  # *1000 because cell count is in cells.mL, *10^12 to obtain the result in fg per cell
  mutate(Unit_conc_pigments = 'fg.cell')

# From there, we calculate the pigment concentrations per biovolume
Pigments_biov <- PE_conc %>%
  left_join(Samples_photoacclim2, by = c('Name', 'Irradiance', 'Species')) %>%
  mutate_at(vars('PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 'Croco',
                 'Monado', 'Fuco', 'Diadino', 'Diato', 'Phaeophor', 'Phaeophy'), list(~ (./Biovolume))) %>%
  mutate(Unit_conc_pigments = 'fg.microm3') %>%
  mutate(Vs_CHN = as.character(Vs_CHN)) # Needed to combine with the Nas in the 2nd dataframe

Pigments_biov2 <- PE_conc2 %>%
  left_join(Samples_temp_17.5, by = c('Name', 'Irradiance', 'Species')) %>%
  mutate_at(vars('PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 'Croco',
                 'Monado', 'Fuco', 'Diadino', 'Diato', 'Phaeophor', 'Phaeophy'), list(~ (./Biovolume))) %>%
  mutate(Unit_conc_pigments = 'fg.microm3')

# Joining the tables of the 2 experiments
Pigments_biov <- bind_rows(Pigments_biov, Pigments_biov2)

# Saving the results as .csv files
# write.csv2(PE_conc, file = 'Pigments_per_cell_photoacclim2.csv', row.names = FALSE)
# write.csv2(PE_conc2, file = 'Pigments_per_cell_photoacclim.csv', row.names = FALSE)
# write.csv2(Pigments_biov, file = 'Pigments_per_biov_photoacclim.csv', row.names = FALSE)

# We can now calculate the ratios of [Pigments] / [Chl a]
Pigment_ratios <- Pigments_biov %>%
  mutate_at(vars('PE', 'Chl.c2', 'Allox', 'alpha.car', 'Croco',
                 'Monado', 'Fuco', 'Diadino', 'Diato', 'Phaeophor', 
                 'Phaeophy', 'Chl.a'), list(~ ./Chl.a)) %>% # Note that Chl a is last, else it f*cks up
  mutate(Unit_conc_pigments = 'g.g[Chl a]-1')

# Saving the results as a .csv file
# write.csv2(Pigment_ratios, file = 'Pigment_ratios_photoacclim.csv', row.names = FALSE) 

### Grouping the data by condition

# Exclude replicates considered as outliers
Pigments_biov <- filter(Pigments_biov, !(Name %in% c('Meso20_1','Meso65_3'))) # Removing outliers in the Meso20 and Meso65 conditions
Pigment_ratios <- filter(Pigment_ratios, !(Name %in% c('Meso20_1','Meso65_3'))) # Removing outliers in the Meso20 and Meso65 conditions

# Pigment concentrations

Pigments_means <- Pigments_biov %>%
  group_by(Species, Temp, Irradiance) %>%
  summarise_at(vars('PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 'Croco',
                    'Monado', 'Fuco', 'Diadino', 'Diato', 'Phaeophor', 'Phaeophy'), list( ~ mean(.)), .groups = 'keep')

Pigments_sds <- Pigments_biov %>%
  group_by(Species, Temp, Irradiance) %>%
  summarise_at(vars('PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 'Croco',
                    'Monado', 'Fuco', 'Diadino', 'Diato', 'Phaeophor', 'Phaeophy'), list( ~ sd(.)), .groups = 'keep')

Pigments <- left_join(Pigments_means, Pigments_sds, by = c('Species', 'Irradiance', 'Temp'), 
                      suffix = c('_mean', '_sd'))
# Pigment ratios

Ratios_means <- Pigment_ratios %>%
  group_by(Species, Temp, Irradiance) %>%
  summarise_at(vars('PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 'Fuco', 'Diadino', 'Phaeophy'), list( ~ mean(.)), .groups = 'keep')

Ratios_sds <- Pigment_ratios %>%
  group_by(Species, Temp, Irradiance) %>%
  summarise_at(vars('PE', 'Chl.a', 'Chl.c2', 'Allox', 'alpha.car', 'Fuco', 'Diadino', 'Phaeophy'), list( ~ sd(.)), .groups = 'keep')

Ratios <- left_join(Ratios_means, Ratios_sds, by = c('Species', 'Irradiance', 'Temp'), 
                    suffix = c('_mean', '_sd'))

#write.csv2(Ratios, 'Ratios_pigments_lab3.csv', row.names = FALSE)

#### Plots : pigment concentrations ####

# PE
PE.plot<-ggplot(Pigments, aes(x = Irradiance, y = PE_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2)) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = PE, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,75, by = 25), limits = c(0,95)) +
  # Title and labels (legend and axis)
  ggtitle("PE") +
  labs(y = c(expression(paste('[PE 545] (fg.',mu,'m'^'-3',')'))), 
       x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')'))),
       color = 'Temperature (°C)',
       shape = 'Species') +
  # Errorbars
  geom_errorbar(aes(ymin=PE_mean-PE_sd, ymax=PE_mean+PE_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

PE.plot

# Chl a
Chla.plot<-ggplot(Pigments, aes(x = Irradiance, y = Chl.a_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2)) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Chl.a, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,15, by = 5), limits = c(0,17.5)) +
  # Title and labels (legend and axis)
  ggtitle(c(expression(paste("Chl ", italic('a'))))) +
  labs(y = c(expression(paste('[Chl ', italic('a'), '] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = 'Temperature (°C)',
       shape = 'Species') +
  # Error bars
  geom_errorbar(aes(ymin=Chl.a_mean-Chl.a_sd, ymax=Chl.a_mean+Chl.a_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Chla.plot

# Chl c2
Chlc2.plot<-ggplot(Pigments, aes(x = Irradiance, y = Chl.c2_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Chl.c2, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Chl ', italic('c')['2'])))) +
  labs(y = c(expression(paste('[Chl ', italic('c')['2'],'] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Chl.c2_mean-Chl.c2_sd, ymax=Chl.c2_mean+Chl.c2_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Chlc2.plot

# Allox
Allox.plot<-ggplot(Pigments, aes(x = Irradiance, y = Allox_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2)) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Allox, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,5, by = 1), limits = c(0,5)) +
  # Title and labels (axis and legend)
  ggtitle('Allox') +
  labs(y = c(expression(paste('[Allox] (fg.',mu,'m'^'-3',')'))), 
       x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')'))),
       color = 'Temperature (°C)',
       shape = 'Species') +
  geom_errorbar(aes(ymin=Allox_mean-Allox_sd, ymax=Allox_mean+Allox_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Allox.plot


# Patchwork plot
ggarrange(Chla.plot, Chlc2.plot, Allox.plot, PE.plot, common.legend = TRUE, legend ='bottom')

# Saving the plot
# ggsave('Pigment_conc_v3.tiff', dpi = 600, width = 164, height = 157.6, units = 'mm')

#### Supp plots : pigment ratios ####

# PE
PE_ratio.plot<-ggplot(Ratios, aes(x = Irradiance, y = PE_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2)) +
  geom_point(size = 4) +
  geom_point(data = Pigment_ratios, aes(x = Irradiance, y = PE, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  scale_y_continuous(breaks = seq(2,5, by = 1), limits = c(1.5,6)) +
  # Format of y axis
  # scale_y_continuous(breaks = seq(0,75, by = 25), limits = c(0,95)) +
  # Title and labels (legend and axis)
  ggtitle(c(expression(paste('PE/Chl ',italic('a'))))) +
  labs(y = c(expression(paste('[PE 545]/[Chl ',italic('a'),'] (g.g'^'-1',')'))), 
       x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')'))),
       color = 'Temperature (°C)',
       shape = 'Species') +
  # Errorbars
  geom_errorbar(aes(ymin=PE_mean-PE_sd, ymax=PE_mean+PE_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

PE_ratio.plot

# Chl c2 / Chl a !!! PROBLEM
Chlc2_ratio.plot<-ggplot(Ratios, aes(x = Irradiance, y = Chl.c2_mean, color = factor(Temp), shape = Species)) +
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = Pigment_ratios, aes(x = Irradiance, y = Chl.c2, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  scale_y_continuous(breaks = seq(0.1,0.2, by = 0.025), limits = c(0.1,0.2)) +
  # Title and labels
  ggtitle(c(expression(paste('Chl ', italic('c')['2'],'/Chl ',italic('a'))))) +
  labs(y = c(expression(paste('[Chl ', italic('c')['2'],']/[Chl ',italic('a'),'] (g.g'^'-1',')'))), x = NULL) +
  # Error bars
  geom_errorbar(aes(ymin=Chl.c2_mean-Chl.c2_sd, ymax=Chl.c2_mean+Chl.c2_sd), width = 3,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.5, 0.3), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Chlc2_ratio.plot

# Allox / Chl a
Allox_ratio.plot<-ggplot(Ratios, aes(x = Irradiance, y = Allox_mean, color = factor(Temp), shape = Species)) +
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  geom_point(size = 4) +
  geom_point(data = Pigment_ratios, aes(x = Irradiance, y = Allox, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Title and labels
  ggtitle(c(expression(paste('Allox/Chl ',italic('a'))))) +
  labs(y = c(expression(paste('[Allox]/[Chl ',italic('a'),'] (g.g'^'-1',')'))), x = NULL) +
  # Format of axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # scale_y_continuous(breaks = seq(0.2,0.4, by = 0.05), limits = c(0.2,0.41)) +
  # Error bars
  geom_errorbar(aes(ymin=Allox_mean-Allox_sd, ymax=Allox_mean+Allox_sd), width = 3,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

Allox_ratio.plot

# Patchwork plot
ggarrange(Chlc2_ratio.plot, Allox_ratio.plot, PE_ratio.plot, align = 'v', ncol = 1, common.legend = TRUE, legend ='bottom')

# Saving the plot
# ggsave('Pigment_ratios2.tiff', dpi = 600, width = 80, height = 173.6, units = 'mm')
#### Supp plots : minor pigments ####

### Plots of the "minor pigments" for supp mat

# Beta,epsilon-carotene
alpha.car.plot<-ggplot(Pigments, aes(x = Irradiance, y = alpha.car_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2)) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = alpha.car, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.6, by = 0.2), limits = c(0,0.65)) +
  # Title and labels (legend and axis)
  ggtitle(c(expression(paste(beta,',',epsilon,'-carotene')))) +
  labs(y = c(expression(paste('[',beta,',',epsilon,'-carotene] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = 'Temperature (°C)',
       shape = 'Species') +
  # Error bars
  geom_errorbar(aes(ymin=alpha.car_mean-alpha.car_sd, ymax=alpha.car_mean+alpha.car_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

alpha.car.plot

# Crocoxanthin
Croco.plot<-ggplot(Pigments, aes(x = Irradiance, y = Croco_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Croco, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,0.6, by = 0.2), limits = c(0,0.65)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Crocoxanthin')))) +
  labs(y = c(expression(paste('[Crocoxanthin] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Croco_mean-Croco_sd, ymax=Croco_mean+Croco_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Croco.plot

# Monadoxanthin
Monado.plot<-ggplot(Pigments, aes(x = Irradiance, y = Monado_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Monado, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Monadoxanthin')))) +
  labs(y = c(expression(paste('[Monadoxanthin] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Monado_mean-Monado_sd, ymax=Monado_mean+Monado_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Monado.plot

# Fucoxanthin
Fuco.plot<-ggplot(Pigments, aes(x = Irradiance, y = Fuco_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Fuco, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Fucoxanthin*')))) +
  labs(y = c(expression(paste('[Fucoxanthin] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Fuco_mean-Fuco_sd, ymax=Fuco_mean+Fuco_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Fuco.plot

# Diatoxanthin
Diato.plot<-ggplot(Pigments, aes(x = Irradiance, y = Diato_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Diato, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Diatoxanthin*')))) +
  labs(y = c(expression(paste('[Diatoxanthin] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Diato_mean-Diato_sd, ymax=Diato_mean+Diato_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Diato.plot

# Diadinoxanthin
Diadino.plot<-ggplot(Pigments, aes(x = Irradiance, y = Diadino_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Diadino, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Diadinoxanthin*')))) +
  labs(y = c(expression(paste('[Diadinoxanthin] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Diadino_mean-Diadino_sd, ymax=Diadino_mean+Diadino_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Diadino.plot

# Phaeophytin
Phaeophy.plot<-ggplot(Pigments, aes(x = Irradiance, y = Phaeophy_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Phaeophy, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Phaeophytin ', italic('a'))))) +
  labs(y = c(expression(paste('[Phaeophytin ', italic('a'),'] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Phaeophy_mean-Phaeophy_sd, ymax=Phaeophy_mean+Phaeophy_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Phaeophy.plot

# Phaeophorbide a
Phaeophor.plot<-ggplot(Pigments, aes(x = Irradiance, y = Phaeophor_mean, color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4','firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Pigments_biov, aes(x = Irradiance, y = Phaeophor, color = factor(Temp), shape = Species), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,220, by = 50), limits = c(0,225)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(0,3, by = 1), limits = c(0,3.25)) +
  # Title and labels (axis and legend)
  ggtitle(c(expression(paste('Phaeophorbide ', italic('a'))))) +
  labs(y = c(expression(paste('[Phaeophorbide ',italic('a'),'] (fg.',mu,'m'^'-3',')'))), 
       x = NULL,
       color = c(expression(paste('Temperature (°C):'))))+
  geom_errorbar(aes(ymin=Phaeophor_mean-Phaeophor_sd, ymax=Phaeophor_mean+Phaeophor_sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1)) +
  guides(shape = guide_legend('Species:'))

Phaeophor.plot

# Patchwork plot
ggarrange(alpha.car.plot, Croco.plot, Monado.plot, Fuco.plot, Diato.plot,
          Diadino.plot, Phaeophy.plot, Phaeophor.plot, align = 'v', nrow = 3, ncol = 3, common.legend = TRUE, legend ='bottom')

# Saving the plot
# ggsave('Pigment_minors_supp.tiff', dpi = 600, width = 270, height = 180, units = 'mm')
