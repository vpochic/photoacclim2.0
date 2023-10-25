### This script is part of the available material associated with the article 
### "Photoacclimation in the kleptoplastidic ciliate Mesodinium rubrum and its 
### cryptophyte prey Teleaulax amphioxeia: phenotypic variability and implications
### for red tide remote sensing", submitted to the Journal of Plankton Research on 2023/09/29

### This script and all the associated data files are made available publicly to follow the
### FAIR (Findable, Accessible, Interoperable, Reusable) principles for scientific data.

### This script was written entirely by the authors of the article. Please give credit if reusing all or parts
### of the script

####### Script Photoacclim 2.0: Flow cytometry and growth rates ###
### Date of last modification: 2023/10/25
# Author: Victor Pochic, ISOMer, Nantes Universite

# This scripts processes and plots flow cytometry data for the photoacclimation 
# experiment presented in the article, and calculates the associated growth rates

### Packages needed:

library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)

#### Import data ####

Tele1 <- read.csv('Cyto_data_photoacclim/Tele1.csv', sep = ';', dec = ',', header = TRUE)
Tele2 <- read.csv('Cyto_data_photoacclim/Tele2.csv', sep = ';', dec = ',', header = TRUE)
Tele3 <- read.csv('Cyto_data_photoacclim/Tele3.csv', sep = ';', dec = ',', header = TRUE)
Tele4 <- read.csv('Cyto_data_photoacclim/Tele4.csv', sep = ';', dec = ',', header = TRUE)
Meso1 <- read.csv('Cyto_data_photoacclim/Meso1.csv', sep = ';', dec = ',', header = TRUE)
Meso2 <- read.csv('Cyto_data_photoacclim/Meso2.csv', sep = ';', dec = ',', header = TRUE)
Meso3 <- read.csv('Cyto_data_photoacclim/Meso3.csv', sep = ';', dec = ',', header = TRUE)
Meso4 <- read.csv('Cyto_data_photoacclim/Meso4.csv', sep = ';', dec = ',', header = TRUE)
Meso5 <- read.csv('Cyto_data_photoacclim/Meso5.csv', sep = ';', dec = ',', header = TRUE)

#### Curate data ####

Cyto_DB <- bind_rows(Tele1, Tele2) %>%
  bind_rows(Tele3) %>%
  bind_rows(Tele4) %>%
  bind_rows(Meso1) %>%
  bind_rows(Meso2) %>%
  bind_rows(Meso3) %>%
  bind_rows(Meso4) %>%
  bind_rows(Meso5) %>%
  mutate(Species = ifelse(grepl('Tele', SID, fixed = TRUE) == TRUE, 'Tele', 'Meso')) %>%
  # Adds the information of the species based on the sample ID
  mutate(Irradiance = ifelse(grepl('200', SID, fixed = TRUE) == TRUE, '200', 
                             ifelse(grepl('80', SID, fixed = TRUE) == TRUE, '80', '20'))) %>%
  # Adds the information of the irradiance based on the sample ID
  mutate(Sample = ifelse(grepl('surn', SID, fixed = TRUE) == TRUE, 'surn', 'culture')) %>%
  # Adds the information of the replicate number based on the sample ID
  mutate(Replicate = ifelse(grepl('_1', SID, fixed = TRUE) == TRUE, 1, 
                            ifelse(grepl('_2', SID, fixed = TRUE) == TRUE, 2, 3))) %>%
  # Adds the info of whether the sample is culture or supernatant
  mutate(PR = ifelse(grepl('PR', SID, fixed = TRUE) == TRUE, 'Yes', 'No')) %>%
  # Adds the info whether the count is done after transplanting of the culture
  mutate(Date = ymd(Date)) %>%
  # modifies the date format
  mutate(Day = day(Date)) %>%
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) # We cut the pipe here so we have a reference that we can call for the 
# next step (calculating the experimental time)

Cyto_data <- Cyto_DB %>%
  mutate(Time_exp = Day - min(Cyto_DB$Day)) %>%
  # Computes the experimental time in days(0 = start of the experiment)
  mutate(Dilution_factor = ifelse((grepl('200', SID, fixed = TRUE) == TRUE |
                                     grepl('80', SID, fixed = TRUE) == TRUE) &
                                    (Time_exp == 3 | Time_exp == 5 | Time_exp == 8) &
                                    (PR == 'No') & (Species == 'Tele'), 2, 1)) %>%
  # There was a dilution of some samples before counting, at certain time points
  mutate(Count.mL = Count.mL*Dilution_factor)
# Taking into account the dilution factors for the cell counts

Cyto_groups <- Cyto_data %>%
  group_by(Species, Irradiance, Time_exp, PR) %>%
  filter(Sample == 'culture') %>% # Excluding the supernatant samples
  filter(grepl('20N', SID, fixed = TRUE) == FALSE) # Excluding the samples of Mesodinium 'clouds'

Cyto_means <- Cyto_groups %>%
  summarise_at(vars('Count.mL', 'FSC.H.Median', 'SSC.H.Median', 'B2.H.Median', 'B3.H.Median'),
               list( ~ mean(.)), .groups = 'keep')
Cyto_sds <- Cyto_groups %>%
  summarise_at(vars('Count.mL', 'FSC.H.Median', 'SSC.H.Median', 'B2.H.Median', 'B3.H.Median'),
               list( ~ sd(.)), .groups = 'keep')

Cyto_stats <- left_join(Cyto_means, Cyto_sds, by = c('Irradiance', 'Species', 'Time_exp', 'PR'), 
                        suffix = c('_mean', '_sd'))

# Creating separated data frames for each species
# Tele
Cyto_stats_Tele <- Cyto_stats %>%
  filter(Species == 'Tele')

Cyto_groups_Tele <- Cyto_groups %>%
  filter(Species == 'Tele')

# Meso
Cyto_stats_Meso <- Cyto_stats %>%
  filter(Species == 'Meso')

Cyto_groups_Meso <- Cyto_groups %>%
  filter(Species == 'Meso')

#### Plotting ####

#### Cell counts ####

### Teleaulax

Cyto_count_Tele.plot <- ggplot(Cyto_stats_Tele, aes(x = Time_exp, y = Count.mL_mean, color = Irradiance)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3', 'darkorange')) +
  geom_point(size = 4, shape = 17) +
  geom_point(data = Cyto_groups_Tele, aes(x = Time_exp, y = Count.mL, color = Irradiance), 
             size = 2, shape = 17) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,12, by = 2)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Cells per mL - ', italic('T. amphioxeia'))))) +
  labs(y = c(expression(paste('Cells.mL'^'-1'))), x = 'Time (days)') +
  # Error bars
  geom_errorbar(aes(ymin=Count.mL_mean-Count.mL_sd, ymax=Count.mL_mean+Count.mL_sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.872, 0.85), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

Cyto_count_Tele.plot

# Saving the plot as image
# ggsave('Cyto_count_Tele.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

### Mesodinium

Cyto_count_Meso.plot <- ggplot(Cyto_stats_Meso, aes(x = Time_exp, y = Count.mL_mean, color = Irradiance)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3', 'darkorange')) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = Cyto_groups_Meso, aes(x = Time_exp, y = Count.mL, color = Irradiance), 
             size = 2, shape = 16) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,13, by = 2)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,30000, by = 10000), limits = c(0,30000)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Cells per mL - ', italic('M. rubrum'))))) +
  labs(y = c(expression(paste('Cells.mL'^'-1'))), x = 'Time (days)') +
  # Error bars
  geom_errorbar(aes(ymin=Count.mL_mean-Count.mL_sd, ymax=Count.mL_mean+Count.mL_sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.19, 0.79), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

Cyto_count_Meso.plot

# Saving the plot as image
# ggsave('Cyto_count_Meso.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

#### B2H - Fluo PE 545 ####

### Tele 

FluoPE_Tele.plot <- ggplot(Cyto_stats_Tele, aes(x = Time_exp, y = B2.H.Median_mean, color = Irradiance)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3', 'darkorange')) +
  geom_point(size = 4, shape = 17) +
  geom_point(data = Cyto_groups_Tele, aes(x = Time_exp, y = B2.H.Median, color = Irradiance), 
             size = 2, shape = 17) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,12, by = 2)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,30, by = 10), limits = c(0,30)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Fluorescence of PE 545 - ', italic('T. amphioxeia'))))) +
  labs(y = 'B2-H', x = 'Time (days)') +
  # Error bars
  geom_errorbar(aes(ymin=B2.H.Median_mean-B2.H.Median_sd, ymax=B2.H.Median_mean+B2.H.Median_sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.8, 0.60), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

FluoPE_Tele.plot

# Saving the plot as image
# ggsave('FluoPE_Tele.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

### Meso

FluoPE_Meso.plot <- ggplot(Cyto_stats_Meso, aes(x = Time_exp, y = B2.H.Median_mean, color = Irradiance)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3', 'darkorange')) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = Cyto_groups_Meso, aes(x = Time_exp, y = B2.H.Median, color = Irradiance), 
             size = 2, shape = 16) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,13, by = 2)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,20, by = 5), limits = c(0,20)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Fluorescence of PE 545 - ', italic('M. rubrum'))))) +
  labs(y = 'B2-H', x = 'Time (days)') +
  # Error bars
  geom_errorbar(aes(ymin=B2.H.Median_mean-B2.H.Median_sd, ymax=B2.H.Median_mean+B2.H.Median_sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.75, 0.2), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

FluoPE_Meso.plot

# Saving the plot as image
# ggsave('FluoPE_Meso.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')


#### B3H - Fluo Chl a ####

### Tele

FluoChla_Tele.plot <- ggplot(Cyto_stats_Tele, aes(x = Time_exp, y = B3.H.Median_mean, color = Irradiance, shape = Species)) +
  #Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3', 'darkorange')) +
  geom_point(size = 4, shape = 17) +
  geom_point(data = Cyto_groups_Tele, aes(x = Time_exp, y = B3.H.Median, color = Irradiance), 
             size = 2, shape = 17) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,12, by = 2)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(15,45, by = 15), limits = c(10,50)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Fluorescence of Chl ', italic('a'), ' - ', italic('T. amphioxeia'))))) +
  labs(y = 'B3-H', x = 'Time (days)') +
  geom_errorbar(aes(ymin=B3.H.Median_mean-B3.H.Median_sd, ymax=B3.H.Median_mean+B3.H.Median_sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.2, 0.2), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

FluoChla_Tele.plot

# Saving the plot as image
# ggsave('FluoChla_Tele.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

### Meso

FluoChla_Meso.plot <- ggplot(Cyto_stats_Meso, aes(x = Time_exp, y = B3.H.Median_mean, color = Irradiance, shape = Species)) +
  #Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3', 'darkorange')) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = Cyto_groups_Meso, aes(x = Time_exp, y = B3.H.Median, color = Irradiance), 
             size = 2, shape = 16) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,13, by = 2)) +
  # Format of y axis
  #scale_y_continuous(breaks = seq(15,45, by = 15), limits = c(10,50)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Fluorescence of Chl ', italic('a'), ' - ', italic('M. rubrum'))))) +
  labs(y = 'B3-H', x = 'Time (days)') +
  geom_errorbar(aes(ymin=B3.H.Median_mean-B3.H.Median_sd, ymax=B3.H.Median_mean+B3.H.Median_sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.2, 0.2), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', size = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

FluoChla_Meso.plot

# Saving the plot as image
# ggsave('FluoChla_Meso.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

#### Calculating p_surn for the PE concentrations ####

p_surn <- Cyto_data %>%
  filter(Time_exp == 12) %>%
  filter(grepl('_2', SID, fixed = TRUE)&(Count != 0)&(Species == 'Tele')) %>%
  # We don't check the percentage of cells in the supernatant for Mesodinium, 
  # because it is assumed that no cells remain in the supernatant after centrifugation for this species
  group_by(Irradiance) %>%
  arrange(as.numeric(Irradiance))

p_surn_Tele <- as.data.frame(as.numeric(unique(p_surn$Irradiance)))
colnames(p_surn_Tele) = 'Irradiance'
p_surn_Tele <- p_surn_Tele %>%
  arrange(Irradiance)

p_surn_Tele$p_surn = NA
p_surn_Tele$p_surn[1] <- (p_surn$Count.mL[2]/p_surn$Count.mL[1])*100
p_surn_Tele$p_surn[2] <- (p_surn$Count.mL[4]/p_surn$Count.mL[3])*100
p_surn_Tele$p_surn[3] <- (p_surn$Count.mL[6]/p_surn$Count.mL[5])*100

p_surn_mean <- p_surn_Tele %>%
  summarise(p_surn = mean(p_surn))

p_surn_sd <- p_surn_Tele %>%
  summarise(p_surn = sd(p_surn))

p_surn_stats <- bind_cols(p_surn_mean, p_surn_sd, by = NULL, suffix = c('.mean', '.sd'))

# Still need to : gather info on Mesodinium nutrition (lab notebook + excel file) and compute prey/predator ratios
#### Calculating instantaneous growth rates ####

# Instantaneous growth rate

GR <- Cyto_groups %>%
  # All the following steps arrange the dataset in order to group replicates
  # and arrange them by experimental time, so that we go from t0 to tf for each replicate
  group_by(Species, Irradiance, Replicate) %>%
  mutate(Irradiance = as.numeric(Irradiance)) %>%
  arrange(PR, .by_group = TRUE) %>%
  arrange(Irradiance, .by_group = TRUE) %>%
  arrange(Time_exp, .by_group = TRUE) %>%
  # We select only the columns that matter for the growth rate -> simplifies the table
  select(SID, Time, Date, Count.mL, Species, Irradiance, Replicate, PR, Time_exp) %>%
  # Then we can compute the instantaneous growth rate from C(t) and C(t-1)
  # Function lag() serves to get the value at t-1 (previous row in the table)
  mutate(GR_inst = ifelse(Time_exp == 0 | PR == 'Yes', 'Na', log(Count.mL/lag(Count.mL))/(Time_exp-lag(Time_exp))))

# Let's have a look : all seems fine

# Computing means and standard deviations for plots

GR_groups <- ungroup(GR) %>%
  # Getting rid of Nas
  filter(GR_inst != 'Na') %>%
  # Transforming the GR_inst values into a numeric format
  mutate(GR_inst = as.numeric(GR_inst)) %>%
  # New grouping -> for calculating stats on the GR
  group_by(Species, Irradiance, Time_exp)

GR_means <- GR_groups %>%
  summarise(GR_inst.mean = mean(GR_inst), .groups = 'keep')

GR_sds <- GR_groups %>%
  summarise(GR_inst.sd = sd(GR_inst), .groups = 'keep')

GR_stats <- left_join(GR_means, GR_sds, by = c('Species', 'Irradiance', 'Time_exp'))

# One table for each species
# Tele
GR_stats_Tele <- filter(GR_stats, Species == 'Tele')
GR_groups_Tele <- filter(GR_groups, Species == 'Tele')

# Meso
GR_stats_Meso <- filter(GR_stats, Species == 'Meso')
GR_groups_Meso <- filter(GR_groups, Species == 'Meso')

#### Plotting GRs ####

# Tele

GR_inst_Tele.plot <- ggplot(GR_stats_Tele, aes(x = Time_exp, y = GR_inst.mean, color = factor(Irradiance))) +
  # Aesthetics
  scale_color_manual(values = c('dodgerblue4', 'darkorange', 'firebrick3')) + 
  geom_point(size = 4, shape = 17) +
  geom_point(data = GR_groups_Tele, aes(x = Time_exp, y = GR_inst, color = factor(Irradiance)), 
             size = 2, shape = 17) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,12, by = 2)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Instantaneous growth rates - ', italic('T. amphioxeia'))))) +
  labs(y = c(expression(paste(mu,' (d'^'-1',')'))), x = 'Time (days)') +
  # Error bars
  geom_errorbar(aes(ymin=GR_inst.mean-GR_inst.sd, ymax=GR_inst.mean+GR_inst.sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.86, 0.52), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

GR_inst_Tele.plot

# Saving the plot as image
# ggsave('GR_inst_Tele.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

# Meso

GR_inst_Meso.plot <- ggplot(GR_stats_Meso, aes(x = Time_exp, y = GR_inst.mean, color = factor(Irradiance))) +
  # Aesthetics
  scale_color_manual(values = c('dodgerblue4', 'darkorange', 'firebrick3')) + 
  geom_point(size = 4, shape = 16) +
  geom_point(data = GR_groups_Meso, aes(x = Time_exp, y = GR_inst, color = factor(Irradiance)), 
             size = 2, shape = 16) +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,12, by = 2)) +
  # Title and axis labels
  ggtitle(c(expression(paste('Instantaneous growth rates - ', italic('M. rubrum'))))) +
  labs(y = c(expression(paste(mu,' (d'^'-1',')'))), x = 'Time (days)') +
  # Error bars
  geom_errorbar(aes(ymin=GR_inst.mean-GR_inst.sd, ymax=GR_inst.mean+GR_inst.sd), width = 0.2,
                position=position_dodge(), linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 13), 
        # Axis
        axis.title.x = element_text(size=12), axis.title.y =element_text(size=12), axis.text = element_text(size=10, color = 'black'), 
        # Legend
        legend.text = element_text(size = 13), legend.title = element_text(size = 14), 
        legend.position = c(0.27, 0.21), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = guide_legend(override.aes = list(size = 7)), shape = 'none')

GR_inst_Meso.plot

# Saving the plot as image
# ggsave('GR_inst_Meso.tiff', dpi = 600, width = 164, height = 147.6, units = 'mm')

#### Calculating GRs on defined growth intervals ####

# Based on the analysis of the semi-continuous growth curves and the instantaneous GRs,
# we chose to take the intervals (5-7 days) and (8-11 days) to calculate GR for T. amphioxeia
# and the interval (6-8 days) for M. rubrum

### Tele

Cyto_Tele_mu1 <- filter(GR, Species == 'Tele' & ((Time_exp == 5 & PR == 'Yes') | Time_exp %in% c(6,7))) %>%
  mutate(log.count = log(Count.mL)) %>%
  summarise(mu = lm(log.count ~ Time_exp)$coefficients['Time_exp'], .groups = 'keep')
# Adding the temperature variable
Cyto_Tele_mu1$Temp <- 21.0
# Growth rates for Tele in the first interval

Cyto_Tele_mu2 <- filter(GR, Species == 'Tele' & ((Time_exp == 8 & PR == 'Yes') | Time_exp %in% c(9,11))) %>%
  mutate(log.count = log(Count.mL)) %>%
  summarise(mu = lm(log.count ~ Time_exp)$coefficients['Time_exp'], .groups = 'keep')
# Adding the temperature variable
Cyto_Tele_mu2$Temp <- 21.0
# Growth rates for Tele in the second interval

## The first interval seems to be the best --> discussion needed nonetheless!

### Meso
Cyto_Meso_mu1 <- filter(GR, Species == 'Meso' & ((Time_exp == 6 & PR == 'Yes') | Time_exp %in% c(7,8))) %>%
  mutate(log.count = log(Count.mL)) %>%
  summarise(mu = lm(log.count ~ Time_exp)$coefficients['Time_exp'], .groups = 'keep')
# Adding the temperature variable
Cyto_Meso_mu1$Temp <- 21.0
# Growth rates for Meso in the only interval selected

### Data for Meso from the experiment with 4 irradiance levels
mu_Meso_2 <- read.csv2(file = 'Experiment_samples/Meso_temp_17.5.csv', header = TRUE) %>%
  group_by(Species, Temp, Irradiance)

mu_Meso_2_means <- mu_Meso_2 %>%
  summarise(mu.mean = mean(mu), .groups = 'keep')

mu_Meso_2_sds <- mu_Meso_2 %>%
  summarise(mu.sd = sd(mu), .groups = 'keep')

mu_Meso_2_stats <- left_join(mu_Meso_2_means, mu_Meso_2_sds, by = c('Species', 'Irradiance', 'Temp'))


# Now that growth rates have been computed separately for both species, we can join the data
mu_groups <- bind_rows(Cyto_Tele_mu1, Cyto_Meso_mu1, mu_Meso_2)

# Stats for plotting
# For both species together
mu_means <- ungroup(mu_groups) %>%
  group_by(Species, Temp, Irradiance) %>%
  summarise(mu.mean = mean(mu), .groups = 'keep')

mu_sds <- ungroup(mu_groups) %>%
  group_by(Species, Temp, Irradiance) %>%
  summarise(mu.sd = sd(mu), .groups = 'keep')

mu_stats <- left_join(mu_means, mu_sds, by = c('Species', 'Temp', 'Irradiance'))
mu_stats <- group_by(mu_stats, Species, Temp, Irradiance)

### Functions for growth models (based on Poisson law, after MacIntyre et al., 2002)
# Functions correspond to equation 1 in the paper

# Tele
Growth_Tele_mean <- function(x){
  0.94*(1 - exp(-(x/45.03)))
}

Growth_Tele_max <- function(x){
  0.98*(1 - exp(-(x/39.15))) #50.91
}

Growth_Tele_min <- function(x){
  0.90*(1 - exp(-(x/50.91))) #39.15
}

# Meso 21 degrees celsius
Growth_Meso21_mean <- function(x){
  0.47*(1 - exp(-(x/51.21)))
}

Growth_Meso21_max <- function(x){
  0.51*(1 - exp(-(x/37.93))) # 64.49
}

Growth_Meso21_min <- function(x){
  0.43*(1 - exp(-(x/64.49))) # 37.93
}

# Meso 17.5 degrees celsius
Growth_Meso17.5_mean <- function(x){
  0.32*(1 - exp(-(x/50.92)))
}

Growth_Meso17.5_max <- function(x){
  0.33*(1 - exp(-(x/46.33))) # 55.51
}

Growth_Meso17.5_min <- function(x){
  0.31*(1 - exp(-(x/55.51))) # 46.33
}

#### Plotting ####

# Growth rates
mu_stats.plot<-ggplot(data = mu_stats, aes(x = Irradiance, y = mu.mean, 
                                           color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3')) +
  scale_shape_discrete(c(1,2))+
  geom_point(size = 4) +
  geom_point(data = mu_groups, aes(x = Irradiance, y = mu, color = factor(Temp)), 
             size = 2, fill = 'white') +
  # Growth models
  # # Tele
  geom_function(fun = Growth_Tele_mean, color = 'firebrick3', linewidth = .5, linetype = 1) +
  geom_function(fun = Growth_Tele_max, color = 'firebrick3', linewidth = .5, linetype = 2) +
  geom_function(fun = Growth_Tele_min, color = 'firebrick3', linewidth = .5, linetype = 2) +
  # Meso 21
  geom_function(fun = Growth_Meso21_mean, color = 'firebrick3', linewidth = .5, linetype = 1) +
  geom_function(fun = Growth_Meso21_max, color = 'firebrick3', linewidth = .5, linetype = 2) +
  geom_function(fun = Growth_Meso21_min, color = 'firebrick3', linewidth = .5, linetype = 2) +
  # Meso 17.5
  geom_function(fun = Growth_Meso17.5_mean, color = 'dodgerblue4', linewidth = .5, linetype = 1) +
  geom_function(fun = Growth_Meso17.5_max, color = 'dodgerblue4', linewidth = .5, linetype = 2) +
  geom_function(fun = Growth_Meso17.5_min, color = 'dodgerblue4', linewidth = .5, linetype = 2) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,250, by = 50), limits = c(0,255)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(0,1, by = 0.25), limits = c(0,1)) +
  # Title and axis labels
  ggtitle("A") +
  labs(y = c(expression(paste(mu,' (d'^'-1',')'))), 
       x = c(expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')')))) +
  # Errorbars
  geom_errorbar(aes(ymin=mu.mean-mu.sd, ymax=mu.mean+mu.sd), width = 7,
                linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=11), axis.title.y =element_text(size=11), axis.text = element_text(size=8, color = 'black'), 
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.9, 0.75), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(color = 'none', shape = 'none')

mu_stats.plot

# Saving the plot
# ggsave('mu_plot2.tiff', dpi = 600, width = 80, height = 80, units = 'mm')

#### Biovolumes ####

Samples_photoacclim2 <- read.csv('Experiment_samples/Samples_photoacclim2.csv', sep = ';', dec = ',', header = TRUE) %>%
  select(Species, Irradiance, Replicate, Diameter)
# Adding the temperature value
Samples_photoacclim2$Temp <- 21.0

Samples_photoacclim1 <- read.csv2(file = 'Experiment_samples/Samples_photoacclim2_temp17.5.csv', header = TRUE) %>%
  mutate(Vs_CHN = as.numeric(Vs_CHN))

Biovolumes <- bind_rows(Samples_photoacclim1, Samples_photoacclim2) %>%
  group_by(Species, Temp, Irradiance) %>%
  mutate(Biovolume = (4/3)*pi*((Diameter/2)^3)) %>%
  # Creating a new variable for plotting widely different biovolumes on the same plot
  mutate(Biovolume_plot = ifelse(Species == 'Meso', Biovolume/100, Biovolume))


# Statistics for plotting

Biovolumes_means <- Biovolumes %>%
  summarise(Biovolume.mean = mean(Biovolume), .groups = 'keep')

Biovolumes_sds <- Biovolumes %>%
  summarise(Biovolume.sd = sd(Biovolume), .groups = 'keep')

Biovolumes_stats <- left_join(Biovolumes_means, Biovolumes_sds, by = c('Species', 'Temp', 'Irradiance')) %>%
  # Creating new variables for plotting widely different biovolumes and errors on the same plot
  mutate(Biovolume_plot.mean = ifelse(Species == 'Meso', Biovolume.mean/100, Biovolume.mean)) %>%
  mutate(Biovolume_plot.sd = ifelse(Species == 'Meso', Biovolume.sd/100, Biovolume.sd))


### Plotting

Biovolumes.plot<-ggplot(Biovolumes_stats, aes(x = Irradiance, y = Biovolume_plot.mean, 
                                              color = factor(Temp), shape = Species)) +
  # Aesthetics
  scale_color_discrete(type = c('dodgerblue4', 'firebrick3')) +
  scale_shape_discrete(c(1,2), labels = c(c(expression(paste(italic('M. rubrum')))),
                                          c(expression(paste(italic('T. amphioxeia')))))) +
  geom_point(size = 4) +
  geom_point(data = Biovolumes, aes(x = Irradiance, y = Biovolume_plot, color = factor(Temp)), 
             size = 2, fill = 'white') +
  geom_path(linewidth = .75) +
  # Format of x axis
  scale_x_continuous(breaks = seq(0,250, by = 50), limits = c(0,255)) +
  # Format of y axis
  scale_y_continuous(breaks = seq(20,50, by = 10), limits = c(20,52),
                     sec.axis = sec_axis(~.*100, 
                                         name = c(expression(paste(italic('M. rubrum'),' biovolume (',mu,'m'^'3',')'))))) + # c(expression(paste(italic('M. rubrum'),' biovolume (',mu,'m'^'3',')'))) -> name of the secondary axis if needed
  # Title and axis labels
  ggtitle("B") +
  labs(y = c(expression(paste(italic('T. amphioxeia'),' biovolume (',mu,'m'^'3',')'))), 
       x = (expression(paste('Irradiance (', mu,'mol photons.m'^'-2','.s'^'-1',')'))),
       color = c(expression(paste('Temperature (°C):')))) +
  # Errorbars
  geom_errorbar(aes(ymin=Biovolume_plot.mean-Biovolume_plot.sd, ymax=Biovolume_plot.mean+Biovolume_plot.sd), 
                width = 7, linewidth = .5, color = 'black') +
  theme(plot.title = element_text(size = 14, face = 'bold'), 
        # Axis
        axis.title.x = element_text(size=11), axis.title.y =element_text(size=11), axis.text = element_text(size=8, color = 'black'), 
        axis.title.y.right = element_text(angle = 90), # putting the secondary axis title in the same direction
        # Legend
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), 
        legend.position = c(0.8, 0.25), legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.key = element_blank(),
        # Panel
        panel.grid.major = element_line(color = 'gray10', linewidth = .5), panel.grid.minor = element_line(color = 'gray70', linewidth = .25),
        panel.background = element_rect(fill = 'white', color = 'white', linetype = 1))+
  guides(shape = guide_legend('Species:'))

Biovolumes.plot


# Patchwork plot
ggarrange(mu_stats.plot, Biovolumes.plot, align = 'h', ncol = 2, common.legend = TRUE, legend ='bottom')

# Saving the plot
# ggsave('mu_biov2.tiff', dpi = 600, width = 164, height = 95, units = 'mm')
