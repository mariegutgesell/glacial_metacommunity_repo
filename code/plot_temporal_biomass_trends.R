##Plotting seasonal insect community biomass

source("code/sp_matrix_arrays/invert_sp_matrix.R")
##Plot time series of family biomass for each stream 


##calculate total biomass
df_total <- df_summary %>%
  ungroup() %>%
  group_by(site_type, StreamID, Month_Year, Date) %>%
  summarise_at(vars(sp_total_biomass), list(total_biomass = sum))


ggplot() +
  geom_line(data = df_summary, aes(x = Date, y = sp_total_biomass, color = Family.x)) +
  geom_line(data = df_total, aes(x = Date, y = total_biomass)) + 
  geom_line() +
 # scale_y_log10()+
  theme_classic() +
  facet_wrap(~site_type)
##in herbert, there were no perlodidae caught in july and august, so line connects them but biomass is 0
##For plotting, add 0's in months where family was not found? 
ggplot(df_total, aes(x = Date, y = total_biomass, color = site_type)) +
  geom_line()+
  scale_y_log10() +
  theme_classic()



##calculate total biomass
df_total_std <- df_summary_std %>%
  ungroup() %>%
  group_by(site_type, StreamID, Month_Year, Date) %>%
  summarise_at(vars(sp_total_biomass_std), list(total_biomass_std = sum))


ggplot() +
  geom_line(data = df_summary_std, aes(x = Date, y = sp_total_biomass_std, color = Family.x)) +
  geom_line(data = df_total_std, aes(x = Date, y = total_biomass_std)) + 
  geom_line() +
  # scale_y_log10()+
  theme_classic() +
  facet_wrap(~site_type)
##in herbert, there were no perlodidae caught in july and august, so line connects them but biomass is 0
##For plotting, add 0's in months where family was not found? 
ggplot(df_total_std, aes(x = Date, y = total_biomass_std, color = site_type)) +
  geom_line()+
#  scale_y_log10() +
  theme_classic()
