
library(lubridate)
library(rsample)
library(tidyverse)


#setwd("/Users/matthewdunkle/Desktop/Manuscripts/01_InProcess/Ch.2\ SEAK_SecondaryProduction/Analyses")
setwd("/project/modelscape/users/mdunkle")

SEAK_Inverts_IndMeasure_AbunTakeMeanBL = read.csv("SEAK_Inverts_IndMeasure_AbundanceTakeMeanBL.csv") %>% as_tibble() %>% 
  filter(!Sample_Number %in% c("Surber_06","Surber_07","Surber_08"))



mydata = 
  SEAK_Inverts_IndMeasure_AbunTakeMeanBL %>% 
    mutate(monthyear = paste0(lubridate::month(as.POSIXlt(Date, format = "%m/%d/%y"), abbr=T, label = T),  "_",lubridate::year(as.POSIXlt(Date, format = "%m/%d/%y")))) %>% 
    nest_by(monthyear,StreamID, Family.x, .keep = T)
  
# mydata$data[1]

bstimes = 1000  #change this to 1000 for true run
maxbs = 100000

results = data.frame()
output <- data.frame()

for(j in 1:length(mydata$data)){
boot = rsample::bootstraps(mydata$data[[j]],breaks = 5, times = bstimes, strata = Sample_Number)
  
  
      for(i in 1:bstimes){
        tmp = boot$splits[[i]] %>% assessment() %>% mutate(bs_num = i)  
        assign(x="results", value=rbind(get("results"), tmp)) 
        #results = rbind(results, tmp)
        
                
        }
if( (j%%2==0) | (j==length(mydata$data)) )  #every 2k sites write out the data to disk to save RAM! 
{
  firstrow <- j-1
  lastrow <- j
  
  #output = results %>% mutate(BodyLength = round(BodyLength/ 0.5)*.5) %>% group_by(StreamID, Date, Family.x, bs_num, BodyLength) %>% 
  # dplyr::summarise(n = n(), Abundance = n/.75, Biomass_ind = mean(Biomass_Ind)) %>% 
  #group_by(StreamID, Date, Family.x, BodyLength) %>% 
  #dplyr::summarise(Abun_mean = mean(Abundance, na.rm=T), Abun_sd = sd(Abundance, na.rm=T),
  #                 Abun_se = Abun_sd/sqrt(n()),
  #                 lowCI = Abun_mean-1.96*Abun_se,
  #                 hiCI = Abun_mean+1.96*Abun_se,
  #                 Biom_mean = mean(Biomass_ind)) %>%
  #mutate(Biomass = Biom_mean*Abun_mean,
  #       Biom_lo =Biom_mean*lowCI,
  #       Biom_hi = Biom_mean*hiCI) %>%
  #distinct(across(everything()))
  
  write.csv(results, file=paste0("SEAK_bootstrap","_rows_",firstrow,"-",lastrow,".csv"), quote=FALSE, row.names=FALSE, sep=",")
  #output <- data.frame()
  results = data.frame()
  
       }
}
# 
# 
# 
# output = results %>% mutate(BodyLength = round(BodyLength/ 0.5)*.5) %>% group_by(StreamID, Date, Family.x, bs_num, BodyLength) %>% 
#   dplyr::summarise(n = n(), Abundance = n/.75, Biomass_ind = mean(Biomass_Ind)) %>% 
#   group_by(StreamID, Date, Family.x, BodyLength) %>% 
#   dplyr::summarise(Abun_mean = mean(Abundance, na.rm=T), Abun_sd = sd(Abundance, na.rm=T),
#                    Abun_se = Abun_sd/sqrt(n()),
#                    lowCI = Abun_mean-1.96*Abun_se,
#                     hiCI = Abun_mean+1.96*Abun_se,
#                    Biom_mean = mean(Biomass_ind)) %>%
#     mutate(Biomass = Biom_mean*Abun_mean,
#                     Biom_lo =Biom_mean*lowCI,
#                   Biom_hi = Biom_mean*hiCI) %>%
#   distinct(across(everything()))
# 
# output = write.csv("SEAK_Bootstrap_AbundanceBiomass.csv")
# 
# 
# output %>% group_by(StreamID, Date, Family.x) %>%
#   dplyr::summarise(Abundance = sum(Abun_mean), Biomass = sum(Biomass), Biom_low = sum(Biom_lo, na.rm=T), Biom_hi = sum(Biom_hi, na.rm=T)) %>%
#   mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
#   ggplot(aes(x=Date, y = Biomass, fill = Family.x))+geom_area()+facet_wrap(~StreamID, scales = "free")+theme(legend.position = "none", panel.background = element_blank())


