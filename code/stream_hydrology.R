##Plotting Stream Hydrology/Physiochemical Data

library(tidyverse)
library(ggplot2)
library(lubridate)
library(ggpubr)

##Read in and clean/format data
##temperature data
temp_df <- read.csv("data/SEAK_2018_TempData.csv") %>%
  select(Stream:Temp) %>%
  mutate(site_type = case_when(
    startsWith(Stream, "Herbert") ~ "Glacier-fed",
    startsWith(Stream, "Steep") ~ "Snow-fed",
    startsWith(Stream, "Peter") ~ "Rain-fed", 
    startsWith(Stream, "Montana") ~ "Mixed",
  )) %>%
  filter(site_type != "Mixed") 

temp_df$date <- as.Date(temp_df$DateTime, "%Y-%m-%d %H:%M:%S") 

temp_df_avg <- temp_df %>%
  select(site_type, date, Temp) %>%
  group_by(site_type, date) %>%
  summarise_at(vars(Temp), list(temp_mean = mean, temp_sd = sd)) %>%
  mutate(Month_Year = format(date, "%Y-%m")) %>%
  filter(Month_Year < "2019-01") %>%
  mutate(month = format(date, "%b"),
         midmonth = as.Date(format(date, "%Y-%m-15")))

temp_month_ticks <- temp_df_avg %>%
  group_by(month) %>%
  summarize(mid_month = unique(midmonth), .groups = "drop")


#temp_df_avg$month <- months(as.Date(paste("01", temp_df_avg$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)


##discharge data
dis_df <- read.csv("data/SEAK_Discharge_Summary_Data[58].csv") %>%
  select(Date, Stream:Discharge_cmday) %>%
  filter(Stream != "Transitional")
  
 dis_df$Date <- as.Date(dis_df$Date, "%Y-%m-%d")  
  
dis_df <- dis_df %>%
  mutate(Month_Year = format(Date, "%Y-%m")) %>%
  filter(Month_Year < "2018-12") %>%
  mutate(month = format(Date, "%b"),
         midmonth = as.Date(format(Date, "%Y-%m-15")))

dis_month_ticks <- dis_df %>%
  group_by(month) %>%
  summarize(mid_month = unique(midmonth), .groups = "drop")

#dis_df$month <- months(as.Date(paste("01", dis_df$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)


##physical data 
p_df <- read.csv("data/SEAK_FoodWebs_BackgroundConditions_PointData.csv") %>%
mutate(site_type = case_when(
  startsWith(Stream, "Herbert") ~ "Glacier-fed",
  startsWith(Stream, "Steep") ~ "Snow-fed",
  startsWith(Stream, "Peter") ~ "Rain-fed", 
  startsWith(Stream, "Montana") ~ "Mixed",
)) %>%
  filter(site_type != "Mixed") %>%
  mutate(date = format(ymd_hms(Date, tz = "UTC"),"%Y-%m-%d")) %>%
  filter(date < "2018-12-01") %>%
  filter(date > "2018-03-30")

p_df$date <- as.Date(p_df$date, "%Y-%m-%d")

p_df$Turbidity <- as.numeric(p_df$Turbidity)
p_df$DOC <- as.numeric(p_df$DOC)
p_df$TDN <- as.numeric(p_df$TDN)
p_df$TP <- as.numeric(p_df$TP)

p_df <- p_df %>%
  mutate(Month_Year = format(date, "%Y-%m")) %>%
  group_by(site_type, Month_Year) %>%
  mutate(month = format(date, "%b"),
         startmonth = as.Date(min(date)))

p_df_avg <- p_df %>%
  group_by(site_type, Month_Year) %>%
  summarise_at(vars(Turbidity, DOC, TDN, TP), list(mean = mean, sd = sd))
  
p_df_avg$month <- months(as.Date(paste("01", p_df_avg$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)




##Plotting 
##what are the most important to show? look at matts stuff - temperature, runoff (discharge), turbidity, DOC, TDN and TDP
temp_df_avg$site_type <- ordered(temp_df_avg$site_type,
                             levels = c("Glacier-fed", "Snow-fed", "Rain-fed"))
temp_plot <- ggplot(temp_df_avg, aes(x = date, y = temp_mean, group = site_type, color = site_type)) +
  geom_line() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Mean Daily Temperature (˚C)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_date(breaks = temp_month_ticks$mid_month, labels = temp_month_ticks$month) 
temp_plot




##discharge 
dis_df$Stream <- ordered(dis_df$Stream,
                                 levels = c("Glacier-fed", "Snow-fed", "Rain-fed"))
dis_plot <- ggplot(dis_df, aes(x = Date, y = log(Discharge_cmday), group = Stream, color = Stream)) +
  geom_line() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Discharge (cm/day)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_date(breaks = dis_month_ticks$mid_month, labels = dis_month_ticks$month) 

dis_plot

##why NAs for some dates for rain-fed? 
p_df_avg$site_type <- ordered(p_df_avg$site_type,
                                 levels = c("Glacier-fed", "Snow-fed", "Rain-fed"))


turbidity <- ggplot(p_df_avg, aes(x = Month_Year, y = log(Turbidity_mean), group = site_type, color = site_type)) +
  geom_line() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Turbidity Log(NTU)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_discrete(labels = p_df_avg$month)
turbidity

doc<- ggplot(p_df_avg, aes(x = Month_Year, y = DOC_mean, group = site_type, color = site_type)) +
  geom_line() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Dissolved Organic Carbon (mg/L)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_discrete(labels = p_df_avg$month)
doc

tdn<- ggplot(p_df_avg, aes(x = Month_Year, y = TDN_mean, group = site_type, color = site_type)) +
  geom_line() +
  #geom_point() +
  #geom_smooth() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Total Dissolved Nitrogen (mg/L)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_discrete(labels = p_df_avg$month)
tdn


tp<- ggplot(p_df_avg, aes(x = Month_Year, y = TP_mean, group = site_type, color = site_type)) +
  geom_line() +
  #geom_point() +
  #geom_smooth() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Total Dissolved Phosphorous (mg/L)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_discrete(labels = p_df_avg$month)
tp


##smoothed lines
p_df$site_type <- ordered(p_df$site_type,
                              levels = c("Glacier-fed", "Snow-fed", "Rain-fed"))

turbidity_2 <- ggplot(p_df, aes(x = date, y = log(Turbidity), group = site_type, color = site_type)) +
  geom_point() +
  geom_smooth() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Turbidity Log(NTU)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_date(breaks = p_df$startmonth, labels = p_df$month) 

turbidity_2

doc_2<- ggplot(p_df, aes(x = date, y = DOC, group = site_type, color = site_type)) +
  geom_point() +
  geom_smooth() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Dissolved Organic Carbon (mg/L)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_date(breaks = p_df$startmonth, labels = p_df$month) 
doc_2

tdn_2<- ggplot(p_df, aes(x =date, y = TDN, group = site_type, color = site_type)) +
  #geom_line() +
  geom_point() +
  geom_smooth() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Total Dissolved Nitrogen (mg/L)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_date(breaks = p_df$startmonth, labels = p_df$month) 
tdn_2


tp_2<- ggplot(p_df, aes(x = date, y = log(TP), group = site_type, color = site_type)) +
  #geom_line() +
  geom_point() +
  geom_smooth() +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic() +
  ylab("Total Dissolved Phosphorous Log(mg/L)") +
  theme(axis.title.y = element_text(size = 12),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank())+
  scale_x_date(breaks = p_df$startmonth, labels = p_df$month) 
tp_2

##combine plots 
##Create figure legend for 3 sites 
data <- data.frame(
  Xdata = rnorm(3),
  Ydata = rnorm(3),
  LegendData = c("Glacial", "Snowmelt", "Rainfed")
)

data$LegendData <- factor(data$LegendData, levels = c("Glacial", "Snowmelt", "Rainfed"))

gplot <- ggplot(data, aes(Xdata, Ydata, color = LegendData)) +
  geom_line(size = 1) +
  scale_colour_manual(values = c("gray", "blue", "brown")) +
  theme_classic()+
  guides(colour = guide_legend(title = "")) +
  theme(legend.title = element_text(family = "Times New Roman", size = 12), legend.text = element_text(size = 12, family = "Times New Roman"), legend.position = "bottom")
gplot

leg_fig <- get_legend(gplot)

hydro_plot<- ggarrange(temp_plot, dis_plot, turbidity_2, doc_2, tdn_2, tp_2,
                     ncol = 2, nrow = 3, labels = c("a)", "b)", "c)", "d)", "e)", "f)"), font.label = list(colour = "black", size = 14, family = "Times New Roman"), legend = "bottom", common.legend = TRUE, legend.grob = leg_fig)
hydro_plot

