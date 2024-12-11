#Maya Powell
#Curacao Environmental data
#DHW Graph
#Summary and time data
#all environmental data stats
#Code for DHW graph from Sarah Solomon
#ALL environmental data from Verena Schoepf's lab:
#2020 data can be found in de Jong et al 2023 preprint and 2021 data from Sarah Solomon
#Graphs of time series and averages
#Stats for daily averages

#Setup
setwd("~/Documents/Castillo Lab/CW_2020/Environmental_Data")
library(dplyr)
#install.packages('gt')
library(gt)
library(ggplot2)
library(rstatix)
library(RColorBrewer)

######DHW and Time Series Temp Data
#copied from Sarah Solomon's data - creating the same graph that she did here
#wooh HUGE shoutout to Sarah for concatenating all NOAA and aqualink data to create this graph!!
dhw_data<-read.csv("dhw_CW2020.csv")

plot2 <- ggplot(dhw_data) +
  geom_line(size=1.1, aes(x = Day, y = Daily.average), colour = "#0000FF") +
  geom_hline(yintercept= 29, linetype="dashed", color="orange",size=0.75) +
  geom_hline(yintercept= 28, linetype="dashed", color="black",size=0.75) +
  geom_vline(xintercept= 72, linetype="dotted", color="gray",size=0.5) +
  geom_vline(xintercept= 81, linetype="dotted", color="gray",size=0.5) +
  geom_vline(xintercept= 305, linetype="dotted", color="gray",size=0.5) +
  geom_vline(xintercept= 336, linetype="dotted", color="gray",size=0.5) +
  geom_vline(xintercept= 663, linetype="dotted", color="gray",size=0.5) +
  geom_vline(xintercept= 695, linetype="dotted", color="gray",size=0.5) +
  scale_x_continuous(breaks = seq(0, 720, by = 60), position = "bottom") +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    panel.border = element_rect(color = "black", fill = NA), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    strip.background = element_blank()
  )
plot2
#removed this code
#geom_segment(aes(x = Day, xend = Day, y = Min, yend = mean), linetype = "dotdash", size = 0.8, alpha = 0.5) +  # Add vertical line for minimum
#geom_segment(aes(x = Day, xend = Day, y = Max, yend = mean), linetype = "dotdash", size = 0.8, alpha = 0.5) +  # Add vertical line for maximum
ggsave(plot2, file = "dailyavgtemp.pdf", dpi="retina", w = 20, h = 8)

plot3 <- ggplot(dhw_data) +
  geom_line(size=1.1, aes(x = Day, y = DHW.per.day), colour = "red") +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  scale_x_continuous(breaks = seq(0, 800, by = 60), position = "bottom") +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    panel.border = element_rect(color = "black", fill = NA), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    strip.background = element_blank()
  )
plot3
ggsave(plot3, file = "dhw.pdf", dpi="retina", w = 20, h = 8)
#then put graphs together in illustrator

#pH data
#creating combined excel sheet with all data across sites and seasons
#March 2020 data from folder "Calibrated 24hr cycle data"
#Nov 2020 data from "Calibrated 24hr cycle data"
#Nov 2021 data from "All useable pH data from Nov 21"

pH_data <-read.csv("all_pH_CW2020.csv")
pH_data$Date_Time <- paste(pH_data$Date, pH_data$Time)
pH_data$Date_Time = as.POSIXct(pH_data$Date_Time, format = "%m/%d/%Y %H:%M")

pH_march20 <- pH_data %>% filter(Timepoint == "March 2020")
pH_nov20 <- pH_data %>% filter(Timepoint == "November 2020")
pH_nov21 <- pH_data %>% filter(Timepoint == "November 2021")

#MARCH 2020
pH_SWB_mar20 <- pH_march20 %>% filter(Site == "Spaanse Water Bay")
pH_SWR_mar20 <- pH_march20 %>% filter(Site == "Spaanse Water Reef")
pH_SMB_mar20 <- pH_march20 %>% filter(Site == "Santa Martha Bay")
pH_SMR_mar20 <- pH_march20 %>% filter(Site == "Santa Martha Reef")

all.march20 <- data.frame(Date_Time=seq(min(pH_march20$Date_Time), 
                                             max(pH_march20$Date_Time), by="5 min"))

pH_SWB_mar20_merge <- merge(all.march20, pH_SWB_mar20, by = "Date_Time", all=TRUE)
pH_SWR_mar20_merge <- merge(all.march20, pH_SWR_mar20, by = "Date_Time", all=TRUE)
pH_SMB_mar20_merge <- merge(all.march20, pH_SMB_mar20, by = "Date_Time", all=TRUE)
pH_SMR_mar20_merge <- merge(all.march20, pH_SMR_mar20, by = "Date_Time", all=TRUE)

mar20_lims <- as.POSIXct(strptime(c("0020-03-12 00:00:00", "0020-03-22 00:00:00"), format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#pH March 2020
pH_line_mar20 <- ggplot(NULL, aes(x = Date_Time, y = pHT, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = pH_SMB_mar20, col = "coral3", size = 0.5)+
  geom_point(data = pH_SMR_mar20, col = "paleturquoise", size = 0.5)+
  geom_point(data = pH_SWB_mar20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = pH_SWR_mar20, col = "deepskyblue4", size = 0.5)+
  scale_x_datetime(limits = mar20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("pHT")+
  scale_y_continuous(limits = c(7.7, 8.3))
pH_line_mar20

ggsave(pH_line_mar20, file = "pH_dot_mar20_final.pdf", w=4, h=6)

#pH NOVEMBER 2020
pH_SWB_nov20 <- pH_nov20 %>% filter(Site == "Spaanse Water Bay")
pH_SWR_nov20 <- pH_nov20 %>% filter(Site == "Spaanse Water Reef")
pH_SMB_nov20 <- pH_nov20 %>% filter(Site == "Santa Martha Bay")
pH_SMR_nov20 <- pH_nov20 %>% filter(Site == "Santa Martha Reef")

all.nov20 <- data.frame(Date_Time=seq(min(pH_nov20$Date_Time), 
                                        max(pH_nov20$Date_Time), by="15 min"))

pH_SWB_nov20_merge <- merge(all.nov20, pH_SWB_nov20, by = "Date_Time", all=TRUE)
pH_SWR_nov20_merge <- merge(all.nov20, pH_SWR_nov20, by = "Date_Time", all=TRUE)
pH_SMB_nov20_merge <- merge(all.nov20, pH_SMB_nov20, by = "Date_Time", all=TRUE)
pH_SMR_nov20_merge <- merge(all.nov20, pH_SMR_nov20, by = "Date_Time", all=TRUE)

nov20_lims <- as.POSIXct(strptime(c("0020-10-30 10:43:00", "0020-12-02 12:48:00"), format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#trying to just combine jeez
pH_line_nov20 <- ggplot(NULL, aes(x = Date_Time, y = pHT, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = pH_SMB_nov20, col = "coral3", size = 0.5)+
  geom_point(data = pH_SMR_nov20, col = "paleturquoise", size = 0.5)+
  geom_point(data = pH_SWB_nov20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = pH_SWR_nov20, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("pHT")+
  scale_y_continuous(limits = c(7.7, 8.3))
pH_line_nov20
ggsave(pH_line_nov20, file = "pH_dot_nov20_final.pdf", w=10, h=6)

##pH NOVEMBER 2021
pH_SWB_nov21 <- pH_nov21 %>% filter(Site == "Spaanse Water Bay")
pH_SWR_nov21 <- pH_nov21 %>% filter(Site == "Spaanse Water Reef")
pH_SMB_nov21 <- pH_nov21 %>% filter(Site == "Santa Martha Bay")
pH_SMR_nov21 <- pH_nov21 %>% filter(Site == "Santa Martha Reef")

all.nov21 <- data.frame(Date_Time=seq(min(pH_nov21$Date_Time), 
                                      max(pH_nov21$Date_Time), by="1 min"))

limits <- c(min(pH_nov21$Date_Time),max(pH_nov21$Date_Time))
nov21_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

#"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#trying to just combine jeez
pH_line_nov21 <- ggplot(NULL, aes(x = Date_Time, y = pHT, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = pH_SMB_nov21, col = "coral3", size = 0.5)+
  geom_point(data = pH_SMR_nov21, col = "paleturquoise", size = 0.5)+
  geom_point(data = pH_SWB_nov21, col = "lightsalmon1", size = 0.5)+
  geom_point(data = pH_SWR_nov21, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov21_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("pHT")+
  scale_y_continuous(limits = c(7.7, 8.3))
pH_line_nov21
ggsave(pH_line_nov21, file = "pH_dot_nov21_final.pdf", w=10, h=6)

pH_time <- ggarrange(pH_line_mar20,pH_line_nov20,pH_line_nov21, nrow = 1, ncol = 3, widths = c(1,3,3))
pH_time
ggsave(pH_time, file = "pH_time.pdf", w=30, h=6)

#
#DO data
#creating combined excel sheet with all data across sites and seasons
#March 2020 data from folder "DO_FINAL" from each site folder
#Nov 2020 data from "DO and Temp - use cut data"
#Nov 2021 data from "All useable DO data from Nov 2021"

do_data <-read.csv("all_do_CW2020.csv")
do_data$Date_Time <- paste(do_data$Date, do_data$Time_mil)
do_data$Date_Time = as.POSIXct(do_data$Date_Time, format = "%m/%d/%Y %H:%M")

#MARCH 2020
do_march20 <- do_data %>% filter(Timepoint == "March 2020")

do_SWB_mar20 <- do_march20 %>% filter(Site == "Spaanse Water Bay")
do_SWR_mar20 <- do_march20 %>% filter(Site == "Spaanse Water Reef")
do_SMB_mar20 <- do_march20 %>% filter(Site == "Santa Martha Bay")
do_SMR_mar20 <- do_march20 %>% filter(Site == "Santa Martha Reef")

"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#do March 2020
limits <- c("2020-03-06 15:10:00","2020-03-21 20:50:00")
march20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

do_point_mar20 <- ggplot(NULL, aes(x = Date_Time, y = DO_Conc, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = do_SMB_mar20, col = "coral3", size = 0.5)+
  geom_point(data = do_SMR_mar20, col = "paleturquoise", size = 0.5)+
  geom_point(data = do_SWB_mar20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = do_SWR_mar20, col = "deepskyblue4", size = 0.5)+
  scale_x_datetime(limits = march20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("Dissolved Oxygen (mg L-1)")+
  scale_y_continuous(limits = c(2, 10))
do_point_mar20
ggsave(do_point_mar20, file = "do_dot_mar20_final.pdf", w=15, h=6)

#do NOVEMBER 2020
do_nov20 <- do_data %>% filter(Timepoint == "November 2020")

do_SWB_nov20 <- do_nov20 %>% filter(Site == "Spaanse Water Bay")
do_SWR_nov20 <- do_nov20 %>% filter(Site == "Spaanse Water Reef")
do_SMB_nov20 <- do_nov20 %>% filter(Site == "Santa Martha Bay")
do_SMR_nov20 <- do_nov20 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(do_nov20$Date_Time),max(do_nov20$Date_Time))
nov20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#trying to just combine jeez
do_line_nov20 <- ggplot(NULL, aes(x = Date_Time, y = DO_Conc, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = do_SMB_nov20, col = "coral3", size = 0.5)+
  geom_point(data = do_SMR_nov20, col = "paleturquoise", size = 0.5)+
  geom_point(data = do_SWB_nov20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = do_SWR_nov20, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("do")+
  scale_y_continuous(limits = c(2, 10))
do_line_nov20
ggsave(do_line_nov20, file = "do_do_nov20_final.pdf", w=20, h=6)

##do NOVEMBER 2021
do_nov21 <- do_data %>% filter(Timepoint == "November 2021")
do_SWB_nov21 <- do_nov21 %>% filter(Site == "Spaanse Water Bay")
do_SWR_nov21 <- do_nov21 %>% filter(Site == "Spaanse Water Reef")
do_SMB_nov21 <- do_nov21 %>% filter(Site == "Santa Martha Bay")
do_SMR_nov21 <- do_nov21 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(do_nov21$Date_Time),max(do_nov21$Date_Time))

nov21_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

#"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#trying to just combine jeez
do_line_nov21 <- ggplot(NULL, aes(x = Date_Time, y = DO_Conc, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = do_SMB_nov21, col = "coral3", size = 0.5)+
  geom_point(data = do_SMR_nov21, col = "paleturquoise", size = 0.5)+
  geom_point(data = do_SWB_nov21, col = "lightsalmon1", size = 0.5)+
  geom_point(data = do_SWR_nov21, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov21_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("do")+
  scale_y_continuous(limits = c(2, 10))
do_line_nov21
ggsave(do_line_nov21, file = "do_do_nov21_final.pdf", w=20, h=6)

#time
do_time <- ggarrange(do_point_mar20,do_line_nov20,do_line_nov21, nrow = 1, ncol = 3, widths = c(1,2.2,2.2))
do_time
ggsave(do_time, file = "do_time.pdf", w=30, h=6)

#sal data
#creating combined excel sheet with all data across sites and seasons

sal_data <-read.csv("all_sal_CW2020.csv")
sal_data$Date_Time <- paste(sal_data$Date, sal_data$Time)
sal_data$Date_Time = as.POSIXct(sal_data$Date_Time, format = "%m/%d/%Y %H:%M")

#MARCH 2020
sal_march20 <- sal_data %>% filter(Timepoint == "March 2020")

sal_SWB_mar20 <- sal_march20 %>% filter(Site == "Spaanse Water Bay")
sal_SWR_mar20 <- sal_march20 %>% filter(Site == "Spaanse Water Reef")
sal_SMB_mar20 <- sal_march20 %>% filter(Site == "Santa Martha Bay")
sal_SMR_mar20 <- sal_march20 %>% filter(Site == "Santa Martha Reef")

"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#sal March 2020
limits <- c(min(sal_march20$Date_Time),max(sal_march20$Date_Time))
march20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

sal_point_mar20 <- ggplot(NULL, aes(x = Date_Time, y = Salinity, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = sal_SMB_mar20, col = "coral3", size = 0.5)+
  geom_point(data = sal_SMR_mar20, col = "paleturquoise", size = 0.5)+
  geom_point(data = sal_SWB_mar20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = sal_SWR_mar20, col = "deepskyblue4", size = 0.5)+
  scale_x_datetime(limits = march20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("Salinity")+
  scale_y_continuous(limits = c(27, 38), n.breaks = 10)
#scale_y_continuous(limits = c(7.8, 8.2))
sal_point_mar20
ggsave(sal_point_mar20, file = "sal_salt_mar20_final.pdf", w=15, h=6)

#sal NOVEMBER 2020
sal_nov20 <- sal_data %>% filter(Timepoint == "November 2020")

sal_SWB_nov20 <- sal_nov20 %>% filter(Site == "Spaanse Water Bay")
sal_SWR_nov20 <- sal_nov20 %>% filter(Site == "Spaanse Water Reef")
sal_SMB_nov20 <- sal_nov20 %>% filter(Site == "Santa Martha Bay")
sal_SMR_nov20 <- sal_nov20 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(sal_nov20$Date_Time),max(sal_nov20$Date_Time))
nov20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#trying to just combine jeez
sal_line_nov20 <- ggplot(NULL, aes(x = Date_Time, y = Salinity, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = sal_SMB_nov20, col = "coral3", size = 0.5)+
  geom_point(data = sal_SMR_nov20, col = "paleturquoise", size = 0.5)+
  geom_point(data = sal_SWB_nov20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = sal_SWR_nov20, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("Salinity")+
  scale_y_continuous(limits = c(27, 38), n.breaks = 10)
#scale_y_continuous(limits = c(7.8, 8.2))
sal_line_nov20
ggsave(sal_line_nov20, file = "sal_sal_nov20_final.pdf", w=20, h=6)

##sal NOVEMBER 2021
sal_nov21 <- sal_data %>% filter(Timepoint == "November 2021")
sal_SWB_nov21 <- sal_nov21 %>% filter(Site == "Spaanse Water Bay")
sal_SWR_nov21 <- sal_nov21 %>% filter(Site == "Spaanse Water Reef")
sal_SMB_nov21 <- sal_nov21 %>% filter(Site == "Santa Martha Bay")
sal_SMR_nov21 <- sal_nov21 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(sal_nov21$Date_Time),max(sal_nov21$Date_Time))

nov21_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

#"coral3", "paleturquoise","lightsalmon1","deepskyblue4"

#trying to just combine jeez
sal_line_nov21 <- ggplot(NULL, aes(x = Date_Time, y = Salinity, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = sal_SMB_nov21, col = "coral3", size = 0.5)+
  geom_point(data = sal_SMR_nov21, col = "paleturquoise", size = 0.5)+
  geom_point(data = sal_SWB_nov21, col = "lightsalmon1", size = 0.5)+
  geom_point(data = sal_SWR_nov21, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov21_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("sal")+
  scale_y_continuous(limits = c(27, 38), n.breaks = 10)
#scale_y_continuous(limits = c(7.8, 8.2))
sal_line_nov21
ggsave(sal_line_nov21, file = "sal_sal_nov21_final.pdf", w=20, h=6)

#sal_time
sal_time <- ggarrange(sal_point_mar20,sal_line_nov20,sal_line_nov21, nrow = 1, ncol = 3, widths = c(1,2.2,2.2))
sal_time
ggsave(sal_time, file = "sal_time.pdf", w=30, h=6)

#par data
#creating combined excel sheet with all data across sites and seasons

par_data <-read.csv("all_par_CW2020.csv")
par_data$Date_Time <- paste(par_data$Date, par_data$Time)
par_data$Date_Time = as.POSIXct(par_data$Date_Time, format = "%m/%d/%Y %H:%M")

#MARCH 2020
par_march20 <- par_data %>% filter(Timepoint == "March 2020")

par_SWB_mar20 <- par_march20 %>% filter(Site == "Spaanse Water Bay")
par_SWR_mar20 <- par_march20 %>% filter(Site == "Spaanse Water Reef")
par_SMB_mar20 <- par_march20 %>% filter(Site == "Santa Martha Bay")
par_SMR_mar20 <- par_march20 %>% filter(Site == "Santa Martha Reef")

"coral3", "paleturquoise","lightparmon1","deepskyblue4"

#par March 2020
limits <- c(min(par_march20$Date_Time),max(par_march20$Date_Time))
march20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

par_point_mar20 <- ggplot(NULL, aes(x = Date_Time, y = PAR_cal, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = par_SMB_mar20, col = "coral3", size = 0.5)+
  geom_point(data = par_SMR_mar20, col = "paleturquoise", size = 0.5)+
  geom_point(data = par_SWB_mar20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = par_SWR_mar20, col = "deepskyblue4", size = 0.5)+
  scale_x_datetime(limits = march20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("PAR")+
  scale_y_continuous(limits = c(0, 2200))
par_point_mar20
ggsave(par_point_mar20, file = "par_dot_mar20_final.pdf", w=15, h=6)

#par NOVEMBER 2020
par_nov20 <- par_data %>% filter(Timepoint == "November 2020")

par_SWB_nov20 <- par_nov20 %>% filter(Site == "Spaanse Water Bay")
par_SWR_nov20 <- par_nov20 %>% filter(Site == "Spaanse Water Reef")
par_SMB_nov20 <- par_nov20 %>% filter(Site == "Santa Martha Bay")
par_SMR_nov20 <- par_nov20 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(par_nov20$Date_Time),max(par_nov20$Date_Time))
nov20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lightparmon1","deepskyblue4"

#trying to just combine jeez
par_line_nov20 <- ggplot(NULL, aes(x = Date_Time, y = PAR_cal, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = par_SMB_nov20, col = "coral3", size = 0.5)+
  geom_point(data = par_SMR_nov20, col = "paleturquoise", size = 0.5)+
  geom_point(data = par_SWB_nov20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = par_SWR_nov20, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("parinity")+
  scale_y_continuous(limits = c(0, 2200))
#scale_y_continuous(limits = c(7.8, 8.2))
par_line_nov20
ggsave(par_line_nov20, file = "par_dot_nov20_final.pdf", w=20, h=6)

##par NOVEMBER 2021
par_nov21 <- par_data %>% filter(Timepoint == "November 2021")
par_SWB_nov21 <- par_nov21 %>% filter(Site == "Spaanse Water Bay")
par_SWR_nov21 <- par_nov21 %>% filter(Site == "Spaanse Water Reef")
par_SMB_nov21 <- par_nov21 %>% filter(Site == "Santa Martha Bay")
par_SMR_nov21 <- par_nov21 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(par_nov21$Date_Time),max(par_nov21$Date_Time))

nov21_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

#"coral3", "paleturquoise","lightparmon1","deepskyblue4"

#trying to just combine jeez
par_line_nov21 <- ggplot(NULL, aes(x = Date_Time, y = PAR_cal, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = par_SMB_nov21, col = "coral3", size = 0.5)+
  geom_point(data = par_SMR_nov21, col = "paleturquoise", size = 0.5)+
  geom_point(data = par_SWB_nov21, col = "lightsalmon1", size = 0.5)+
  geom_point(data = par_SWR_nov21, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov21_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("par")+
  scale_y_continuous(limits = c(0, 2200))
#scale_y_continuous(limits = c(7.8, 8.2))
par_line_nov21
ggsave(par_line_nov21, file = "par_dot_nov21_final.pdf", w=20, h=6)

##par_time
par_time <- ggarrange(par_point_mar20,par_line_nov20,par_line_nov21, nrow = 1, ncol = 3, widths = c(1,2.2,2.2))
par_time
ggsave(par_time, file = "par_time.pdf", w=30, h=6)


#TEMP DATA TWO WAYS - either with pH temperature, OR with dissolved oxygen temperature
#ENDED UP USING DO DATA AS FINAL TEMP DATA BC OF MORE DATAPOINTS/COMPREHENSIVE SAMPLING

#temp data from pH data
#creating combined excel sheet with all data across sites and seasons

temp_data <-read.csv("all_pH_CW2020.csv")
temp_data$Date_Time <- paste(temp_data$Date, temp_data$Time)
temp_data$Date_Time = as.POSIXct(temp_data$Date_Time, format = "%m/%d/%Y %H:%M")

#MARCH 2020
temp_march20 <- temp_data %>% filter(Timepoint == "March 2020")

temp_SWB_mar20 <- temp_march20 %>% filter(Site == "Spaanse Water Bay")
temp_SWR_mar20 <- temp_march20 %>% filter(Site == "Spaanse Water Reef")
temp_SMB_mar20 <- temp_march20 %>% filter(Site == "Santa Martha Bay")
temp_SMR_mar20 <- temp_march20 %>% filter(Site == "Santa Martha Reef")

"coral3", "paleturquoise","lighttempmon1","deepskyblue4"

#temp March 2020
limits <- c(min(temp_march20$Date_Time),max(temp_march20$Date_Time))
march20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

temp_point_mar20 <- ggplot(NULL, aes(x = Date_Time, y = Temp_C, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = temp_SMB_mar20, col = "coral3", size = 0.5)+
  geom_point(data = temp_SMR_mar20, col = "paleturquoise", size = 0.5)+
  geom_point(data = temp_SWB_mar20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = temp_SWR_mar20, col = "deepskyblue4", size = 0.5)+
  scale_x_datetime(limits = march20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("temp")+
  scale_y_continuous(limits = c(26, 31), n.breaks = 6)
temp_point_mar20
ggsave(temp_point_mar20, file = "temp_dot_mar20_pH_logger_final.pdf", w=15, h=6)

#temp NOVEMBER 2020
temp_nov20 <- temp_data %>% filter(Timepoint == "November 2020")

temp_SWB_nov20 <- temp_nov20 %>% filter(Site == "Spaanse Water Bay")
temp_SWR_nov20 <- temp_nov20 %>% filter(Site == "Spaanse Water Reef")
temp_SMB_nov20 <- temp_nov20 %>% filter(Site == "Santa Martha Bay")
temp_SMR_nov20 <- temp_nov20 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(temp_nov20$Date_Time),max(temp_nov20$Date_Time))
nov20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lighttempmon1","deepskyblue4"

#trying to just combine jeez
temp_line_nov20 <- ggplot(NULL, aes(x = Date_Time, y = Temp_C, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = temp_SMB_nov20, col = "coral3", size = 0.5)+
  geom_point(data = temp_SMR_nov20, col = "paleturquoise", size = 0.5)+
  geom_point(data = temp_SWB_nov20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = temp_SWR_nov20, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("tempinity")+
  scale_y_continuous(limits = c(26, 31), n.breaks = 6)
temp_line_nov20
ggsave(temp_line_nov20, file = "temp_dot_nov20_pH_logger_final.pdf", w=20, h=6)

##temp NOVEMBER 2021
temp_nov21 <- temp_data %>% filter(Timepoint == "November 2021")
temp_SWB_nov21 <- temp_nov21 %>% filter(Site == "Spaanse Water Bay")
temp_SWR_nov21 <- temp_nov21 %>% filter(Site == "Spaanse Water Reef")
temp_SMB_nov21 <- temp_nov21 %>% filter(Site == "Santa Martha Bay")
temp_SMR_nov21 <- temp_nov21 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(temp_nov21$Date_Time),max(temp_nov21$Date_Time))

nov21_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

#"coral3", "paleturquoise","lighttempmon1","deepskyblue4"

#trying to just combine jeez
temp_line_nov21 <- ggplot(NULL, aes(x = Date_Time, y = Temp_C, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = temp_SMB_nov21, col = "coral3", size = 0.5)+
  geom_point(data = temp_SMR_nov21, col = "paleturquoise", size = 0.5)+
  geom_point(data = temp_SWB_nov21, col = "lightsalmon1", size = 0.5)+
  geom_point(data = temp_SWR_nov21, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov21_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("temp")+
  scale_y_continuous(limits = c(26, 31), n.breaks = 6)
temp_line_nov21
ggsave(temp_line_nov21, file = "temp_dot_nov21_pH_logger_final.pdf", w=20, h=6)

##temp_pH_time
temp_pH_time <- ggarrange(temp_point_mar20,temp_line_nov20,temp_line_nov21, nrow = 1, ncol = 3, widths = c(1,3,3))
temp_pH_time
ggsave(temp_pH_time, file = "temp_pH_time.pdf", w=30, h=6)


#temp data from DO data
temp_data <-read.csv("all_DO_CW2020.csv")
temp_data$Date_Time <- paste(temp_data$Date, temp_data$Time_mil)
temp_data$Date_Time = as.POSIXct(temp_data$Date_Time, format = "%m/%d/%Y %H:%M")
#this data contains one clear outlier of 34C - removing in the plot

#MARCH 2020
temp_march20 <- temp_data %>% filter(Timepoint == "March 2020")

temp_SWB_mar20 <- temp_march20 %>% filter(Site == "Spaanse Water Bay")
temp_SWR_mar20 <- temp_march20 %>% filter(Site == "Spaanse Water Reef")
temp_SMB_mar20 <- temp_march20 %>% filter(Site == "Santa Martha Bay")
temp_SMR_mar20 <- temp_march20 %>% filter(Site == "Santa Martha Reef")

"coral3", "paleturquoise","lighttempmon1","deepskyblue4"

#temp March 2020
limits <- c(min(temp_march20$Date_Time),max(temp_march20$Date_Time))
march20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

temp_point_mar20 <- ggplot(NULL, aes(x = Date_Time, y = Temp, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = temp_SMB_mar20, col = "coral3", size = 0.5)+
  geom_point(data = temp_SMR_mar20, col = "paleturquoise", size = 0.5)+
  geom_point(data = temp_SWB_mar20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = temp_SWR_mar20, col = "deepskyblue4", size = 0.5)+
  scale_x_datetime(limits = march20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("temp")+
  #scale_y_continuous(limits = c(25, 34), n.breaks = 6)
  scale_y_continuous(limits = c(26, 31), n.breaks = 6)
temp_point_mar20
ggsave(temp_point_mar20, file = "temp_dot_mar20_DO_logger_final.pdf", w=15, h=6)

#temp NOVEMBER 2020
temp_nov20 <- temp_data %>% filter(Timepoint == "November 2020")

temp_SWB_nov20 <- temp_nov20 %>% filter(Site == "Spaanse Water Bay")
temp_SWR_nov20 <- temp_nov20 %>% filter(Site == "Spaanse Water Reef")
temp_SMB_nov20 <- temp_nov20 %>% filter(Site == "Santa Martha Bay")
temp_SMR_nov20 <- temp_nov20 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(temp_nov20$Date_Time),max(temp_nov20$Date_Time))
nov20_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

"coral3", "paleturquoise","lighttempmon1","deepskyblue4"

#trying to just combine jeez
temp_line_nov20 <- ggplot(NULL, aes(x = Date_Time, y = Temp, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = temp_SMB_nov20, col = "coral3", size = 0.5)+
  geom_point(data = temp_SMR_nov20, col = "paleturquoise", size = 0.5)+
  geom_point(data = temp_SWB_nov20, col = "lightsalmon1", size = 0.5)+
  geom_point(data = temp_SWR_nov20, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov20_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("tempinity")+
  scale_y_continuous(limits = c(26, 31), n.breaks = 6)
temp_line_nov20
ggsave(temp_line_nov20, file = "temp_dot_nov20_do_logger_final.pdf", w=20, h=6)

##temp NOVEMBER 2021
temp_nov21 <- temp_data %>% filter(Timepoint == "November 2021")
temp_SWB_nov21 <- temp_nov21 %>% filter(Site == "Spaanse Water Bay")
temp_SWR_nov21 <- temp_nov21 %>% filter(Site == "Spaanse Water Reef")
temp_SMB_nov21 <- temp_nov21 %>% filter(Site == "Santa Martha Bay")
temp_SMR_nov21 <- temp_nov21 %>% filter(Site == "Santa Martha Reef")

limits <- c(min(temp_nov21$Date_Time),max(temp_nov21$Date_Time))

nov21_lims <- as.POSIXct(strptime(limits, format = "%Y/%m/%d %H:%M:%S"))

#"coral3", "paleturquoise","lighttempmon1","deepskyblue4"

#trying to just combine jeez
temp_line_nov21 <- ggplot(NULL, aes(x = Date_Time, y = Temp, group = 1)) +
  theme_classic(base_size = 22)+
  geom_point(data = temp_SMB_nov21, col = "coral3", size = 0.5)+
  geom_point(data = temp_SMR_nov21, col = "paleturquoise", size = 0.5)+
  geom_point(data = temp_SWB_nov21, col = "lightsalmon1", size = 0.5)+
  geom_point(data = temp_SWR_nov21, col = "deepskyblue4", size = 0.5)+
  #scale_color_manual(values = c("deepskyblue4"))+
  scale_x_datetime(limits = nov21_lims, date_breaks = "3 day", date_labels = "%m/%d/%Y")+
  theme(text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Date")+
  ylab("temp")+
  scale_y_continuous(limits = c(26, 31), n.breaks = 6)
temp_line_nov21
ggsave(temp_line_nov21, file = "temp_dot_nov21_do_logger_final.pdf", w=20, h=6)

##temp_do_time
temp_do_time <- ggarrange(temp_point_mar20,temp_line_nov20,temp_line_nov21, nrow = 1, ncol = 3, widths = c(1,2.2,2.2))
temp_do_time
ggsave(temp_do_time, file = "temp_do_time.pdf", w=30, h=6)

#AVERAGES

pH_data <-read.csv("all_pH_CW2020.csv")
do_data <-read.csv("all_DO_CW2020.csv")
sal_data <-read.csv("all_sal_CW2020.csv")
par_data <-read.csv("all_PAR_CW2020.csv")
temp_data <-read.csv("all_DO_CW2020.csv")

library(see)

#pH
pH.mean <- pH_data %>% ggplot(aes(x = Site, y=pHT, fill = Site)) +
  theme_classic(base_size = 22)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  #stat_summary(geom = 'text', label = c("a", "b","c","a", "e","f","g","h"), fun = max, vjust = -0.5, size=6)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  geom_violinhalf()+
  stat_summary(fun = mean, geom = "point", size = 4) + 
  #stat_summary(fun = range, geom = "line", linewidth = 2, alpha = 0.6)+
  #geom_jitter(width = 0.2, alpha = 0.8)+
  #geom_boxplot(alpha = 0)+
  ylab("pH")+
  facet_wrap(Timepoint~.)+
  scale_y_continuous(limits = c(7.7, 8.3))
pH.mean
ggsave(pH.mean,file="pH.mean.range.pdf",h=6,w=12)
ggsave(pH.mean,file="pH.mean.violinhalf.pdf",h=6,w=12)

#do
do.mean <- do_data %>% ggplot(aes(x = Site, y=DO_Conc, fill = Site)) +
  theme_classic(base_size = 22)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  #stat_summary(geom = 'text', label = c("a", "b","c","a", "e","f","g","h"), fun = max, vjust = -0.5, size=6)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  geom_violinhalf()+
  stat_summary(fun = mean, geom = "point", size = 4) + 
  #geom_violin()+
  #geom_boxplot(alpha = 0)+
  #geom_jitter(width = 0.2, aldoa = 0.8)+
  ylab("do")+
  facet_wrap(Timepoint~.)
do.mean
ggsave(do.mean,file="do.mean.box.pdf",h=6,w=12)
ggsave(do.mean,file="do.mean.violinhalf.pdf",h=6,w=12)

#sal
sal.mean <- sal_data %>% ggplot(aes(x = Site, y=Salinity, fill = Site)) +
  theme_classic(base_size = 22)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  #stat_summary(geom = 'text', label = c("a", "b","c","a", "e","f","g","h"), fun = max, vjust = -0.5, size=6)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  geom_violinhalf()+
  stat_summary(fun = mean, geom = "point", size = 4) + 
  #geom_violin()+
  #geom_boxplot(alpha = 0)+
  #geom_jitter(width = 0.2, alsala = 0.8)+
  ylab("sal")+
  facet_wrap(Timepoint~.)+
  scale_y_continuous(limits = c(27, 38), n.breaks = 10)
sal.mean
ggsave(sal.mean,file="sal.mean.box.pdf",h=6,w=12)
ggsave(sal.mean,file="sal.mean.violinhalf.pdf",h=6,w=12)

#par average
par.mean <- par_data %>% ggplot(aes(x = Site, y=PAR_cal, fill = Site)) +
  theme_classic(base_size = 22)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  #stat_summary(geom = 'text', label = c("a", "b","c","a", "e","f","g","h"), fun = max, vjust = -0.5, size=6)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  geom_violinhalf()+
  stat_summary(fun = mean, geom = "point", size = 4) + 
  #geom_violin()+
  #geom_boxplot(alpha = 0)+
  #geom_jitter(width = 0.2, alpara = 0.8)+
  ylab("par")+
  facet_wrap(Timepoint~.)
par.mean
ggsave(par.mean,file="par.mean.box.pdf",h=6,w=12)
ggsave(par.mean,file="par.mean.violinhalf.pdf",h=6,w=12)

#temp using pH data
temp.mean <- pH_data %>% ggplot(aes(x = Site, y=Temp_C, fill = Site)) +
  theme_classic(base_size = 22)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  #stat_summary(geom = 'text', label = c("a", "b","c","a", "e","f","g","h"), fun = max, vjust = -0.5, size=6)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  geom_violinhalf()+
  stat_summary(fun = mean, geom = "point", size = 4) + 
  #geom_violin()+
  #geom_boxplot(alpha = 0)+
  #geom_jitter(width = 0.2, altempa = 0.8)+
  ylab("temp")+
  facet_wrap(Timepoint~.)+
  ylim(26,31)
temp.mean
ggsave(temp.mean,file="temp.mean.box.pH.pdf",h=6,w=12)
ggsave(temp.mean,file="temp.mean.violinhalf.pH.pdf",h=6,w=12)

#temp using DO data
#NAs in the data are making things weird
any(is.na(do_data$Temp))
do_nona <- do_data %>% drop_na(Temp)
#remove the outlier as well
do_nona <- do_nona %>% subset(Temp < 34.00)
#do_nona <- do_nona[-c(30393),] 

temp.mean <- do_nona %>% ggplot(aes(x = Site, y=Temp, fill = Site)) +
  theme_classic(base_size = 22)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  #stat_summary(geom = 'text', label = c("a", "b","c","a", "e","f","g","h"), fun = max, vjust = -0.5, size=6)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  geom_violinhalf()+
  stat_summary(fun = mean, geom = "point", size = 4) + 
  #geom_violin()+
  #geom_boxplot(alpha = 0)+
  #geom_jitter(width = 0.2, altempa = 0.8)+
  ylab("temp")+
  facet_wrap(Timepoint~.)+
  ylim(26, 31)
temp.mean
ggsave(temp.mean,file="temp.mean.box.do.pdf",h=6,w=12)
ggsave(temp.mean,file="temp.mean.violinhalf.do.pdf",h=6,w=12)

####STATS ####

pH_data <-read.csv("all_pH_CW2020.csv")
do_data <-read.csv("all_DO_CW2020.csv")
sal_data <-read.csv("all_sal_CW2020.csv")
par_data <-read.csv("all_PAR_CW2020.csv")
temp_data <-read.csv("all_DO_CW2020.csv")

#running Kruskal-Wallis Tests
#pH
#summary stats
#make functions to calculate standard error and range
se <- function(x) sqrt(var(x)/length(x))
ra <- function(x) max(x)-min(x)

#variables across sites that we want to test are:
# Habitat - Habitat
# Timepoint - Timepoint
# Habitat:Location - Site
# Habitat:Timepoint - Habitat_Time
# Habitat:Location:Timepoint - Site_Time

#stats table of very minimal stats to share:
temp_data %>% drop_na(Temp) %>%
  group_by(Site_Time) %>% 
  summarize(Mean = mean(Temp),
            Min = min(Temp),
            Max = max(Temp),
            Range = ra(Temp),
            Days = length(unique(Date)))
#decided not to include this info in the paper because it's not as relevant/not cited, whereas stats on the differences are cited
#so including stats - see below

library(FSA) #for dunn tests - here doing BH p-value correction for multiple comparisons - stronger than bonferroni
#just doing K-W tests for general comparisons and collapsing factors

#pH
#create separate data frames for each site
#then compute daily averages
pH_SWB <- pH_data %>% filter(Site == "Spaanse Water Bay")
pH_SWR <- pH_data %>% filter(Site == "Spaanse Water Reef")
pH_SMB <- pH_data %>% filter(Site == "Santa Martha Bay")
pH_SMR <- pH_data %>% filter(Site == "Santa Martha Reef")

#filter rows of a dataframe to get only rows with unique values for one variable in r
pH_SWB_da <- distinct(pH_SWB, Date, .keep_all = TRUE)
da_pH_SWB <- pH_SWB %>%
  group_by(Date) %>%
  summarize(mean_daily_pH = mean(pHT))
pH_SWB_da <- merge(pH_SWB_da, da_pH_SWB, by = "Date", all=TRUE)

pH_SWR_da <- distinct(pH_SWR, Date, .keep_all = TRUE)
da_pH_SWR <- pH_SWR %>%
  group_by(Date) %>%
  summarize(mean_daily_pH = mean(pHT))
pH_SWR_da <- merge(pH_SWR_da, da_pH_SWR, by = "Date", all=TRUE)

pH_SMB_da <- distinct(pH_SMB, Date, .keep_all = TRUE)
da_pH_SMB <- pH_SMB %>%
  group_by(Date) %>%
  summarize(mean_daily_pH = mean(pHT))
pH_SMB_da <- merge(pH_SMB_da, da_pH_SMB, by = "Date", all=TRUE)

pH_SMR_da <- distinct(pH_SMR, Date, .keep_all = TRUE)
da_pH_SMR <- pH_SMR %>%
  group_by(Date) %>%
  summarize(mean_daily_pH = mean(pHT))
pH_SMR_da <- merge(pH_SMR_da, da_pH_SMR, by = "Date", all=TRUE)

pH_da <- rbind(pH_SWB_da,pH_SWR_da,pH_SMB_da,pH_SMR_da)

scheirerRayHare(mean_daily_pH ~ Habitat+Timepoint,data = pH_da)
dunnTest(mean_daily_pH ~ Timepoint, data = pH_da, method="bh")

scheirerRayHare(mean_daily_pH ~ Site+Timepoint,data = pH_da)
dunnTest(mean_daily_pH ~ Site, data = pH_da, method="bh")
dunnTest(mean_daily_pH ~ Site_Time, data = pH_da, method="bh")

##Temperature
#create separate data frames for each site
#then compute daily averages
temp_data$Habitat_Time <- paste(temp_data$Habitat, temp_data$Timepoint)
temp_data$Site_Time <- paste(temp_data$Site, temp_data$Timepoint)
temp_SWB <- temp_data %>% filter(Site == "Spaanse Water Bay")
temp_SWR <- temp_data %>% filter(Site == "Spaanse Water Reef")
temp_SMB <- temp_data %>% filter(Site == "Santa Martha Bay")
temp_SMR <- temp_data %>% filter(Site == "Santa Martha Reef")

#filter rows of a dataframe to get only rows with unique values for one variable in r
temp_SWB_da <- distinct(temp_SWB, Date, .keep_all = TRUE)
da_temp_SWB <- temp_SWB %>%
  group_by(Date) %>%
  summarize(mean_daily_temp = mean(Temp))
temp_SWB_da <- merge(temp_SWB_da, da_temp_SWB, by = "Date", all=TRUE)

temp_SWR_da <- distinct(temp_SWR, Date, .keep_all = TRUE)
da_temp_SWR <- temp_SWR %>%
  group_by(Date) %>%
  summarize(mean_daily_temp = mean(Temp))
temp_SWR_da <- merge(temp_SWR_da, da_temp_SWR, by = "Date", all=TRUE)

temp_SMB_da <- distinct(temp_SMB, Date, .keep_all = TRUE)
da_temp_SMB <- temp_SMB %>%
  group_by(Date) %>%
  summarize(mean_daily_temp = mean(Temp))
temp_SMB_da <- merge(temp_SMB_da, da_temp_SMB, by = "Date", all=TRUE)

temp_SMR_da <- distinct(temp_SMR, Date, .keep_all = TRUE)
da_temp_SMR <- temp_SMR %>%
  group_by(Date) %>%
  summarize(mean_daily_temp = mean(Temp))
temp_SMR_da <- merge(temp_SMR_da, da_temp_SMR, by = "Date", all=TRUE)

temp_da <- rbind(temp_SWB_da,temp_SWR_da,temp_SMB_da,temp_SMR_da)

scheirerRayHare(mean_daily_temp ~ Habitat+Timepoint,data = temp_da)
dunnTest(mean_daily_temp ~ Timepoint, data = temp_da, method="bh")

scheirerRayHare(mean_daily_temp ~ Site+Timepoint,data = temp_da)
dunnTest(mean_daily_temp ~ Site, data = temp_da, method="bh")

##DO dissolved oxygen
#create separate data frames for each site
#then compute daily averages
do_data$Habitat_Time <- paste(do_data$Habitat, do_data$Timepoint)
do_data$Site_Time <- paste(do_data$Site, do_data$Timepoint)
do_SWB <- do_data %>% filter(Site == "Spaanse Water Bay")
do_SWR <- do_data %>% filter(Site == "Spaanse Water Reef")
do_SMB <- do_data %>% filter(Site == "Santa Martha Bay")
do_SMR <- do_data %>% filter(Site == "Santa Martha Reef")

#filter rows of a dataframe to get only rows with unique values for one variable in r
do_SWB_da <- distinct(do_SWB, Date, .keep_all = TRUE)
da_do_SWB <- do_SWB %>%
  group_by(Date) %>%
  summarize(mean_daily_do = mean(DO_Conc))
do_SWB_da <- merge(do_SWB_da, da_do_SWB, by = "Date", all=TRUE)

do_SWR_da <- distinct(do_SWR, Date, .keep_all = TRUE)
da_do_SWR <- do_SWR %>%
  group_by(Date) %>%
  summarize(mean_daily_do = mean(DO_Conc))
do_SWR_da <- merge(do_SWR_da, da_do_SWR, by = "Date", all=TRUE)

do_SMB_da <- distinct(do_SMB, Date, .keep_all = TRUE)
da_do_SMB <- do_SMB %>%
  group_by(Date) %>%
  summarize(mean_daily_do = mean(DO_Conc))
do_SMB_da <- merge(do_SMB_da, da_do_SMB, by = "Date", all=TRUE)

do_SMR_da <- distinct(do_SMR, Date, .keep_all = TRUE)
da_do_SMR <- do_SMR %>%
  group_by(Date) %>%
  summarize(mean_daily_do = mean(DO_Conc))
do_SMR_da <- merge(do_SMR_da, da_do_SMR, by = "Date", all=TRUE)

do_da <- rbind(do_SWB_da,do_SWR_da,do_SMB_da,do_SMR_da)

scheirerRayHare(mean_daily_do ~ Habitat+Timepoint,data = do_da)
dunnTest(mean_daily_do ~ Timepoint, data = do_da, method="bh")
dunnTest(mean_daily_do ~ Habitat_Time, data = do_da, method="bh")

scheirerRayHare(mean_daily_do ~ Site+Timepoint,data = do_da)
dunnTest(mean_daily_do ~ Site, data = do_da, method="bh")
dunnTest(mean_daily_do ~ Site_Time, data = do_da, method="bh")

##salinty
#create separate data frames for each site
#then compute daily averages
sal_data$Habitat_Time <- paste(sal_data$Habitat, sal_data$Timepoint)
sal_data$Site_Time <- paste(sal_data$Site, sal_data$Timepoint)
sal_SWB <- sal_data %>% filter(Site == "Spaanse Water Bay")
sal_SWR <- sal_data %>% filter(Site == "Spaanse Water Reef")
sal_SMB <- sal_data %>% filter(Site == "Santa Martha Bay")
sal_SMR <- sal_data %>% filter(Site == "Santa Martha Reef")

#filter rows of a dataframe to get only rows with unique values for one variable in r
sal_SWB_da <- distinct(sal_SWB, Date, .keep_all = TRUE)
da_sal_SWB <- sal_SWB %>%
  group_by(Date) %>%
  summarize(mean_daily_sal = mean(Salinity))
sal_SWB_da <- merge(sal_SWB_da, da_sal_SWB, by = "Date", all=TRUE)

sal_SWR_da <- distinct(sal_SWR, Date, .keep_all = TRUE)
da_sal_SWR <- sal_SWR %>%
  group_by(Date) %>%
  summarize(mean_daily_sal = mean(Salinity))
sal_SWR_da <- merge(sal_SWR_da, da_sal_SWR, by = "Date", all=TRUE)

sal_SMB_da <- distinct(sal_SMB, Date, .keep_all = TRUE)
da_sal_SMB <- sal_SMB %>%
  group_by(Date) %>%
  summarize(mean_daily_sal = mean(Salinity))
sal_SMB_da <- merge(sal_SMB_da, da_sal_SMB, by = "Date", all=TRUE)

sal_SMR_da <- distinct(sal_SMR, Date, .keep_all = TRUE)
da_sal_SMR <- sal_SMR %>%
  group_by(Date) %>%
  summarize(mean_daily_sal = mean(Salinity))
sal_SMR_da <- merge(sal_SMR_da, da_sal_SMR, by = "Date", all=TRUE)

sal_da <- rbind(sal_SWB_da,sal_SWR_da,sal_SMB_da,sal_SMR_da)

scheirerRayHare(mean_daily_sal ~ Habitat+Timepoint,data = sal_da)
dunnTest(mean_daily_sal ~ Timepoint, data = sal_da, method="bh")
dunnTest(mean_daily_sal ~ Habitat_Time, data = sal_da, method="bh")

scheirerRayHare(mean_daily_sal ~ Site+Timepoint,data = sal_da)
dunnTest(mean_daily_sal ~ Site, data = sal_da, method="bh")
dunnTest(mean_daily_sal ~ Site_Time, data = sal_da, method="bh")

##parerature
#create separate data frames for each site
#then compute daily averages
par_data$Habitat_Time <- paste(par_data$Habitat, par_data$Timepoint)
par_data$Site_Time <- paste(par_data$Site, par_data$Timepoint)
par_SWB <- par_data %>% filter(Site == "Spaanse Water Bay")
par_SWR <- par_data %>% filter(Site == "Spaanse Water Reef")
par_SMB <- par_data %>% filter(Site == "Santa Martha Bay")
par_SMR <- par_data %>% filter(Site == "Santa Martha Reef")

#filter rows of a dataframe to get only rows with unique values for one variable in r
par_SWB_da <- distinct(par_SWB, Date, .keep_all = TRUE)
da_par_SWB <- par_SWB %>%
  group_by(Date) %>%
  summarize(mean_daily_par = mean(PAR_cal))
par_SWB_da <- merge(par_SWB_da, da_par_SWB, by = "Date", all=TRUE)

par_SWR_da <- distinct(par_SWR, Date, .keep_all = TRUE)
da_par_SWR <- par_SWR %>%
  group_by(Date) %>%
  summarize(mean_daily_par = mean(PAR_cal))
par_SWR_da <- merge(par_SWR_da, da_par_SWR, by = "Date", all=TRUE)

par_SMB_da <- distinct(par_SMB, Date, .keep_all = TRUE)
da_par_SMB <- par_SMB %>%
  group_by(Date) %>%
  summarize(mean_daily_par = mean(PAR_cal))
par_SMB_da <- merge(par_SMB_da, da_par_SMB, by = "Date", all=TRUE)

par_SMR_da <- distinct(par_SMR, Date, .keep_all = TRUE)
da_par_SMR <- par_SMR %>%
  group_by(Date) %>%
  summarize(mean_daily_par = mean(PAR_cal))
par_SMR_da <- merge(par_SMR_da, da_par_SMR, by = "Date", all=TRUE)

par_da <- rbind(par_SWB_da,par_SWR_da,par_SMB_da,par_SMR_da)

scheirerRayHare(mean_daily_par ~ Habitat+Timepoint,data = par_da)
dunnTest(mean_daily_par ~ Timepoint, data = par_da, method="bh")

scheirerRayHare(mean_daily_par ~ Site+Timepoint,data = par_da)
dunnTest(mean_daily_par ~ Site, data = par_da, method="bh")

#EXCEL CODE FOR MAKING p-values into significance values
#=IF(A2<0.001, "***", IF(A2<=0.01, "**", IF(A2<=0.05, "*", "ns")))
