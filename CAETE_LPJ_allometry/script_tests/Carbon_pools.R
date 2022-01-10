#library(dplyr)
#library(reshape2)
#library(data.table)
library(tidyr)
library (ggplot2)

CarbonPools <- read.csv("carbon_pools_time_5PLS.csv", header = TRUE)

### LEAF CARBON #############

Cleaf <- CarbonPools[,c(1,4,5)]

Cleaf_by_pls <- data.frame(pivot_wider(Cleaf, names_from = pls, values_from = leaf)) 
colnames(Cleaf_by_pls) <- c("time", "pls1", "pls2", "pls3", "pls4", "pls5")

#plot Cleaf of all pls
plotCleaf <- ggplot(Cleaf_by_pls, aes(time)) +  
  geom_line(aes(y = pls1), color = "purple") +
  geom_line(aes(y = pls2), color = "red") +
  geom_line(aes(y = pls3), color = "grey") +
  geom_line(aes(y = pls4), color = "blue") +
  geom_line(aes(y = pls5), color = "pink") +
  labs(y= "Cleaf", x = "Time")
plotCleaf

#calculate mean
Cleaf_by_pls$mean <- rowMeans(Cleaf_by_pls[,2:6])

#Adding Min and Max and some index for the x axis
#Where, in "apply", 1 means to apply FUN to each row of dataframef, 2 would mean to apply FUN to columns.
Cleaf_by_pls <- transform(Cleaf_by_pls, Min = apply(Cleaf_by_pls[,2:6], 1, FUN = min), Max = apply(Cleaf_by_pls[,2:6], 1, FUN = max), indx = seq_len(dim(Cleaf_by_pls)[1]))


#plot mean and amplitude of data 
plot_cleaf_mean <- ggplot(Cleaf_by_pls) +
  geom_line(aes(indx, mean), group = 1) +
  geom_ribbon(aes(x = indx, ymax = Max, ymin = Min), 
              alpha = 0.6, 
              #linetype=1,      #solid, dashed or other line types
              #colour="green",  #border line color
              #size=1,          #border line size 
              fill = "green")+
  labs(y= "Cleaf", x = "Time")
plot_cleaf_mean

#combine two plots
two_plots_cleaf <- ggplot(Cleaf_by_pls) +
  geom_line(aes(x=time, y = pls1), color = "purple") +
  geom_line(aes(x=time, y = pls2), color = "red") +
  geom_line(aes(x=time, y = pls3), color = "grey") +
  geom_line(aes(x=time, y = pls4), color = "blue") +
  geom_line(aes(x=time, y = pls3), color = "pink") +
  geom_line(aes(x=time, y = mean), group = 1) +
  geom_ribbon(aes(x = time, ymax = Max, ymin = Min), 
              alpha = 0.6,     #transparency
              #linetype=1,      #solid, dashed or other line types
              #colour="green",  #border line color
              #size=1,          #border line size
              fill = "green") +
  labs(y= "Cleaf", x = "Time") 
two_plots_cleaf

##### WOOD CARBON ########

Cwood <- CarbonPools[,c(2,4,5)]

Cwood_by_pls <- data.frame(pivot_wider(Cwood, names_from = pls, values_from = wood)) 
colnames(Cwood_by_pls) <- c("time", "pls1", "pls2", "pls3", "pls4", "pls5")

#plot Cwood of all pls
plotCwood <- ggplot(Cwood_by_pls, aes(time)) +  
  geom_line(aes(y = pls1), color = "purple") +
  geom_line(aes(y = pls2), color = "red") +
  geom_line(aes(y = pls3), color = "grey") +
  geom_line(aes(y = pls4), color = "blue") +
  geom_line(aes(y = pls5), color = "pink") +
  labs(y= "Cwood", x = "Time")
plotCwood

#calculate mean
Cwood_by_pls$mean <- rowMeans(Cwood_by_pls[,2:6])

#Adding Min and Max and some index for the x axis
#Where, in "apply", 1 means to apply FUN to each row of dataframef, 2 would mean to apply FUN to columns.
Cwood_by_pls <- transform(Cwood_by_pls, Min = apply(Cwood_by_pls[,2:6], 1, FUN = min), Max = apply(Cwood_by_pls[,2:6], 1, FUN = max), indx = seq_len(dim(Cwood_by_pls)[1]))


#plot mean and amplitude of data 
plot_cwood_mean <- ggplot(Cwood_by_pls) +
  geom_line(aes(indx, mean), group = 1) +
  geom_ribbon(aes(x = indx, ymax = Max, ymin = Min), 
              alpha = 0.6, 
              #linetype=1,      #solid, dashed or other line types
              #colour="green",  #border line color
              #size=1,          #border line size 
              fill = "brown")+
  labs(y= "Cwood", x = "Time")
plot_cwood_mean

#combine two plots
two_plots_cwood <- ggplot(Cwood_by_pls) +
  geom_line(aes(x=time, y = pls1), color = "purple") +
  geom_line(aes(x=time, y = pls2), color = "red") +
  geom_line(aes(x=time, y = pls3), color = "grey") +
  geom_line(aes(x=time, y = pls4), color = "blue") +
  geom_line(aes(x=time, y = pls3), color = "pink") +
  geom_line(aes(x=time, y = mean), group = 1) +
  geom_ribbon(aes(x = time, ymax = Max, ymin = Min), 
              alpha = 0.6,     #transparency
              #linetype=1,      #solid, dashed or other line types
              #colour="green",  #border line color
              #size=1,          #border line size
              fill = "brown") +
  labs(y= "Cwood", x = "Time") 
two_plots_cwood

#### ROOT CARBON #####

Croot <- CarbonPools[,c(3,4,5)]

Croot_by_pls <- data.frame(pivot_wider(Croot, names_from = pls, values_from = root)) 
colnames(Croot_by_pls) <- c("time", "pls1", "pls2", "pls3", "pls4", "pls5")

#plot Croot of all pls
plotCroot <- ggplot(Croot_by_pls, aes(time)) +  
  geom_line(aes(y = pls1), color = "purple") +
  geom_line(aes(y = pls2), color = "red") +
  geom_line(aes(y = pls3), color = "grey") +
  geom_line(aes(y = pls4), color = "blue") +
  geom_line(aes(y = pls5), color = "pink") +
  labs(y= "Croot", x = "Time")
plotCroot

#calculate mean
Croot_by_pls$mean <- rowMeans(Croot_by_pls[,2:6])

#Adding Min and Max and some index for the x axis
#Where, in "apply", 1 means to apply FUN to each row of dataframef, 2 would mean to apply FUN to columns.
Croot_by_pls <- transform(Croot_by_pls, Min = apply(Croot_by_pls[,2:6], 1, FUN = min), Max = apply(Croot_by_pls[,2:6], 1, FUN = max), indx = seq_len(dim(Croot_by_pls)[1]))


#plot mean and amplitude of data 
plot_croot_mean <- ggplot(Croot_by_pls) +
  geom_line(aes(indx, mean), group = 1) +
  geom_ribbon(aes(x = indx, ymax = Max, ymin = Min), 
              alpha = 0.6, 
              #linetype=1,      #solid, dashed or other line types
              #colour="green",  #border line color
              #size=1,          #border line size 
              fill = "skyblue")+
  labs(y= "Croot", x = "Time")
plot_croot_mean

#combine two plots
two_plots_croot <- ggplot(Croot_by_pls) +
  geom_line(aes(x=time, y = pls1), color = "purple") +
  geom_line(aes(x=time, y = pls2), color = "red") +
  geom_line(aes(x=time, y = pls3), color = "grey") +
  geom_line(aes(x=time, y = pls4), color = "blue") +
  geom_line(aes(x=time, y = pls3), color = "pink") +
  geom_line(aes(x=time, y = mean), group = 1) +
  geom_ribbon(aes(x = time, ymax = Max, ymin = Min), 
              alpha = 0.6,     #transparency
              #linetype=1,      #solid, dashed or other line types
              #colour="green",  #border line color
              #size=1,          #border line size
              fill = "skyblue") +
  labs(y= "Croot", x = "Time") 
two_plots_croot

