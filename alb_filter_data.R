library(data.table)
library(tidyverse)
library(magrittr)
library(sp)

gg.theme <- theme_bw() + 
  theme(axis.line = element_line(color="black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"))  

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

############################################
setwd('C:/Users/nany/Documents/')
############################################

# Load data from: prep.ops.data.R
load("ops.trim.dt_2021_full.Rdata")
names(ops.trim.dt)

# Filtering is based off the work done in 2018 
#   Repeating for 2021 analysis

ops.filt = ops.trim.dt[latd>=-50 & latd<=0] # Spatial filter 0 - -50S
ops.filt = ops.filt[lond>=130 ] # Spatial filter > 140  Last assessment: 140 -230
#ops.filt %<>% mutate(drop.ll = ifelse(lond>210 & latd>-5,1,0)) %>% filter(drop.ll==0) %>% select(-drop.ll)

ops.filt = ops.filt[!is.na(logdate)] # No missing logdate
ops.filt = ops.filt[hook!=0] # Hooks >0
ops.filt = ops.filt[year>=1960] # Year >= 1960
ops.filt = ops.filt[alb_cpue>0 | bet_cpue>0 | yft_cpue>0 | swo_cpue>0] #At least 1 fish of one of the 4 main species 


# Check on missing vessel id
#source('./RCode/check.missing.vessel.R')

# Identify strata with insufficient data
# Must have at least 50 sets per year.qtr and 20 sets per 5x5 cell
ops.filt %<>% group_by(yr.qtr) %>% mutate(drop.yr.qtr = ifelse(n()<50,1,0)) %>% group_by(cell.5x5) %>% mutate(drop.cell.5x5 = ifelse(n()<20,1,0))
ops.filt %<>% filter(drop.cell.5x5==0 & drop.yr.qtr==0) # Keep only yr-qtrs with 50+ sets; 5x5 cells with 20+ sets 
ops.filt %<>% select(-c(drop.yr.qtr, drop.cell.5x5))

ops.filt %<>% mutate(YrQtr = year + quarter/4 - 0.25)
ops.filt %<>% mutate(date = as.Date(logdate, format='%m/%d/%Y'))

# Retain only JP in the EPO
load("E:/2021_SC/CPUE/ALB/Data/alb_regions_poly_2021.Rdata")

ops.filt$Region = NA
##change to 1:4
for(i in 1:4){
  coords = regions[[i]]@polygons[[1]]@Polygons[[1]]@coords
  ops.filt$Region = ifelse(point.in.polygon(ops.filt$lond, ops.filt$latd, coords[,1], coords[,2]) %in% c(1,2),i,ops.filt$Region)
}
##only use 4
##skip this step to use all the flag in the analysis 
##ops.filt %<>% mutate(drop = ifelse(Region %in% c(4) & flag_id!='JP',1,0)) %>% filter(drop==0) %>% select(-drop)

rm(coords, ops.trim.dt, regions)

save(ops.filt, file='ops.filt_2021_full.Rdata')

