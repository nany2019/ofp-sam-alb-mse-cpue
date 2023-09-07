library(VAST)
library(tidyverse)
library(magrittr)
library(data.table)
library(splines)

# Remove all previous objects
rm(list=ls())

setwd('C:/Users/nany/Documents/')
working_dir = paste0(getwd(),'/vast_cpp/')

# Read in filtered data - from filter_data.r
load('ops.model_2023.Rdata')

source('./alb_cpue_prep.R')
source('./alb_add_regions.R')

# Correct logdate month, quarter - as format changed for original file
ops.filt %<>% mutate(month = as.numeric(substr(logdate,4,5)), quarter=(ceiling(month/3)/4)-0.25, YrQtr = year+quarter )
ops.filt %<>%  mutate(drop = ifelse(Region %in% c(4) & flag_id!='JP',1,0)) %>% filter(drop==0) %>% select(-drop)

##filter for core fleet, can't do core fleet, but can try flag, at least 20 set for 4 quarters in a year 
keep <- ops.filt %>% group_by(flag_id,YrQtr,Region) %>% summarise(freq = n()) %>% ungroup() %>%
  filter(freq >= 20) %>% select(-freq) %>% mutate(keep = 1)
ops.filt %<>% left_join(keep, by=c("flag_id", "YrQtr","Region")) %>% filter(!is.na(keep))


##############################################################################
# Settings
knots = 150
samp.rate = 0.05
yr_range = c(1960,2021)
grid.deg = 2

# Covariates
q.vars = c('targ.spec')
Q1_formula = ~ factor(targ.spec) + factor(flag_id) 
Q2_formula = ~  factor(targ.spec) +factor(flag_id) 
Overdispersion = c('Eta1'=0, 'Eta2'=0)


# Subsampling routine 
nsamps = round(nrow(ops.filt)*samp.rate) # Number of samples
nn = length(unique(ops.filt$cell.5x5)) # How many 5x5 cells actually exist 
ny = length(unique(ops.filt$YrQtr))
N =ceiling(nsamps/(nn*ny))*2 # Number of samples to draw, evenly across knots (or max nobs per knot if < N)
setDT(ops.filt) 

seed=123
set.seed(seed) 
sub = ops.filt[year>=yr_range[1] & year<=yr_range[2],.SD[sample(.N, min(N,.N))], by = list(cell.5x5, YrQtr)] %>% data.frame() 

rm(ops.filt)
#########################################################################################################################
pdat = prep.alb(sub, knots, grid.deg) # Prep data while defining knots

#########################################################################################################################

DG = pdat[[1]]; Region = pdat[[2]]; strata.limits = pdat[[3]]; grid.df = pdat[[4]]

# Select target cluster
DG$targ.spec = as.factor(DG$k3.h)

## Remove combinations of levels only sampled once,
## otherwise Hessian non-positive (e.g. knot x yrqtr)
count.kyq <- with(DG, table(knot_i, YrQtr))
knot.1yq <- rownames(count.kyq)[rowSums(count.kyq>=1)<2] ## get knots that only occur in one yrqtr
nrk <- nrow(DG)
DG %<>% filter(!(knot_i %in% as.numeric(knot.1yq)))
message(sprintf('%s knot(s) only sampled in three quarters, removed %s observations',
                length(knot.1yq), nrk-nrow(DG)))

## Do this for lognormal component too
count.kyq <- with(filter(DG, alb_cpue>0), table(knot_i, YrQtr))
knot.1yq <- rownames(count.kyq)[rowSums(count.kyq>=1)<2] ## get knots that only occur in one yrqtr
nrk <- nrow(DG)
DG %<>% filter(!(knot_i %in% as.numeric(knot.1yq)))
message(sprintf('%s knot(s) only sampled in three quarters, removed %s observations',
                length(knot.1yq), nrk-nrow(DG)))

# Assign effort offset and rename
DG.fit = DG %>% data.frame() %>% mutate(Area=1) # effort offset
rm(DG)
DG.fit %<>% filter(!is.na(Region))
# Create extrapolation grid
Extrapolation_List = make_extrapolation_info(Region=Region, strata.limits=strata.limits, input_grid=grid.df, 
                                             projargs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Identify region associated with each observation - assign targeting cluster accordingly
regD = add.region(DG.fit, Extrapolation_List, strata.limits, grid.df)
DG.fit = regD[[1]]; Extrapolation_List = regD[[2]]
rm(regD)




Spatial_List = make_spatial_info(grid_size_km=grid.deg*110, n_x=knots,
                                 Lon_i=DG.fit[,'lond'], Lat_i=DG.fit[,'latd'], 
                                 Extrapolation_List=Extrapolation_List, 
                                 DirPath=working_dir, Save_Results=FALSE,
                                 fine_scale=F, knot_method='grid')

##data plot 
 # plot_data(Extrapolation_List = Extrapolation_List,Spatial_List = Spatial_List, Data_Geostat = DG.fit, Lat_i = DG.fit[, "latd"],
 #           Lon_i = DG.fit[, "lond"],
 #           Year_i = DG.fit[, "year"],)
# DG.fit %>% 
#   filter(!is.na(Region)) %>%
#   group_by(year, Region) %>%
#   summarise(Nominal_CPUE=mean(alb_cpue)) %>%
#   ggplot(aes(x=year,y=Nominal_CPUE)) +geom_point()+geom_line() + geom_vline(xintercept = 2018, linetype="dashed") +facet_grid(Region~.)+theme_bw() +
#   theme(legend.title = element_text(size=16),legend.text = element_text(size=16),axis.title =element_text(size=20), axis.text = element_text(size=16),strip.text.y = element_text(size = 16),legend.position = "bottom")


#............................ 
# Model settings 
#............................ 
Version = 'VAST_v12_0_0' #get_latest_version( package="VAST" ) 

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) 
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=2, "Epsilon2"=2) 

ObsModel = c(2,3) # Delta-lognormal c(4,3) or(4,4) - Delta-Gamma c(2,3) or(2,4) (we usually use 2,4)
Options =  c("SD_site_density"=0, "SD_site_logdensity"=0, "Calculate_Range"=0, "Calculate_evenness"=0, 
             "Calculate_effective_area"=0, "Calculate_Cov_SE"=0, 'Calculate_Synchrony'=0, 'Calculate_Coherence'=0) 

# Renumber years to match covariate indexing 
Year_Set <- sort(unique(DG.fit$year))
Years2Include <- 1:length(Year_Set)
YearIndex <- Years2Include
names(YearIndex) <- as.character(Year_Set)

## Add separate YQ effect
YQ_Set <- sort(unique(DG.fit$YrQtr))
YQ_Index <- 1:length(YQ_Set)
names(YQ_Index) <- as.character(YQ_Set)


yrs = seq(min(DG.fit$YrQtr),max(DG.fit$YrQtr),.25) 
t_i = match(DG.fit[,'YrQtr'],  yrs ) 

# Note - as of 2021 there was a default of 2000 cells in the extrap grid
#     max_cells argument adjusts that, as we have more in the South Pacific
# Ideally, we'd use fine_scale=T, but the model was running very slow on the VM, so voila...

settings = make_settings(n_x=knots, Region='User', purpose='index2', max_cells=Inf,
                         strata.limits=strata.limits,bias.correct=F, 
                         fine_scale=F, Version=Version) 

#............................ 
# Now fit Model 
#............................ 
A = proc.time()
sink("vast_iters.txt") # suppress output 
fit = fit_model.tvc( settings=settings, 
                 Lat_i=DG.fit[,'latd'], 
                 Lon_i=DG.fit[,'lond'], 
                 t_i=t_i, 
                 c_iz = rep(0,nrow(DG.fit)), 
                 b_i=DG.fit[,'alb_cpue'], 
                 a_i=DG.fit[,'Area'], 
                 OverdispersionConfig = Overdispersion,
                 catchability_data = DG.fit,
                 Q1_formula = Q1_formula,
                 Q2_formula = Q2_formula,
                 input_grid=grid.df, 
                 FieldConfig=FieldConfig, 
                 RhoConfig=RhoConfig, 
                 ObsModel_ez=ObsModel, 
                 extrapolation_list = Extrapolation_List,
                 spatial_list = Spatial_List,
                 CheckForErrors=T,
                 getReportCovariance=T,
                 Aniso=T,
                 working_dir = working_dir) 
sink() 
B=proc.time()
run.time = (B[[3]]-A[[3]])/3600

# Dharma residuals
## rm(dharmaRes)
## dharmaRes = summary(fit, what="residuals", working_dir=NA , n_samples=1500)

out = list(fit, DG.fit)  ##, dharmaRes, run.time)

save(out, file=paste0('P:/SLOTH/alb_MSE/Geo_CPUE/2023update_onlyJP_raneff.Rdata'))

# DG.fit %>% 
#   group_by(year, Region) %>%
#   summarise(Nominal_CPUE=mean(alb_cpue)) %>%
#   ggplot(aes(x=year,y=Nominal_CPUE)) +geom_point()+geom_line() + geom_vline(xintercept = 2018, linetype="dashed") +facet_grid(Region~.)+theme_bw() +
#   theme(legend.title = element_text(size=16),legend.text = element_text(size=16),axis.title =element_text(size=20), axis.text = element_text(size=16),strip.text.y = element_text(size = 16),legend.position = "bottom")
# 

