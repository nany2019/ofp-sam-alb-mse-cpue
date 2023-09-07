# Fill in and test HBF for full time series and split time series
library(randomForest)
library(tidyverse)
library(magrittr)
library(data.table)

# Plotting theme
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


confMat = function(pred.output){
  pred.output %<>% group_by(hbf,pred) %>% summarise(n = n()) %>% spread(key=hbf, value=n, fill=0) %>% data.frame()
  accuracy = t(t(pred.output[,-1])/colSums(pred.output[,-1]))
  row.names(accuracy) = pred.output[,1]
  
  overall.accuracy = diag(as.matrix(pred.output[,-1]))/as.numeric(colSums(pred.output[,-1]))
  names(overall.accuracy) = colnames(pred.output)[-1]
  
  return(list(pred.output, accuracy, overall.accuracy))
}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

################################################
# Load data
setwd("C:/Users/nany/Documents/")


# Define time series
ts = c('full')#, 'early','recent')
load('ops.filt_2023.Rdata')

##################################################
# Predict missing HBF values 
##################################################

  
  
     tmp=ops.filt

  
  
  tmp %<>% mutate(hbf = ifelse(hk_bt_flt==0,NA,hk_bt_flt)) %>% ungroup()
  tmp$hbf =as.numeric(tmp$hbf)
  tmp$hbf = factor(ceiling(tmp$hbf/5)*5, levels= sort(unique(ceiling(tmp$hbf/5)*5)) )
  tmp$hbf.label = as.numeric(tmp$hbf)-1
  
  tmp$total.catch = tmp$alb_n + tmp$bet_n + tmp$yft_n + tmp$swo_n 
  tmp %<>% mutate(alb =  alb_n/total.catch,
                  bet = bet_n/total.catch,
                  yft = yft_n/total.catch,
                  swo = swo_n/total.catch)
  
# Now fill in missing values 

hbf.good = tmp %>% filter(!is.na(hbf)) %>% select(hbf,year,flag_id, month, lond, latd, cell.5x5, hhook,total.catch, alb, bet, yft, swo)
setDT(hbf.good)

hbf.miss = tmp %>% filter(is.na(hbf)) %>% select(hbf,year,flag_id, month, lond, latd, cell.5x5, hhook,total.catch, alb, bet, yft, swo)

##look at the missing part of the data set
flag_good <- as.data.frame(table(hbf.good$flag_id))
colnames(flag_good)<-c("Flag","hbf.good")

flag_miss <- as.data.frame(table(hbf.miss$flag_id))
colnames(flag_miss)<-c("Flag","hbf.miss")


table_ratio <- full_join(flag_good,flag_miss) %>%
rowwise() %>%
mutate(total=sum(hbf.good,hbf.miss,na.rm=TRUE)) %>%
mutate(miss_ratio=hbf.miss/total) 

#png("missing_ratio.png",width = 1200, height = 1000)
#print(full_join(flag_good,flag_miss) %>%
#rowwise() %>% pivot_longer(-Flag) %>%
#ggplot(aes(fill=name, y=value, x=Flag)) + 
#    geom_bar(position="fill", stat="identity"))
#dev.off()


# Repeat with 10 bootstrap samples
# define subsampling

nsamp = (nrow(hbf.good)*0.2) # about 20% of all samples
nstrata = length(unique(tmp$year))*length(unique(tmp$cell.5x5))
N = round(nsamp/nstrata)

pred.hbf= list()

for(j in 1:10){
  
  train= hbf.good[,.SD[sample(.N, min(N,.N))],by = list(year,cell.5x5)] %>% data.frame() %>% select(-cell.5x5)

  rf = randomForest(hbf~., data=train )
  rf.pred = predict(rf, hbf.miss, predict.all=TRUE)
  
  pred.hbf[[j]] = hbf.miss %>% mutate(pred = rf.pred$aggregate)
}

##########################################
# Get estimates of HBF and compare
##########################################
hbf.estimates = matrix()

for (j in 1:10){
  hbf.estimates = cbind(hbf.estimates, pred.hbf[[j]][,'pred']  )
}
hbf.estimates = hbf.estimates[,-1] %>% data.frame()
hbf.estimates %<>% mutate(row = 1:nrow(hbf.estimates))

hbf.estimates$mode = as.numeric(apply(hbf.estimates[,-11], 1,mode))

rm(hbf.good, hbf.miss)
##########################################
# Fill in data set
##########################################
tmp %<>% mutate(hbf = ifelse(hk_bt_flt==0,NA,hk_bt_flt)) %>% ungroup()
tmp$hbf = as.numeric(tmp$hbf)
tmp$hbf = factor(ceiling(tmp$hbf/5)*5, levels= sort(unique(ceiling(tmp$hbf/5)*5)) )

fill.hbf = tmp %>% filter(is.na(hbf)) %>% mutate(hbf = hbf.estimates$mode, hbf.fill=1)
good.hbf =  tmp %>% filter(!is.na(hbf)) %>% mutate(hbf.fill=0)

ops.filt.hbf = rbind(fill.hbf, good.hbf) %>% mutate(hbf = as.numeric(hbf))

save(ops.filt.hbf, file=paste0('ops.filt.hbf_2023.Rdata'))
##}

