library(tidyverse)
library(magrittr)
library(ggplot2)
library(sp)
##compare data 
##compare of operational data 
load('E:/2021_SC/CPUE/ops.trim.dt.RData')
ops.filt_SA <- ops.trim.dt[,c("year","logdate","lond","latd","flag_id","yr.qtr","alb_cpue")]
load("E:/2021_SC/CPUE/ALB/Data/alb_regions_poly_2021.Rdata")
ops.filt_SA$Region = NA
##change to 1:4
for(i in 1:4){
  coords = regions[[i]]@polygons[[1]]@Polygons[[1]]@coords
  ops.filt_SA$Region = ifelse(point.in.polygon(ops.filt_SA$lond, ops.filt_SA$latd, coords[,1], coords[,2]) %in% c(1,2),i,ops.filt_SA$Region)
}
ops.filt_SA %<>% mutate(month = as.numeric(substr(logdate,4,5)), quarter=(ceiling(month/3)/4)-0.25, YrQtr = year+quarter )
ops.filt_SA$cate<-"SA"
##the latest updated 2023 data

load('C:/Users/nany/Documents/ops.trim.dt_2023.Rdata')
ops.filt_new<-ops.trim.dt[,c("year","logdate","lond","latd","flag_id","yr.qtr","alb_cpue")]
ops.filt_new <-subset(ops.filt_new,alb_cpue!="Inf")
ops.filt_new$Region = NA
##change to 1:4
for(i in 1:4){
  coords = regions[[i]]@polygons[[1]]@Polygons[[1]]@coords
  ops.filt_new$Region = ifelse(point.in.polygon(ops.filt_new$lond, ops.filt_new$latd, coords[,1], coords[,2]) %in% c(1,2),i,ops.filt_new$Region)
}

ops.filt_new$cate<- "2023"
ops.filt_new %<>% mutate(month = as.numeric(substr(logdate,4,5)), quarter=(ceiling(month/3)/4)-0.25, YrQtr = year+quarter )


png("P:/SLOTH/alb_MSE/Geo_CPUE/compare_operational_LL_SAvs2023.png", width = 1400, height = 600)
rbind(ops.filt_SA,ops.filt_new) %>%
   filter(Region %in% c(1:4)) %>%
  group_by(YrQtr, Region,cate) %>%
  summarise(Nominal_CPUE=mean(alb_cpue,na.rm = TRUE)) %>%
  ggplot(aes(x=YrQtr,y=Nominal_CPUE,color=cate)) +geom_point()+geom_line() + geom_vline(xintercept = 2018, linetype="dashed") +theme_bw() + facet_grid(Region~.)+
  theme(legend.title = element_text(size=16),legend.text = element_text(size=16),axis.title =element_text(size=20), axis.text = element_text(size=16),strip.text.y = element_text(size = 16),legend.position = "bottom")
dev.off()  











load('E:/2021_SC/CPUE/ALB/Data/ops.model.Rdata')
ops.filt_SA <- ops.filt[,c("year","logdate","lond","latd","flag_id","yr.qtr","alb_cpue")]
load("E:/2021_SC/CPUE/ALB/Data/alb_regions_poly_2021.Rdata")
ops.filt_SA$Region = NA
##change to 1:4
for(i in 1:4){
  coords = regions[[i]]@polygons[[1]]@Polygons[[1]]@coords
  ops.filt_SA$Region = ifelse(point.in.polygon(ops.filt_SA$lond, ops.filt_SA$latd, coords[,1], coords[,2]) %in% c(1,2),i,ops.filt_SA$Region)
}


ops.filt_SA %<>% mutate(month = as.numeric(substr(logdate,4,5)), quarter=(ceiling(month/3)/4)-0.25, YrQtr = year+quarter )
ops.filt_SA$cate<-"SA"
##the latest updated 2023 data

load('C:/Users/nany/Documents/ops.model_2023.Rdata')
ops.filt_new<-ops.filt[,c("year","logdate","lond","latd","Region","flag_id","yr.qtr","alb_cpue")]
ops.filt_new$cate<- "JP_only"
ops.filt_new %<>% mutate(month = as.numeric(substr(logdate,4,5)), quarter=(ceiling(month/3)/4)-0.25, YrQtr = year+quarter )

png("P:/SLOTH/alb_MSE/Geo_CPUE/compare_operational_LL_SAvsupdate.png", width = 1400, height = 600)
rbind(ops.filt_SA,ops.filt_new) %>%
  filter(Region %in% c(1:4),YrQtr<2020) %>%
  group_by(YrQtr, Region,cate) %>%
  summarise(Nominal_CPUE=mean(alb_cpue)) %>%
  ggplot(aes(x=YrQtr,y=Nominal_CPUE,color=cate)) +geom_point()+geom_line() + geom_vline(xintercept = 2018, linetype="dashed") +theme_bw() + facet_grid(Region~.)+
  theme(legend.title = element_text(size=16),legend.text = element_text(size=16),axis.title =element_text(size=20), axis.text = element_text(size=16),strip.text.y = element_text(size = 16),legend.position = "bottom")
dev.off()  
  

##make hyper data
load('C:/Users/nany/Documents/ops.model_2021_full.Rdata')

add.data <- ops.filt %>% filter(Region==4 & year %in% c(2018,2019), flag_id =="JP") 
add.data %<>%
  select(-c(k5,k6,k7,Region))

load('E:/2021_SC/CPUE/ALB/Data/ops.model.Rdata')
load("E:/2021_SC/CPUE/ALB/Data/alb_regions_poly_2021.Rdata")
ops.filt$Region = NA
##change to 1:4
for(i in 1:4){
  coords = regions[[i]]@polygons[[1]]@Polygons[[1]]@coords
  ops.filt$Region = ifelse(point.in.polygon(ops.filt$lond, ops.filt$latd, coords[,1], coords[,2]) %in% c(1,2),i,ops.filt$Region)
}

ops.filt %<>%  mutate(drop = ifelse(Region %in% c(4) & year>=2018,1,0)) %>% filter(drop==0) %>% select(-c(drop, Region)) 
ops.filt %<>% select(-c(vess_grt,vess_loa,vess_callsign))

ops.filt <-rbind(ops.filt,add.data)
 
save(ops.filt,file="hybird.ops.model.Rdata")

##to plot out the locations of the data in region 4
library(sf)
plot_data<-ops.filt_new %>% filter(Region==4&year >2017)

world_map <- map_data("world")
 
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

#Add map to base plot
base_world <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour="gray", fill="gray") +theme_bw()

png("C:/Users/nany/Documents/data_map_2018.png", width = 1200, height = 1000)
base_world+ ylim(-50,0) +xlim(160,280) +
  geom_point(data=plot_data,aes(x=lond,y=latd,color=as.factor(year))) +geom_segment(aes(x=230,y=-50,xend=230,yend=-5),linetype="dashed",color="red")+
  geom_segment(aes(x=210,y=-5,xend=210,yend=0),linetype="dashed",color="red")+ geom_segment(aes(x=210,y=-5,xend=230,yend=-5),linetype="dashed",color="red")
dev.off()

plot_data<-ops.filt_new %>% filter(Region==4&year==2017)

png("C:/Users/nany/Documents/data_map_2017.png", width = 1200, height = 1000)
base_world+ ylim(-50,0) +xlim(160,280) +
  geom_point(data=plot_data,aes(x=lond,y=latd,color=as.factor(year))) +geom_segment(aes(x=230,y=-50,xend=230,yend=-5),linetype="dashed",color="red")+
  geom_segment(aes(x=210,y=-5,xend=210,yend=0),linetype="dashed",color="red")+ geom_segment(aes(x=210,y=-5,xend=230,yend=-5),linetype="dashed",color="red")
dev.off()


##try alb cpue
png("C:/Users/nany/Documents/nominal_CPUE_compare.png", width = 1200, height = 1000)

print(ops.filt %>%
  filter(Region %in% c(1:4)) %>%
  group_by(YrQtr,Region, cate) %>%
  summarise(cpue=mean(alb_cpue)) %>%
  ggplot(aes(x=YrQtr,y=cpue,color=cate)) + geom_line()+facet_grid(Region~.) + geom_vline(xintercept=2018,linetype = "dashed")+theme_bw())

dev.off()

png("P:/SLOTH/alb_MSE/Geo_CPUE/Region_4_data_allflag.png", width = 1600, height = 1000)

ops.filt_new %>% filter(Region==3) %>%
  group_by(year) %>%
  summarise(set.no=n()) %>%
  ggplot(aes(x=as.factor(year),y=set.no)) + geom_bar(stat="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title = element_text(size=16),legend.text = element_text(size=16),axis.title =element_text(size=20), axis.text = element_text(size=16),strip.text.y = element_text(size = 16),legend.position = "bottom")

dev.off()


Region_4_record <- ops.filt_new %>% filter(Region==4) %>%
  group_by(year) %>%
  summarise(set.no=n())

write.csv(Region_4_record,"no_region_4.csv")


