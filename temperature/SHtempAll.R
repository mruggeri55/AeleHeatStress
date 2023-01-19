library(ggplot2)
library(dplyr)
library(gridExtra)
library(lme4)
library(MASS)
library(EnvStats)
library(multcomp)
library(lubridate)

setwd('/Users/maria/Desktop/Kenkel_lab/Anthopleura/temp/Shark Harbor')

L1=read.csv('SH_low_Nov2019July2020.csv')
L2=read.csv('SH_low_Nov2020July2021.csv')
L3=read.csv('SH_low_Aug2021Sep2021.csv')
L4=read.csv('SH_low_Dec2021.csv')
L5=read.csv('SH_low_JuneAugust2022.csv')
L6=read.csv('SHlow_AugOct2022.csv')
L7=read.csv('SH_low_OctNov2022.csv')
low=as.data.frame(rbind(L1[,1:2],L2[,1:2],L3[,1:2],L4[,1:2],L5[,1:2],L6[,1:2],L7[,1:2]))
low$Date_time=mdy_hm(low$Date_Time)
low$month=month(low$Date_time)
low$day=mday(low$Date_time)
low$mday=paste(low$month,low$day)
low$zone=rep('low',length(low$Date_Time))

H1=read.csv('SH_High_Nov2019Jan2020.csv')
H2=read.csv('SH_high_Aug2020Nov2020.csv')
H3=read.csv('SH_high_Aug2021Sep2021.csv')
H4=read.csv('SH_high_Dec2021.csv')
H5=read.csv('SH_high_JuneAugust2022.csv')
H6=read.csv('SHhigh_AugOct2022.csv')
H7=read.csv('SH_high_OctNov2022.csv')
high=as.data.frame(rbind(H1[,1:2],H2[,1:2],H3[,1:2],H4[,1:2],H5[,1:2],H6[,1:2],H7[,1:2]))
high$Date_time=mdy_hm(high$Date_Time)
high$month=month(high$Date_time)
high$day=mday(high$Date_time)
high$mday=paste(high$month,high$day)
high$zone=rep('high',length(high$Date_Time))

all=as.data.frame(rbind(low,high))
all=na.omit(all)
all$year=year(all$Date_time)
all$my=paste(all$month,all$year,sep = '_')

#write.csv(all,'SH_Nov2019to2022.csv')
#all=read.csv('SH_Nov2019to2022.csv',row.names = 1)

ggplot(all,aes(x=Date_time,y=Temp,color=zone))+geom_line()

#subset by year
ggplot(all[all$year==c(2019,2020),],aes(x=Date_time,y=Temp,color=zone))+geom_line()
ggplot(all[all$year==2020,],aes(x=Date_time,y=Temp,color=zone))+geom_line()
ggplot(all[all$year==2021,],aes(x=Date_time,y=Temp,color=zone))+geom_line()
ggplot(all[all$year==2022,],aes(x=Date_time,y=Temp,color=zone))+geom_line()

#subset times with high and low data
ggplot(all[all$my %in% c("11_2019","12_2019","1_2020"),],aes(x=Date_time,y=Temp,color=zone))+geom_line()
ggplot(all[all$my %in% c("8_2021","9_2021"),],aes(x=Date_time,y=Temp,color=zone))+geom_line()
ggplot(all[all$my %in% c("12_2021","1_2022"),],aes(x=Date_time,y=Temp,color=zone))+geom_line()
ggplot(all[all$my %in% c("6_2022","7_2022","8_2022","9_2022","10_2022","11_2022"),],aes(x=Date_time,y=Temp,color=zone))+
  geom_line()+
  scale_color_manual(values=c("#BF812D","#80CDC1"))+
  xlab('')+ylab('Temperature (\u00B0C)')+
  theme_bw(base_size = 20)


# compare high and low
# may not be relevant because do not have high and low at same timepoints
mean(all$Temp)
all %>% group_by(zone) %>% summarise(temp=mean(Temp),sd=sd(Temp),n=n())
# zone   temp    sd      n
# <chr> <dbl> <dbl>  <int>
# 1 high   20.4  6.05  90863
# 2 low    18.0  3.44 191730
# so way more data points for low intertidal

ggplot(all,aes(x=zone,y=Temp))+geom_boxplot()

all %>% group_by(zone,month) %>% summarise(temp=mean(Temp),sd=sd(Temp),n=n())
# zone  month  temp    sd     n
# <chr> <dbl> <dbl> <dbl> <int>
# 1 high      1  15.7  5.45  8699
# 2 high      6  21.4  5.35  3350
# 3 high      7  22.0  4.73  8928
# 4 high      8  23.4  4.77 14879
# 5 high      9  23.1  5.50 20604
# 6 high     10  22.0  5.97 11126
# 7 high     11  17.5  5.32 11836
# 8 high     12  15.4  3.85 11441
# 9 low       1  15.8  3.39 18663
# 10 low       2  16.4  4.05 16416
# 11 low       3  16.3  3.45 17856
# 12 low       4  17.0  3.17 17280
# 13 low       5  18.1  2.25 17856
# 14 low       6  18.8  1.92 20630
# 15 low       7  20.5  1.54 20063
# 16 low       8  21.5  1.51 11033
# 17 low       9  22.2  1.78  9054
# 18 low      10  21.0  1.75  2192
# 19 low      11  17.8  3.15  5700
# 20 low      12  15.8  2.73 26340

# so we have data from all months for low intertidal, but missing Feb-May for high intertidal
# let's subset June-Jan
months=c(6,7,8,9,10,11,12,1)
ggplot(all[all$month %in% months,],aes(x=factor(month, levels=months),y=Temp,fill=zone))+geom_boxplot()

# plot months we have data for both
sub=all[all$month %in% months,]
subsum=as.data.frame(all %>% group_by(month,zone) %>% summarise(mean=mean(Temp),sd=sd(Temp)))
pd=position_dodge(0.3)
ggplot(subsum,aes(x=factor(month, levels= c(6,7,8,9,10,11,12,1,2,3,4,5)),y=mean,color=zone,group=zone))+
  geom_point(position = pd)+
  geom_line(position = pd)+
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),position = pd, width = 0.2)+
  ylab('temperature')+
  xlab('month')+
  theme_classic(base_size=20)

subyear=as.data.frame(all %>% group_by(month,year) %>% summarise(mean=mean(Temp),sd=sd(Temp)))
pd=position_dodge(0.3)
ggplot(subyear,aes(x=factor(month, levels= c(6,7,8,9,10,11,12,1,2,3,4,5)),y=mean,color=as.factor(year),group=as.factor(year)))+
  geom_point(position = pd)+
  geom_line(position = pd)+
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),position = pd, width = 0.2)+
  ylab('temperature')+
  xlab('month')+
  theme_classic(base_size=20)


Nov=all[all$month==11,]
mean(Nov$Temp) #17.74
Nov %>% group_by(zone) %>% summarise(mean=mean(Temp),sd=sd(Temp))
# zone   mean    sd
# <chr> <dbl> <dbl>
# 1 high   17.6  5.46
# 2 low    17.9  3.20

Fall=all[all$month %in% c(9,10,11),]
mean(Fall$Temp) #20.78

##### center vs edge ######
inside=read.csv('SH_inside_(low).csv')
inside$pos=rep('center',times=length(inside$Date_Time))
out=read.csv('SH_edge_(low).csv')
out$pos=rep('edge',times=length(out$Date_Time))

pos_temp=as.data.frame(rbind(inside,out))
pos_temp$Date_time=mdy_hm(pos_temp$Date_Time)
pos_temp$month=month(pos_temp$Date_time)
pos_temp$day=mday(pos_temp$Date_time)
pos_temp$pos=factor(pos_temp$pos,levels=c('edge','center'))

#write.csv(pos_temp,'SHtemp_CenterEdge.csv')

ggplot(pos_temp,aes(x=Date_time,y=Temp,color=pos))+
  geom_line()+
  scale_color_manual(values=c("#BF812D","#80CDC1"))+
  xlab('')+ylab('Temperature (\u00B0C)')+
  theme_bw(base_size = 20)

pos_temp %>% group_by(pos) %>% summarise(temp=mean(Temp),sd=sd(Temp),n=n())
# pos     temp    sd     n
# 1 edge    16.5  2.02 47975
# 2 center  16.5  1.93 47975
t.test(pos_temp$Temp~pos_temp$pos) 
# ns

pos_temp$mday=paste(pos_temp$month,pos_temp$day)
pos_range=pos_temp %>% group_by(pos,mday) %>% summarise(min=min(Temp),max=max(Temp),range=(max(Temp)-min(Temp))) 
pos_range %>% group_by(pos) %>% summarise(mean=mean(range),max=max(max),min=min(min))
# pos     mean   max   min
# <fct>  <dbl> <dbl> <dbl>
# 1 edge    6.22  29.1  9.26
# 2 center  5.76  29.4  9.99
t.test(pos_range$range~pos_range$pos)
#ns but edge on average ~0.5 C greater variation than center

##### sub 2022 #######
sub2022=all[all$my %in% c("6_2022","7_2022","8_2022","9_2022","10_2022","11_2022"),]
sub2022 %>% group_by(zone) %>% summarise(temp=mean(Temp),sd=sd(Temp),n=n())
# zone   temp    sd     n
# <chr> <dbl> <dbl> <int>
# 1 high   22.1  5.34 40697
# 2 low    20.9  2.14 40681
t.test(sub2022$Temp~sub2022$zone) # use Welch's t-test because variances unequal
# t = 42.934, df = 53455, p-value < 2.2e-16

drange=sub2022 %>% group_by(zone,mday) %>% summarise(min=min(Temp),max=max(Temp),range=(max(Temp)-min(Temp))) 

drange %>% group_by(zone) %>% summarise(mean=mean(range),max=max(max),min=min(min))
# zone   mean   max   min
# 1 high  16.1   42.6  8.91
# 2  low   6.34  33.0  9.24
t.test(drange$range~drange$zone)
# t = 18.28, df = 251.58, p-value < 2.2e-16
ggplot(drange,aes(x=zone,y=range,fill=zone,group=zone))+geom_boxplot()+theme_bw(base_size=20)+ylab('daily temp range (\u00B0C)')




