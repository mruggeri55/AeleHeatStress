library(ggplot2)
library(dplyr)
library(gridExtra)
library(lme4)
library(MASS)
library(EnvStats)
library(multcomp)
library(lubridate)

HI=read.csv('/Users/maria/Desktop/Kenkel_lab/Anthopleura/temp/Point Dume/PD_High_In.csv')
HO=read.csv('/Users/maria/Desktop/Kenkel_lab/Anthopleura/temp/Point Dume/PD_high_out.csv')

HI$date_time=mdy_hm(HI$Date_Time)
HO$date_time=mdy_hm(HO$Date_Time)

HI=HI[25:47974,]
HO=HO[1:47950,]
HI$position=rep('in',times=47950)
HO$position=rep('out',times=47950)

# inside
mean(HI$Temp) #15.11866 
sd(HI$Temp) #1.930977
max(HI$Temp) #25.95
min(HI$Temp) #6.95
# outside
mean(HO$Temp) #15.15947
sd(HO$Temp) #2.208346
max(HO$Temp) #28.95
min(HO$Temp) #6.26

#concatenate inside and outside loggers
conc=rbind(HI,HO)
conc$position=as.factor(conc$position)
conc$position=relevel(conc$position, 'out')

conc %>% group_by(position) %>% summarise(temp=mean(Temp),sd=sd(Temp),n=n())
# position  temp    sd     n
# <fct>    <dbl> <dbl> <int>
# 1 out       15.2  2.21 47950
# 2 in        15.1  1.93 47950
t.test(conc$Temp~conc$pos) 
# sig p=0.002316

conc$month=month(conc$date_time)
conc$day=mday(conc$date_time)
conc$mday=paste(conc$month,conc$day)
pos_range=conc %>% group_by(position,mday) %>% summarise(min=min(Temp),max=max(Temp),range=(max(Temp)-min(Temp))) 
pos_range %>% group_by(position) %>% summarise(mean=mean(range),max=max(max),min=min(min))
# position  mean   max   min
# <fct>    <dbl> <dbl> <dbl>
# 1 out       6.47  29.0  6.26
# 2 in        5.53  26.0  6.95
t.test(pos_range$range~pos_range$position)
# sig p=0.001135

#quartz()
ggplot(conc, aes(x=date_time,y=Temp,group=position,color=position))+geom_line()+theme_bw(base_size = 20)+xlab('')+ylab('Temperature °C')
#quartz.save('/Users/maria/Desktop/Kenkel_lab/PD_inVout_temp.pdf',type='pdf')
#dev.off()
ggplot(conc,aes(x=Temp,fill=position))+geom_density(alpha=0.5)+theme_classic()
anova=aov(conc$Temp~conc$position)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# conc$position     1     40   39.94   9.281 0.00232 **
# Residuals     95898 412623    4.30                   

conc$month_day=paste(month(conc$date_time),day(conc$date_time))
daily_range = group_by(conc,position,month_day) %>% summarise(range=max(Temp)-min(Temp))
daily_range_sum=group_by(daily_range,position) %>% summarise(mean=mean(range),sd=sd(range))
# position  mean    sd
# <fct>    <dbl> <dbl>
# 1 out       6.47  2.92
# 2 in        5.53  2.32
anova=aov(range~position,data=daily_range)
summary(anova)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# position      1   75.1   75.05   10.79 0.00113 **
# Residuals   334 2323.3    6.96

drange=as.data.frame(daily_range)
ggplot(drange,aes(x=position,y=range))+geom_boxplot()
quartz()
ggplot(drange,aes(x=range,fill=position))+geom_density(alpha=0.5)+theme_classic()+xlab('Daily Range')+
  geom_vline(data=daily_range_sum,aes(xintercept=mean,color=position),linetype="dashed",size=1)
quartz.save('/Users/maria/Desktop/Kenkel_lab/PD_inVout_TempDrangeDensity.pdf',type='pdf')
dev.off()

daily_max = group_by(conc,position,month_day) %>% summarise(max=max(Temp))
daily_max_sum=group_by(daily_max,position) %>% summarise(mean=mean(max),sd=sd(max))
#   position  mean    sd
# 1 out       18.3  3.17
# 2 in        17.7  2.51

max=as.data.frame(daily_max)
ggplot(max,aes(x=max,fill=position))+geom_density(alpha=0.5)+theme_classic()

daily_min = group_by(conc,position,month_day) %>% summarise(min=min(Temp))
daily_min_sum=group_by(daily_min,position) %>% summarise(mean=mean(min),sd=sd(min))
# position  mean    sd
# 1 out       11.9  2.03
# 2 in        12.2  1.89
min=as.data.frame(daily_min)
ggplot(min,aes(x=min,fill=position))+geom_density(alpha=0.5)+theme_classic()

