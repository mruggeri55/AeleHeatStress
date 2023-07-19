setwd('/Users/maria/Desktop/Kenkel_lab/Anthopleura/Nov2020_AeleHeatStress/')

SymToHost=read.csv('qPCR/qPCR_master.csv',row.names = 1)
SymToHostSub=subset(SymToHost,select=c('Sample','fold_change'))
chl=read.csv("physiology/analysis/AHS_chl_quant.csv")
chl=subset(chl,select=c('sampleID','Chla_norm','tot_chl_norm'))
PAM=read.csv('PAM/PAMexpMaster.csv',row.names=1)
PAM$sample_rep=paste(PAM$sample,PAM$tank,sep = '-')
PAM_EndPt=PAM[PAM$Date == '11/30/20',]
PAM_EndPtSub=subset(PAM_EndPt,select=c('sample_rep','MQY'))
size=read.csv('physiology/size.csv')
colnames(size)=c('sampleID','weight')

meta=subset(PAM_EndPt,select=c('sample_rep','genotype','pos','Zone','treatment','tank'))
colnames(meta)[2]='agg'
meta$agg=factor(meta$agg,levels=c('1','2','3','4','5','11','7','8','9','10','12','13'))
agg=c('1','2','3','4','5','11','7','8','9','10','12','13')
genet=c('G1','G1','G1','G4','G5','G10','G1','G1','G9','G10','G12','G13')
agg_genet=as.data.frame(cbind(agg,genet))
meta2=merge(meta,agg_genet,by='agg')
meta2$genet=factor(meta2$genet,levels=c('G4','G5','G1','G10','G9','G12','G13'))
meta=meta2
colnames(meta)[2]='sampleID'

TideHeight=read.csv('TidalHeight_7Nov2022.csv')
colnames(TideHeight)[1]='agg'
#TideHeightMid=TideHeight[,c(2,4)]
meta2=merge(meta,TideHeight, by='agg')
meta=meta2
#write.csv(meta,'meta.csv')

colnames(SymToHostSub)=c('sampleID','SymHostRatio')
colnames(PAM_EndPtSub)[1]='sampleID'

library(tidyverse)
#put all data frames into list
df_list <- list(meta,PAM_EndPtSub,chl,SymToHostSub,size)
#merge all data frames in list
trait_master = df_list %>% reduce(full_join, by='sampleID')
trait_master$pos=as.factor(trait_master$pos)
trait_master$Zone=as.factor(trait_master$Zone)
trait_master$agg=as.factor(trait_master$agg)
trait_master$treatment=as.factor(trait_master$treatment)
trait_master$tank=as.factor(trait_master$tank)
trait_master$genet=as.factor(trait_master$genet)

trait_master=trait_master[!is.na(trait_master$pos),] # remove extra samples not in expt

# this worked but for some reason 4 entries for 12-I-3
# 12-I-3 2 MQY values -- look at tank layout from last day (I think I just input sample name wrong?)
# let's remove dups for 12-I-3 and replace MQY with NA
trait_master=trait_master[-c(43,44,49),]
trait_master[trait_master$sampleID=='12-I-3',9]<-NA

#write.csv(trait_master,'trait_master.csv')


#### normality ####
trait=trait_master$MQY # set trait of interest here

# is trait normally distributed?
hist(trait)
boxplot(trait)
outliers=boxplot(trait)$out
trait_master[trait %in% outliers,] # outliers for sym to host all have CV <10%, no consistent outliers across traits
clean=trait_master[!trait %in% outliers,]
trait_clean=clean$MQY #change to trait of interest

library(e1071)
skewness(na.omit(trait))
skewness(na.omit(trait_clean))
# for right skew -- tail is on the right 
skewness(log(na.omit(trait)))
skewness(sqrt(na.omit(trait)))
# for left skew -- tail is on the left
skewness(na.omit(trait)^2)

library(car)
qqPlot(na.omit(trait)) # check if residuals normal
qqPlot(na.omit(trait_clean))
qqPlot(na.omit(log(trait)))
qqPlot(na.omit(trait)^2)
       
shapiro.test(na.omit(trait)) # stat test for normality
shapiro.test(na.omit(trait_clean))
shapiro.test(na.omit(log(trait)))
shapiro.test(na.omit(trait)^2)

# so MQY more normal if we remove outliers but shap test still sig, looks best when squared but still off
# Chl looks pretty normal, two outliers are skewing results but I do not see good reason to remove them
# SymHost most normal when log transformed, but clean also looks better, let's use log transform to keep all data
# weight needs log transform
trait_master$logSymHost=log(trait_master$SymHostRatio)
trait_master$logWeight=log(trait_master$weight)


########## start here ##########
trait_master=read.csv('trait_master.csv',row.names = 1)
trait_master$agg=factor(trait_master$agg,levels=c('1','2','3','4','5','11','7','8','9','10','12','13'))
trait_master$genet=factor(trait_master$genet,levels=c('G4','G5','G1','G10','G9','G12','G13'))
trait_master$tank=as.factor(trait_master$tank)
trait_master$pos=as.factor(trait_master$pos)
trait_master$Zone=as.factor(trait_master$Zone)
trait_master$treatment=as.factor(trait_master$treatment)
trait_master$logSymHost=log(trait_master$SymHostRatio)
trait_master$logWeight=log(trait_master$weight)

#trait=trait_master$logSymHost
#trait=trait_master$logWeight # note for weight, treatment and tank not relevant bc only initial
#trait=trait_master$Chla_norm
trait=trait_master$MQY

#### MIXED MODELS ####
library(lme4)
library(lmerTest)
library(emmeans)
library(sjPlot)

MixModel=lmer(trait ~ treatment*Zone + pos + (1|genet) + (1|tank:treatment), data=trait_master)

# checking residuals
qqPlot(residuals(MixModel))
plot(MixModel,which=1)

library(piecewiseSEM) 
# marginal R2 considers only the variance by the fixed effects, and the conditional R2 by both the fixed and random effects
rsquared(MixModel)
# MQY 
# r2=0.554 treatment*Zone*pos + (1|genet) + (1|tank:treatment)
# adding in tidal height does not improve model r2=0.546 

# Sym to host (log)
#r2 does not really change when removing terms, conditional ~55%
# adding in tidal height improves model r2=0.614 compared to r2=0.565 with Zone only

# Chla_norm 
# r2=0.389 treatment*Zone*pos + (1|genet) + (1|tank:treatment)
# using tidal height instead of Zone does not change r2

rand(MixModel)
anova(MixModel,ddf='Kenward-Roger')
pairs(emmeans(MixModel,c('Zone','treatment')))
pairs(emmeans(MixModel,'treatment'))
emmeans(MixModel,c('Zone','treatment'))
emmeans(MixModel,'treatment')

# export tables
MQY=lmer(MQY ~ treatment*Zone + pos + (1|genet) + (1|tank:treatment), data=trait_master)
SymHost=lmer(logSymHost ~ treatment*Zone + pos + (1|genet) + (1|tank:treatment), data=trait_master)
Chla_norm=lmer(Chla_norm ~ treatment*Zone + pos + (1|genet) + (1|tank:treatment), data=trait_master)
tab_model(MQY,SymHost,Chla_norm,show.df = TRUE)

# using continuous tidal height measurements
trait_master$tidal_height=trait_master$actual.mid
MQY=lmer(MQY ~ treatment*tidal_height + pos + (1|genet) + (1|tank:treatment), data=trait_master)
SymHost=lmer(logSymHost ~ treatment*tidal_height + pos + (1|genet) + (1|tank:treatment), data=trait_master)
Chla_norm=lmer(Chla_norm ~ treatment*tidal_height + pos + (1|genet) + (1|tank:treatment), data=trait_master)
tab_model(MQY,SymHost,Chla_norm,show.df = TRUE)

# binning aggregation 3 in high intertidal
trait_master$Zone_v2=trait_master$Zone
trait_master$Zone_v2[trait_master$agg=='3']='high'
MQY=lmer(MQY ~ treatment*Zone_v2 + pos + (1|genet) + (1|tank:treatment), data=trait_master)
SymHost=lmer(logSymHost ~ treatment*Zone_v2 + pos + (1|genet) + (1|tank:treatment), data=trait_master)
Chla_norm=lmer(Chla_norm ~ treatment*Zone_v2 + pos + (1|genet) + (1|tank:treatment), data=trait_master)
tab_model(MQY,SymHost,Chla_norm,show.df = TRUE)


### Chl ####
# treatment and genet sig (trait ~ treatment*Zone + pos + (1|genet) + (1|tank:treatment))
# treatment and genet sig (trait ~ treatment*Zone*pos + genet + (1|tank:treatment))
# trmt sig, genet marginally sig 0.09 (trait ~ treatment*Zone + genet + (1|tank:treatment))
# no fixed effects sig (trait ~ treatment*genet + (1|tank:treatment))
# if using tidal height, only random effects significant i.e. genet and tank:treatment


#### if residuals not good 
ModelClean=lmer(trait_clean ~ treatment*Zone + pos + (1|genet) + (1|tank:treatment), data=clean)
#ModelRed=lmer(trait_clean ~ treatment*Zone + (1|genet) + (1|tank:treatment), data=clean)

clean$Zone_v2=clean$Zone
clean$Zone_v2[clean$agg=='3']='high'
ModelClean=lmer(trait_clean ~ treatment*Zone_v2 + pos + (1|genet) + (1|tank:treatment), data=clean)

# checking residuals
qqPlot(residuals(ModelClean))
plot(ModelClean,which=1)

rand(ModelClean)
anova(ModelClean,ddf='Kenward-Roger')

####### does acclimation effect response to treatment within a genet #####

Genet=trait_master[trait_master$genet=='G10',] #change to genet of interest
trait=Genet$MQY #change to trait of interest
#Genet=Genet[!is.na(trait),] #remove any NAs
#trait=log(Genet$SymHostRatio) #change to trait of interest
#Genet=droplevels(Genet)

GenetMod=lmer(trait ~ treatment*agg + (1|tank:treatment), data=Genet)
#GenetMod=lmer(trait ~ treatment*Zone + (1|tank:treatment), data=Genet)
#GenetMod=lmer(trait ~ treatment*actual.mid + (1|tank:treatment), data=Genet)
#GenetMod=lmer(trait ~ treatment*agg*pos + (1|tank:treatment), data=Genet)


qqPlot(residuals(GenetMod))
plot(GenetMod,which=1)

rand(GenetMod)
anova(GenetMod,ddf='Kenward-Roger')

pairs(emmeans(GenetMod,c('agg','treatment')))

emmeans(GenetMod,c('agg','treatment'))

MQYmod=lmer(MQY ~ treatment*agg + (1|tank:treatment), data=Genet)
MQY=anova(MQYmod,ddf='Kenward-Roger')
rand1=as.data.frame(rand(MQYmod))
logSymHostMod=lmer(logSymHost ~ treatment*agg + (1|tank:treatment), data=Genet)
logSymHost=anova(logSymHostMod,ddf='Kenward-Roger')
rand2=as.data.frame(rand(logSymHostMod))
ChlaMod=lmer(Chla_norm ~ treatment*agg + (1|tank:treatment), data=Genet)
Chla_norm=anova(ChlaMod,ddf='Kenward-Roger')
rand3=as.data.frame(rand(ChlaMod))

tab_dfs(list(MQY,rand1,logSymHost,rand2,Chla_norm,rand3),show.rownames = T,
        titles=c('MQY','','logSymHost','','Chla_norm',''))

# exporting pairwise data
pMQY=pairs(emmeans(MQY,c('agg','treatment')))
pSymHost=pairs(emmeans(logSymHost,c('agg','treatment')))
pChla_norm=pairs(emmeans(Chla_norm,c('agg','treatment')))

pairs_MQY=cbind(as.data.frame(pMQY),'trait'=rep('MQY',times=length(pMQY)))
pairs_SymHost=cbind(as.data.frame(pSymHost),'trait'=rep('SymHost',times=length(pSymHost)))
pairs_Chla_norm=cbind(as.data.frame(pChla_norm),'trait'=rep('Chla_norm',times=length(pChla_norm)))
pMaster=as.data.frame(rbind(pairs_MQY,pairs_SymHost,pairs_Chla_norm))
write.csv(pMaster,'tables/Genet1_AggTrmt_pw.csv')

#### CHLA ####
# genet 1 no agg*treatment int, only trmt sig
# again genet 10 no terms sig

### MQY ####
# residuals still skewed but better without outliers  
# treatment, zone, and the interaction consistently significant,no treatment*genet interaction, no effect of pos, no treatment*pos interaction, genet as random effect ns
# Genet 1: sig agg, tmt, agg*trmt interaction, no tidal height*treatment int or Zone*treatment though
# Genet 10 singular fit

### sym to host ###
# trmt, trmt*Zone sig if log transformed, and genet sig as fixed and random, no genet*trmt int
# Genet 1: sig agg, sig trmt, no agg*trmt or Zone*trmt for raw or log, no TideHeight*trmt int
# Genet 10: sig trmt, marg sig agg*trmt (0.09) for log transform, marg sig Zone*trmt and TideHeight*trmt

# checking which genets differ
library(emmeans)
trait=trait_master$MQY
MixModel=lmer(trait ~ treatment*Zone + pos + genet + (1|tank:treatment), data=trait_master) # no treatment*genet effect for any traits
anova(MixModel,ddf='Kenward-Roger')
genets=emmeans(MixModel, 'genet')
pw_df=as.data.frame(pairs(genets))
# change file output
write.csv(pw_df,'pw_MQY_by_genet.csv')

### sym to host
# G1 vs G12 p=0.0536
# G1 vs G13 p=0.006

### Chla_norm
# G1 vs G10 p=0.0203
# G1 vs G4 p=0.0713

ggplot(trait_master,aes(x=genet,y=MQY,fill=treatment))+geom_boxplot()+theme_bw(base_size = 20)
ggplot(trait_master,aes(x=genet,y=logSymHost,fill=treatment))+geom_boxplot()+theme_bw(base_size = 20)
ggplot(trait_master,aes(x=genet,y=Chla_norm,fill=treatment))+geom_boxplot()+theme_bw(base_size = 20)+ylab('ug Chla/mg protein')



#### PLOTS ####
###### by trait
# color option 1 c('#EF6548','#4EB3D3')

MQY=ggplot(trait_master[!is.na(trait_master$MQY),], aes(x=treatment,y=MQY,fill=Zone))+
  geom_boxplot()+
  scale_fill_manual(values=c("#BF812D","#80CDC1"))+
  theme_classic(base_size = 15)+
  ylim(0.45,0.8) #extend the limits a bit so I can fit in pvals
SymToHost=ggplot(trait_master[!is.na(trait_master$SymHostRatio),], aes(x=treatment,y=log(SymHostRatio),fill=Zone))+
  geom_boxplot()+
  scale_fill_manual(values=c("#BF812D","#80CDC1"))+
  theme_classic(base_size = 15)+
  ylab('log(sym/host ratio)')+
  ylim(-5,0)
Chla=ggplot(trait_master[!is.na(trait_master$Chla_norm),], aes(x=treatment,y=Chla_norm,fill=Zone))+
  geom_boxplot()+
  scale_fill_manual(values=c("#BF812D","#80CDC1"))+
  theme_classic(base_size = 15)+
  ylab('ug Chla/mg protein')+
  ylim(0,7)

library(cowplot)
legend=get_legend(Chla)

library(gridExtra)
library(ggpubr)
grid.arrange(MQY+theme(legend.position = 'none')+xlab(''),
             SymToHost+theme(legend.position = 'none'),
             Chla+theme(legend.position = 'none')+xlab(''),
             legend,nrow=1,widths=c(3,3,3,1))

# agg rxn norms
sum_stats=trait_master %>% group_by(treatment,agg,actual.mid) %>% summarise(mean=mean(na.omit(MQY)))

ggplot(as.data.frame(sum_stats), aes(x=treatment,y=mean,color=actual.mid))+geom_boxplot()+ 
  stat_summary(fun=mean, geom="line", aes(group = agg, color=actual.mid))+
  scale_color_gradientn(colors=rev(colorRampPalette(brewer.pal(11, "BrBG"))(12)))+
  theme_bw(base_size=20)+ylab('MQY')

ggplot(trait_master,aes(x=MQY,y=logSymHost))+geom_point()+geom_smooth(method='lm')
summary(lm(trait_master$logSymHost~trait_master$MQY))


# relationship of sym to host and MQY by tidal height
trait_master$Zone_trmt=paste(trait_master$Zone,trait_master$treatment)

ggplot(trait_master,aes(x=MQY,y=logSymHost,color=Zone_trmt,fill=Zone_trmt))+
  geom_point(size=3,shape=21,stroke=1)+
  scale_fill_manual(name='zone treatment',
                    labels=c('high control','high heat','low control','low heat'),
                    values=c(alpha("#BF812D"),alpha("#BF812D",0.3),alpha("#35978F"),alpha("#35978F",0.3)),
                    breaks=c('high C','high H','low C','low H'),
                    guide='legend') +
  scale_colour_manual(name='zone treatment',
                      labels=c('high control','high heat','low control','low heat'),
                      values=c("#BF812D","#BF812D","#35978F","#35978F"),
                      breaks=c('high C','high H','low C','low H')) +
  theme_bw(base_size=20)+ylab('log(sym/host)')

# change to MQY or logSymHost
ggplot(trait_master,aes(x=MQY,fill=Zone_trmt,color=Zone_trmt))+
  geom_density()+
  scale_fill_manual(values=c(alpha("#BF812D",0.8),alpha("#BF812D",0.3),alpha("#35978F",0.8),alpha("#35978F",0.3)))+
  scale_color_manual(values=c("#BF812D","#BF812D","#35978F","#35978F"))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

mod=summary(lm(logSymHost~MQY+Zone,trait_master))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -9.0788     0.9387  -9.672  < 2e-16 ***
# MQY           9.8437     1.4104   6.979  1.3e-10 ***
# Zonelow       0.3650     0.1409   2.590   0.0107 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7561 on 132 degrees of freedom
# (8 observations deleted due to missingness)
# Multiple R-squared:  0.2696,	Adjusted R-squared:  0.2585 
# F-statistic: 24.36 on 2 and 132 DF,  p-value: 9.908e-10
tab_model(mod)


###### by genet
#trait_master$ChlSymRatio=trait_master$Chla_norm/trait_master$SymHostRatio
Genet=trait_master[trait_master$genet=='G1',] # or G10
trait='MQY'
genet='genet 1'

library(ggplot2)
ggplot(Genet,aes(x=agg,y=MQY,fill=treatment))+geom_boxplot()+
  ylab(trait)+xlab('aggregation')+
  ggtitle(genet)+
  theme_bw(base_size = 15)

library(dplyr)
library(RColorBrewer)

Genet$agg=droplevels(Genet$agg)
Genet$actual.mid=as.factor(Genet$actual.mid)
brewer.pal(9,'Set1')
display.brewer.pal(9,'BrBG')
brewer.pal(9,'BrBG')
#color=c("#984EA3","#377EB8","#4DAF4A","#FF7F00","#E41A1C") # genet 1 rainbow
color=c("#01665E","#35978F","#DFC27D","#BF812D","#8C510A") # genet 1 tidal
#color=c("#984EA3","#E41A1C") # genet 10 rainbow
#color=c("#35978F","#BF812D") # genet 10 tidal

#### MQY
MQY_sum = Genet[!is.na(Genet$MQY),] %>% group_by(actual.mid,treatment) %>% summarise(mean=mean(MQY),sd=sd(MQY),n=n())
MQY_sum$se = MQY_sum$sd/sqrt(MQY_sum$n)
pd=position_dodge(0.4)
MQY=ggplot(MQY_sum,aes(x=treatment,y=mean,color=actual.mid))+geom_point(position=pd,size=1.5)+
  geom_line(aes(group=actual.mid),position=pd,size=1)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=actual.mid),position=pd,width=0.3,size=1)+
  scale_color_manual(values=color)+
  ylab('MQY')+
  theme_bw(base_size = 15)

#### SymToHost
SymToHost_sum = Genet[!is.na(Genet$SymHostRatio),] %>% group_by(actual.mid,treatment) %>% summarise(mean=mean(log(SymHostRatio)),sd=sd(log(SymHostRatio)),n=n())
SymToHost_sum$se = SymToHost_sum$sd/sqrt(SymToHost_sum$n)
pd=position_dodge(0.4)
SymToHost=ggplot(SymToHost_sum,aes(x=treatment,y=mean,color=actual.mid))+geom_point(position=pd,size=1.5)+
  geom_line(aes(group=actual.mid),position=pd,size=1)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=actual.mid),position=pd,width=0.3,size=1)+
  scale_color_manual(values=color)+
  ylab('log(Sym/Host)')+
  theme_bw(base_size = 15)

#### CHL
Chl_sum = Genet[!is.na(Genet$Chla_norm),] %>% group_by(actual.mid,treatment) %>% summarise(mean=mean(Chla_norm),sd=sd(Chla_norm),n=n())
Chl_sum$se = Chl_sum$sd/sqrt(Chl_sum$n)
pd=position_dodge(0.4)
Chla=ggplot(Chl_sum,aes(x=treatment,y=mean,color=actual.mid))+geom_point(position=pd,size=1.5)+
  geom_line(aes(group=actual.mid),position=pd,size=1)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=actual.mid),position=pd,width=0.3,size=1)+
  scale_color_manual(values=color)+
  ylab('ug Chla/mg protein')+
  theme_bw(base_size = 15)

library(cowplot)
legend=get_legend(Chla)

library(gridExtra)
library(ggpubr)
grid.arrange(MQY+theme(legend.position = 'none')+xlab(''),
             SymToHost+theme(legend.position = 'none')+xlab(''),
             Chla+theme(legend.position = 'none')+xlab(''),
             legend,nrow=1,widths=c(3,3,3,1),top=text_grob(genet,size=18))

grid.arrange(MQY+theme(legend.position = 'none')+xlab('')+ylab(''),
             SymToHost+theme(legend.position = 'none')+xlab('')+ylab(''),
             Chla+theme(legend.position = 'none')+xlab('')+ylab(''),
             legend,nrow=3,ncol=2,widths=c(3,1),top=text_grob(genet,size=18), 
             layout_matrix=cbind(c(1,2,3),c(4,4,4)))



#### PCA all traits #####
pc_df=na.omit(trait_master)
pc=prcomp(pc_df[,8:12],scale.=T)
biplot(pc)

library(ggbiplot)
ggbiplot(pc, groups=pc_df$Zone)
ggbiplot(pc, groups=pc_df$treatment)

#### by aggregation rank ####
trait_master$agg=factor(trait_master$agg, levels=c(2,4,1,5,11,3,10,13,8,12,7,9))
ggplot(trait_master,aes(x=agg,y=logSymHost))+geom_boxplot()
ggplot(trait_master,aes(x=agg,y=Chla_norm))+geom_boxplot()
ggplot(trait_master,aes(x=agg,y=MQY))+geom_boxplot()
ggplot(trait_master,aes(x=agg,y=logWeight))+geom_boxplot()+theme_bw(base_size = 20)+ylab('log(mg wet weight)')+xlab('aggregation')

# plot of tide height aggregation limits
melt=gather(TideHeight,'pos','height',3:5)
ggplot(melt,aes(x=agg,y=height,color=zone))+geom_point()+geom_line()+theme_bw(base_size=20)+xlab('aggregation')+ylab('tidal height (m)')


##### position plots ####
MQY=ggplot(trait_master[!is.na(trait_master$MQY),], aes(x=treatment,y=MQY,fill=pos))+
  geom_boxplot()+
  scale_fill_manual(values=c("#80CDC1","#BF812D"))+
  theme_classic(base_size = 15)+
  ylim(0.45,0.8) #extend the limits a bit so I can fit in pvals
SymToHost=ggplot(trait_master[!is.na(trait_master$SymHostRatio),], aes(x=treatment,y=log(SymHostRatio),fill=pos))+
  geom_boxplot()+
  scale_fill_manual(values=c("#80CDC1","#BF812D"))+
  theme_classic(base_size = 15)+
  ylab('log(sym/host ratio)')+
  ylim(-5,0)
Chla=ggplot(trait_master[!is.na(trait_master$Chla_norm),], aes(x=treatment,y=Chla_norm,fill=pos))+
  geom_boxplot()+
  scale_fill_manual(values=c("#80CDC1","#BF812D"))+
  theme_classic(base_size = 15)+
  ylab('ug Chla/mg protein')+
  ylim(0,7)

library(cowplot)
legend=get_legend(Chla)

library(gridExtra)
library(ggpubr)
grid.arrange(MQY+theme(legend.position = 'none')+xlab(''),
             SymToHost+theme(legend.position = 'none'),
             Chla+theme(legend.position = 'none')+xlab(''),
             legend,nrow=1,widths=c(3,3,3,1))


###### models for weight ######
trait=trait_master$logWeight
ggplot(trait_master,aes(x=genet,y=logWeight))+geom_boxplot()+theme_bw(base_size = 20)+ylab('log(mg wet weight)')
ggplot(trait_master,aes(x=pos,y=logWeight))+geom_boxplot()+theme_bw(base_size = 20)+
  ylab('log(mg wet weight)')+xlab('position')+
  scale_x_discrete(labels=c('center','edge'))

MixModel=lmer(trait ~ Zone + pos + (1|genet), data=trait_master)
#MixModel=lm(trait ~ Zone + pos + genet, data=trait_master)

# checking residuals
qqPlot(residuals(MixModel))
plot(MixModel,which=1)
# marginal R2 considers only the variance by the fixed effects, and the conditional R2 by both the fixed and random effects
rsquared(MixModel)

rand(MixModel)
anova(MixModel,ddf='Kenward-Roger')
pairs(emmeans(MixModel,c('Zone','pos')))
#pairs(emmeans(MixModel,'genet'))
emmeans(MixModel,c('Zone','treatment'))

##### broad sense heritability #####
library(MCMCglmm)
prior1=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002)))
prior2=prior1=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002),G2 = list(V = 1,nu=0.002)))

model1=MCMCglmm(MQY ~ treatment*Zone,random=~genet + tank:treatment,data=trait_master,prior=prior2,nitt=50000,thin=20,burnin=10000)
summary(model1)
plot(model1)
vcv=colnames(model1$VCV)
autocorr.diag(model1$Sol) #for fixed effects
autocorr.diag(model1$VCV) #for random effects
mean(model1$VCV[,1]) # Vg
mean(model1$VCV[,2]+model1$VCV[,3]) # Ve

H1=(model1$VCV[,1])/(model1$VCV[,1]+model1$VCV[,2]+model1$VCV[,3])
mean(H1)
HPDinterval(H1) 
MQYH2=c(mean(H1),HPDinterval(H1))

model2=MCMCglmm(logSymHost ~ treatment*Zone,random=~genet + tank:treatment,data=trait_master,prior=prior2,nitt=50000,thin=20,burnin=10000)
H2=(model2$VCV[,1])/(model2$VCV[,1]+model2$VCV[,2]+model2$VCV[,3])
mean(H2)
HPDinterval(H2) 
SymToHostH2=c(mean(H2),HPDinterval(H2))

model3=MCMCglmm(Chla_norm ~ treatment*Zone,random=~genet + tank:treatment,data=trait_master,prior=prior2,nitt=50000,thin=20,burnin=10000)
H3=(model3$VCV[,1])/(model3$VCV[,1]+model3$VCV[,2]+model3$VCV[,3])
mean(H3)
HPDinterval(H3) 
ChlaH2=c(mean(H3),HPDinterval(H3))

#### now weight using reduced model because only initial weight measured
prior1=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002)))
model4=MCMCglmm(logWeight ~ Zone,random=~genet,data=trait_master,prior=prior1,nitt=50000,thin=20,burnin=10000)
H4=(model4$VCV[,1])/(model4$VCV[,1]+model4$VCV[,2])
mean(H4)
weightH2=c(mean(H4),HPDinterval(H4))

## make df with all traits
H2all=rbind(MQYH2,SymToHostH2,ChlaH2,weightH2)
colnames(H2all)=c("mean","lwr","upper")
H2all=data.frame(H2all)
H2all$trait=c("MQY","Sym/Host","Chla","weight")
H2all$trait=factor(H2all$trait,levels=c("MQY","Sym/Host","Chla","weight"))
write.csv(H2all,'H2alltraits.csv')

ggplot(H2all, aes(x=trait, y=mean)) + 
  geom_bar(stat="identity",colour="black",size=.3) +
  geom_errorbar(aes(ymin=lwr, ymax=upper),size=.3,width=.2) +
  ylim(0,1)+
  theme_bw(base_size = 20)+
  labs(y=expression(H^2))+
  theme(axis.title.x = element_blank()) 

## let's calculate 'bleaching' control-heat for traits and explore their heritability
#### cannot actually do this because now only one entry per genet
## this will only work for MQY using the initial vs final for each individual


#### let's do broad sense heritability of change in MQY ####
PAM_master=read.csv('/Users/maria/Desktop/Kenkel_lab/Anthopleura/Nov2020_AeleHeatStress/PAM/PAMexpMaster.csv',row.names = 1)
initial=PAM_master[PAM_master$Date=='11/20/20',]
final=PAM_master[PAM_master$Date=='11/30/20',]
initial$sample_rep==final$sample_rep
order=initial$sample_rep
final_order=final[match(order,final$sample_rep),]
initial$sample_rep==final_order$sample_rep

delMQY=initial$MQY-final_order$MQY
delMQYdf=as.data.frame(cbind(initial[,13:17],'tank'=initial$tank, delMQY))

delMQYdf$sampleID=gsub(' ','-',delMQYdf$sample_rep)
# there are two 12-I-3's, remove both
delMQYdf=delMQYdf[!delMQYdf$sampleID=='12-I-3',]
order=delMQYdf$sampleID
trait_order=trait_master[match(order,trait_master$sampleID),]
delMQYdf$genet=trait_order$genet
delMQYdf$tank=as.factor(delMQYdf$tank)
write.csv(delMQYdf,'PAM/delMQY.csv')

prior1=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002)))
prior2=prior1=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002),G2 = list(V = 1,nu=0.002)))

m1=MCMCglmm(delMQY ~ treatment*Zone,random=~genet+tank:treatment,data=delMQYdf,prior=prior2,nitt=50000,thin=20,burnin=10000)
summary(m1)
plot(m1)
vcv=colnames(m1$VCV)
autocorr.diag(m1$Sol) #for fixed effects
autocorr.diag(m1$VCV) #for random effects
mean(m1$VCV[,1]) # Vg
mean(m1$VCV[,2]+m1$VCV[,3]) # Ve
H=(m1$VCV[,1])/(m1$VCV[,1]+m1$VCV[,2]+m1$VCV[,3]) 
mean(H)
HPDinterval(H)

delMQYH2=as.data.frame(cbind(mean(H),HPDinterval(H)))
colnames(delMQYH2)=c('mean','lwr','upper')
delMQYH2$trait='MQY_initial - MQY_final'

ggplot(delMQYH2, aes(x=trait,y=mean)) + 
  geom_bar(stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=lwr, ymax=upper),
                size=.3,    # Thinner lines
                width=.2) +
  theme_bw(base_size = 20)+
  labs(y=expression(H^2),x=expression(MQY_initial))+
  theme(axis.title.x = element_blank())+
  ylim(0,1)
