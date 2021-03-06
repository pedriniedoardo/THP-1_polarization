#load the data
setwd("C:/Users/edoardo pedrini/Desktop/New folder (2)")

#read the dfs
df1<-read.csv("161119 CCL22.csv",skip = 19,stringsAsFactors = F)
#remove form this df the data not relative to CCL22
df1<-df1[df1$Target=="CCL22",]
#remove the negative control of CCL22
df1<-df1[df1$Sample!="",]
df2<-read.csv("161119 GAPDH CXCL10.csv",skip = 19,stringsAsFactors = F)
df3<-read.csv("161119 MRC1 TNF.csv",skip = 19,stringsAsFactors = F)
#keep holy the CCL19
df4<-read.csv("161226 CCL19 CTSC.csv",skip = 19,stringsAsFactors = F) %>% .[.$Target=="CCL19",]
df5<-read.csv("161226 CCL5 MSR1.csv",skip = 19,stringsAsFactors = F) %>% .[.$Target=="MSR1",]
#keep only the HSD11B1
#df5<-read.csv("161226 HSD HST.csv",skip = 19,stringsAsFactors = F) %>% .[.$Target=="HSD11B1" & (.$Sample != "NEG"&.$Sample != "EMPTY"),]
#merge the dfs in one df
df<-rbind.data.frame(df1,df2,df3,df4,df5)
#df<-rbind.data.frame(df1,df2,df3,df4)
#summarize the table
table(df$Target,df$Sample)


#look at the variability of the technical replicated
library(ggplot2)
ggplot(df,aes(x=Sample,y=Cq))+geom_boxplot()+facet_wrap(~Target)
ggplot(df,aes(x=Sample,y=Cq))+geom_jitter()+facet_wrap(~Target)
ggplot(df,aes(x=Sample,y=Cq))+geom_boxplot()+facet_wrap(~Target,scales = "free")

#collapse the tech variability calculationg the mean
df_mean<-data.frame(as.list(aggregate(df$Cq~df$Sample+df$Target,FUN = mean)),stringsAsFactors = F)
#give a more appropirte col label
colnames(df_mean)<-c("sample","target","Cq")
#duplicate the sample column and split the sample type
type<-unlist(t(data.frame(strsplit(df_mean$sample,split = " "),stringsAsFactors = F))[,1])
df_mean$type<-type

#look at the biological variability Cp
ggplot(df_mean,aes(x=type,y=Cq))+geom_boxplot()+facet_wrap(~target)

#define the limit of detection of the assay
lod<-35
#transfor all the value in log2(EX)
df_mean$log2EX<-lod-df_mean$Cq

#look at the biological variability log2(EX)
ggplot(df_mean,aes(x=type,y=log2EX))+geom_boxplot()+facet_wrap(~target)

#subtract all the value of GAPDH from the tech
library(tidyr)
df_spread<-spread(data = df_mean[,-c(3,4)],key = target,value = log2EX)
df_spread$type<-unlist(t(data.frame(strsplit(df_spread$sample,split = " "),stringsAsFactors = F))[,1])
df_dct<-df_spread[2:8]-unlist(df_spread["GAPDH"])
df_dct_wide<-cbind.data.frame(df_dct,"sample"=df_spread$sample,"type"=df_spread$type)

#reform the long format and remove the GAPDH
df_dct_long<-gather(data = df_dct_wide,key = target,value = log2EX,-c(type,sample))
df_dct_long<-df_dct_long[df_dct_long$target!="GAPDH",]
#make the genes as factor
df_dct_long$target<-factor(df_dct_long$target,levels = c("TNF","CXCL10","CCL19","CCL22","MRC1","MSR1"))

#boxplot log2EX

ggplot(df_dct_long,aes(x=type,y=log2EX))+
  geom_boxplot()+
  labs(y="normalized log2(EX)")+
  facet_wrap(~target,scales = "free")

#make the stripplot

ggplot(df_dct_long,aes(x=type,y=log2EX))+
  geom_jitter()+
  labs(y="normalized log2(EX)")+
  facet_wrap(~target,scales = "free")

#manova
#transform the df in matrix
m<-as.matrix(df_dct_wide[,c("CXCL10","TNF","CCL19","MRC1","CCL22","MSR1")])
#define the factros
f<-factor(df_dct_wide$type)
man<-manova(m~f)
summary(man)
sink("manova.txt")
summary(man)
sink()
#the manova is significant test what variable contribute to this feature
summary(aov(m~f))
sink("anova.txt")
summary(aov(m~f))
sink()
#tukey comparison to see wich group is different from the others
library(purrr)
t<-colnames(m)
a<-map(t,function(x){
  TukeyHSD(aov(m[,x]~f))
})

names(a)<-t
sink("TukeyHSD.txt")
a
sink()

#compute the eta squared for each variable

library(purrr)
map(data.frame(m),function(x){
  y<-summary(aov(x~f))[[1]][2][[1]]
  y[1]/sum(y)
})

sink("eta_squared.txt")
map(data.frame(m),function(x){
  y<-summary(aov(x~f))[[1]][2][[1]]
  y[1]/sum(y)
})
sink()

#try to use the levelplot
library(lattice)
#I perform the normalization of each variable by the mean of expression variable wise and divide by its sd
r2<-apply(m,MARGIN = 2,FUN = mean)
s2<-apply(m,MARGIN = 2,FUN = sd)
m2<-t(m)
#m3<-t(m2-r2)
m4<-t((m2-r2)/s2)
#rownames(m3)<-f
rownames(m4)<-f
color=colorRampPalette(c("green","black","red"),space="rgb") 
#levelplot(m3,col.regions=color(100),xlab="samples",ylab="genes")
levelplot(m4,col.regions=color(100),xlab="samples",ylab="genes")


#PCA
library(ggfortify)
#maka a the wide format
autoplot(prcomp(df_dct[,c("CXCL10","TNF","CCL19","MRC1","CCL22","MSR1")]),data=df_spread,colour="type",size=5)
