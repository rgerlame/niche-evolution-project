setwd("C:/Users/rgermain/Dropbox/2 - Manuscripts in prep/Niche evolution project/Trait expt")

require(lmerTest);require(visreg);require(ggplot2);require(car); require(vegan)

###############
#data prep
###############

df<-read.csv("trait experiment - parallel.csv")
df<-subset(df,origin!="sp")
rownames(df)<-paste(df$pop,df$id,sep="-")

df$pop<-factor(df$pop,levels=c("6","19","26","3","20","29","11","5","27","13","16","4","15","10","2","18","23","25","17","9","S2","S31","S33"))
envt<-read.csv("env.means.all.csv")
df.envt<-merge(df,envt,by.x="pop",by.y="site")
df.envt$origin2<-df.envt$origin
df.envt$origin2<-ifelse(df.envt$origin=="s","s","a")

df.seed<-read.csv("seed mass.csv")
seed.df<-merge(df.seed,df.envt[,c(1:5,26,27)],by="pop",all.x=TRUE,all.y=FALSE)
#seed.df<-subset(seed.df,origin!="ad")

#############################
#making composite trait axes
#############################

df.pca<-df

df.pca$est.height<-log(df.pca$est.height)
df.pca$ratio<-log(df.pca$ratio)
df.pca$ratio<-1
df.pca$avg<-log(df.pca$avg)

df.sub<-na.omit(df.pca)
df.sub$label<-paste(df.sub$species,df.sub$pop,sep="")

pca<-prcomp(df.sub[,c(6:8,10)],scale.=TRUE,center=TRUE)
biplot(pca,xlabs=df.sub$label)
pca$x
df.sub$pca1<-pca$x[,1]
df.sub$pca2<-pca$x[,2]

###############
#basic stats
###############

#legend: contrast = which of the 8 Bromus populations the sympatric (s), allopatric from same environment (as), and allopatric from different environment (ad) Vulpia populations are competed against
#two options that I've tried but turned off -- feel free to try: adding in random slopes ((origin|contrast/pop)) and adding in origin*poly(log(productivity+1),2) to test if how productive the site a given Vulpia population originates from affects how strongly coevolutionary history matters
#main result: generally no difference between sympatric and allopatric-same environment treatments, signfificant effect is mostly driven by allopatric-different (aka evolution to abiotic environment)

#seed size
lm1<-lmer(log(seed.wt) ~ origin + (1|contrast/pop),data=subset(seed.df,species.x=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin")

#germ time
lm1<-lmer(log(avg) ~ origin + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#prop germ, more dormancy? 
lm1<-lmer(tot.germ ~ origin + (1|contrast/pop),data=subset(df,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#this model is same response variable as above but has productivity added just to show how it might matter...
lm1<-lmer(tot.germ~origin*log(productivity) + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="productivity",by="origin",overlay=T)

#root depth
lm1<-lmer(root.depth ~ origin + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#emg height
lm1<-lmer(log(emg.height) ~ origin + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#est height
lm1<-lmer(log(est.height) ~ origin + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#ratio
lm1<-lmer(log(ratio) ~ origin + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#flowering
lm1<-lmer(flowering~origin + (1|contrast/pop),data=subset(df.envt,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#pca1
lm1<-lmer(pca1~origin + (1|contrast/pop),data=subset(df.sub,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)

#pca2
lm1<-lmer(pca2~origin + (1|contrast/pop),data=subset(df.sub,species=="vm"))
Anova(lm1)
visreg(lm1,xvar="origin",overlay=T)


##########
#plots
##########

#plotting grouped by 'contrast'
p <- ggplot(subset(df,origin!="sp"), aes(x=origin, y=flowering,color=species)) + geom_boxplot() 
p + facet_wrap( ~ contrast, nrow = 1)

p <- ggplot(subset(df.sub,origin!="sp"), aes(x=origin, y=pca1,color=species)) + geom_boxplot() 
p + facet_wrap( ~ contrast, nrow = 1)

p <- ggplot(subset(df.sub,origin!="sp"), aes(x=origin, y=pca2,color=species)) + geom_boxplot() 
p + facet_wrap( ~ contrast, nrow = 1)



