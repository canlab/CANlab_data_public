


#https://cran.r-project.org/web/packages/robustBLME/robustBLME.pdf 
https://github.com/brentp/rstuff/blob/master/rmodels.R
file:///home/pavel/Downloads/welchADF.pdf

#----confidance intervals for rlmerMod object 
confint.rlmerMod <- function(object,parm,level=0.95) {
  beta <- fixef(object)
  if (missing(parm)) parm <- names(beta)
  se <- sqrt(diag(vcov(object)))
  z <- qnorm((1+level)/2)
  ctab <- cbind(beta-z*se,beta+z*se)
  colnames(ctab) <- stats:::format.perc(c((1-level)/2,(1+level)/2),
                                        digits=3)
  return(ctab[parm,])
}


#calculates p-value for rlmer model
pvalue.rlmerMod <- function(object)  {
  
  2 * (1 - pnorm(abs(fixef(object)/sqrt(diag(vcov(object))))))  
}

#calculates Cohen's D effect size with confidance intervals 95% for lmer effects
cohenD.lmerMod <- function(object)  {
  
  sum_lmer<-data.frame(coef(summary(object)))
  sum_lmer$D<-sum_lmer$t.value/sqrt(sum_lmer$df)
  sum_lmer$Dlow<-(sum_lmer$Estimate- qt(c(.975), df=sum_lmer$df)*sum_lmer$Std..Error)/sum_lmer$Std..Error/sqrt(sum_lmer$df)
  sum_lmer$Dhigh<-(sum_lmer$Estimate+ qt(c(.975), df=sum_lmer$df)*sum_lmer$Std..Error)/sum_lmer$Std..Error/sqrt(sum_lmer$df)
  
  return(sum_lmer)
}



#----reading the movement data------------
setwd("/home/pavel/Projects/MEA_doc_pat") #setup your working directory here
load("data/combined2_1.RData")

save(agg_all_data, file="/home/pavel/Projects/MEA_doc_pat/data/agg_all_data_100219.RData")
load("data/agg_all_data_100219.RData")
names(agg_all_data)

library(zoo)
library(package = "hydromad")
library(tidyr)
library(pracma)


#----plot of the rolling correlation----------
library(lattice)


#---------------------Figure 3--------------------

i="1M_GDGP"
dlog<-(combined2$V1+1)[combined2$id==i]
plog<-(combined2$V2+1)[combined2$id==i]
for_an<-cbind(dlog,plog)
for_an1<-for_an[round(nrow(for_an)*0.1):round(nrow(for_an)*0.9),]
x<-rollccf(for_an1,width = list( 300),by=300)
names(x$rolls)<-"running correlation"
names(x$data)<-c("clinician movement","patient movement")
c_pl<-xyplot(x, main="Concordant dyad",xlab="Interaction time (sec)", scales=list(x=list(at=c(0,1800,3600,5400,7200,9000,10800,12600,14400), labels=c(0,60,120,180,240,300,360, 420,480)))
             ,, panel=function(...) {
               panel.xyplot(...)
               panel.abline(h=.0)
             })


j="7M_YDGP"
dlog2<-(combined2$V1+1)[combined2$id==j]
plog2<-(combined2$V2+1)[combined2$id==j]
for_an2<-cbind(dlog2,plog2)
for_an12<-for_an2[round(nrow(for_an2)*0.1):round(nrow(for_an2)*0.9),]
y<-rollccf(for_an12,width = list( 300),by=300)
names(y$rolls)<-"running correlation"
names(y$data)<-c("clinician movement","patient movement")
d_pl<-xyplot(y, main="Discordant dyad",xlab="Interaction time (sec)", scales=list(x=list(at=c(0,1800,3600,5400,7200,9000,10800,12600,14400,16200), labels=c(0,60,120,180,240,300,360, 420,480,540)))
             ,, panel=function(...) {
               panel.xyplot(...)
               panel.abline(h=.0)
             })
tiff("/home/pavel/Projects/MEA_doc_pat/Figures/Figure 3.tiff", width = 10, height = 10, units = 'in', res = 300)
png("/home/pavel/Projects/MEA_doc_pat/Figures/Figure 3.png", width = 7, height = 7, units = 'in', res = 300)

print(c_pl, position = c(0, 0.5, 0.7, 1), more = TRUE)
print(d_pl, position = c(0, 0, 0.7, 0.5), more = F)
panel.text(40, 50, "A.", cex = 1.5, font = 2)
panel.text(40, 1100, "B.", cex = 1.5, font = 2)

#panel.text(22, 18, "A.", cex = 1.5, font = 2)
#panel.text(22, 350, "B.", cex = 1.5, font = 2)

dev.off() 

#------computing correlations between movement time series 
og<-numeric()
plog<-numeric()
for_an<-numeric()
corr_data_run<-numeric()
corr_data<-numeric()
corr<-numeric()
doc_move<-numeric()
pat_move<-numeric()
library(devtools)
library(plotflow)
table(combined2$id)

library(lattice)

pdf(file="/home/pavel/Projects/MEA_doc_pat/output/plots_corr_011618.pdf")
for (i in levels(as.factor(combined2$id)) ) {
  dlog<-scale(log10(combined2$V1+1))[combined2$id==i]
  plog<-scale(log10(combined2$V2+1))[combined2$id==i]
  
  attr(dlog,"scaled:center")<-NULL 
  attr(dlog,"scaled:scale")<-NULL 
  attr(plog,"scaled:center")<-NULL 
  attr(plog,"scaled:scale")<-NULL 
  for_an<-cbind(dlog,plog)
  for_an1<-for_an[round(nrow(for_an)*0.1):round(nrow(for_an)*0.9),]
  
  #----10 seconds rolling window------
  x<-rollccf(for_an1,width = list( 300),by=150)
  
  
print(xyplot(x, main=paste0(i), panel=function(...) {
   panel.xyplot(...)
 panel.abline(h=.0)
 }) )
  
  #-----window-level  data
  corr<-data.frame(x$rolls)
  corr$m_corr<-ave(corr[,1],corr[,2],corr[,3])
  time <- seq.int(nrow(corr))
  
  #------------amount of movements--------------
  doc_move<-mean(combined2$V1[combined2$id==i])
  pat_move<-mean(combined2$V2[combined2$id==i])
  doc_move_sd<-sd(combined2$V1[combined2$id==i])
  pat_move_sd<-sd(combined2$V2[combined2$id==i])
  
  
  corr_data<-data.frame(corr$m_corr, i,mean(x$lags),doc_move,pat_move,doc_move_sd,pat_move_sd,time)
  
  corr_data_run<-rbind(corr_data,corr_data_run)
  
  summary(corr_data)
  
  
}
dev.off() 


library(psych)

#---names +Fisher Z
names(corr_data_run)<-c("runcorr", "dyad", "max_lag","doc_move","pat_move","doc_move_sd","pat_move_sd","time")

corr_data_run$ind<-with(corr_data_run, ifelse(corr_data_run$runcorr>0.3,  1, 0))
corr_data_run$runcorrz<-fisherz(corr_data_run$runcorr)


library(ggplot2)
corr_data_run1<-separate(corr_data_run, dyad, c("group","cond"), sep = "_")
corr_data_run1$congruent<-with(corr_data_run1, ifelse(corr_data_run1$cond=="YDYP" | corr_data_run1$cond=="GDGP",  1, 0))
corr_data_run1$congruent1 <- factor(corr_data_run1$congruent, levels = c(0,1),labels = c("Discordant", "Concordant"))

tiff("/home/pavel/Projects/MEA_doc_pat/output/runningcorr_011718.tiff", width = 10, height = 10, units = 'in', res = 300)

h <- ggplot(corr_data_run1, aes(time, runcorrz,colour =congruent1 ),na.rm = TRUE)
h  +  stat_smooth(fill = "grey50", span = 0.5,size = 2, alpha = 0.5)+theme(text = element_text(size=18),legend.title = element_blank())+xlab("Running window") +
  ylab("Fisher-Z correlation") 
dev.off()



#------ dyad-level data---------------,

dyad_corr<-aggregate(.~dyad, corr_data_run, FUN = quantile, probs  = 0.80)
dyad_corr<-aggregate(.~dyad, corr_data_run, FUN =  mean)
dyad_corr1<-separate(dyad_corr, dyad, c("group","cond"), sep = "_")
dyad_corr1$doc<-substr(dyad_corr1$cond, start = 1, stop = 2)
dyad_corr1$pat<-substr(dyad_corr1$cond, start = 3, stop = 4)
dyad_corr1$pat_id<-paste0(dyad_corr1$group,"_",dyad_corr1$pat)
dyad_corr1$doc_id<-paste0(dyad_corr1$group,"_",dyad_corr1$doc)
dyad_corr1$pat_id_char<-dyad_corr1$pat_id
dyad_corr1$doc_id_char<-dyad_corr1$doc_id

dyad_corr1$congruent<-with(dyad_corr1, ifelse(dyad_corr1$cond=="YDYP" | dyad_corr1$cond=="GDGP",  1, 0))
dyad_corr1$max_lag_b<-with(dyad_corr1, ifelse(dyad_corr1$max_lag>0,  1, 0))



#Testing path a


library(lme4)
sn<-lmer(runcorrz ~ congruent+pat_move +pat_move_sd+(1 | pat_id)+(1 | doc_id), data=dyad_corr1)
summary(sn)

library(r2glmm)
r2.m1 = r2beta(sn, method = 'nsj', partial = T)



tiff(filename="/home/pavel/Projects/MEA_doc_pat/output/runningcorr_082818.tiff",  width=1000, height=800, res=200)

p <- ggplot(dyad_corr1, aes(factor(congruent), runcorrz))
p<-p+geom_violin(fill = "pink", colour = "#3366FF")
p<-p+geom_point(size=2,colour="#000000",position=position_jitter(width=0.025))
p<-p+geom_line(aes(group=pat_id))
p<-p+ scale_x_discrete(breaks=c("0", "1"),labels=c("Discordant", "Concordant"))
p<-p+xlab("")
p<-p+ylab("Movement synchrony")
p<-p+theme_bw()
p<-p+theme(axis.text.x = element_text(colour="black",size=14))
p<-p+theme(axis.title.y = element_text(colour="black",size=14))
p
dev.off()


#--------pain mediation model-----------------
install.packages("MC")
library(welchADF)
library(robustlmm)
library(lmerTest)
library(emmeans)
library(WRS2)
library(pbkrtest)
library(r2glmm)
library(ggplot2) 
library(MASS)
library(ggpubr)

#unloadNamespace("lmerTest")
names(agg_all_data)


agg_all_data$trust<-(agg_all_data[,41]+agg_all_data[,42])/2 #---mean of two trust measurements


#---control variables of movement
corr_move<- lmer(runcorrz~+pat_move+ doc_move+ (1|doc_id)+(1|pat_id),data=agg_all_data)
summary(corr_move)
cohenD.lmerMod(corr_move)
confint(corr_move)
anova(corr_move)


pain_move<- lmer(CP_PEAK~+pat_move+ doc_move+ (1|doc_id)+(1|pat_id),data=agg_all_data)
summary(pain_move)
cohenD.lmerMod(pain_move)
confint(pain_move)
anova(pain_move)

trust_move<- lmer(trust~+pat_move+ doc_move+ (1|doc_id)+(1|pat_id),data=agg_all_data)
summary(trust_move)
cohenD.lmerMod(trust_move)
confint(trust_move)
anova(trust_move)



#----Mediation pain--------

pain.dir1 <- lmer(CP_PEAK~ congruent+ pat_move+  (1|pat_id)+  (1|doc_id),data=agg_all_data,REML = FALSE, method="DASvar")
summary(pain.dir1)
cohenD.lmerMod(pain.dir1)
anova(pain.dir1)
confint(pain.dir1)


pain.med1 <- rlmer(runcorrz~congruent+trust +pat_move+ (1|doc_id)+(1|pat_id),data=agg_all_data)
summary(pain.med1)
confint(pain.med1)
pvalue.rlmerMod(pain.med1)

pain.med1_1 <- lmer(runcorrz~congruent +pat_move+(1|pat_id),data=agg_all_data)
anova(pain.med1_1)
summary(pain.med1_1)
cohenD.lmerMod(pain.med1_1)
confint(pain.med1_1)

pain.med2 <- rlmer(CP_PEAK~ congruent+trust+runcorrz  +pat_move+ (1|doc_id)+ (1|pat_id),data=agg_all_data)
summary(pain.med2)
confint(pain.med2)
pvalue.rlmerMod(pain.med2)

pain.med2_1 <- lmer(CP_PEAK~ congruent +pat_move+ +runcorrz+(1|pat_id),data=agg_all_data)
anova(pain.med2_1)
summary(pain.med2_1)
cohenD.lmerMod(pain.med2_1)
confint(pain.med2_1)

r2beta(pain.med2_1, method = 'nsj', partial = T)


## robust testing of mediating effect (indirect effect)
#with(agg_all_data, ZYmediate(congruent, CP_PEAK, runcorrz,nboot = 10000))

library(mediation)
med.pain <- mediate(pain.med1_1 , pain.med2_1, treat = "congruent", mediator = "runcorrz", sims = 5000,dropobs=T )
summary(med.pain)



#--------------mediation model for trust----------------

trust.dir <- rlmer(trust~ CP_PEAK+congruent+ pat_move+ (1|doc_id)+ (1|pat_id),data=agg_all_data,REML = FALSE, method="DASvar")
summary(trust.dir)
confint(trust.dir)
pvalue.rlmerMod(trust.dir)

trust.dir1 <- lmer(trust~ congruent+ pat_move+  (1|pat_id),data=agg_all_data,REML = FALSE, method="DASvar")
summary(trust.dir1)
cohenD.lmerMod(trust.dir1)
confint(trust.dir1)
anova(trust.dir1)

summary(agg_all_data$trust)

trust.med <- rlmer(trust~ CP_PEAK+congruent+runcorrz+ pat_move+ (1|doc_id)+ (1|pat_id),data=agg_all_data)
summary(trust.med)
confint(trust.med)
pvalue.rlmerMod(trust.med)

trust.med1 <- lmer(trust~ congruent+runcorrz+ pat_move+  (1|pat_id),data=agg_all_data)
summary(trust.med1)
cohenD.lmerMod(trust.med1)
confint(trust.med1)
anova(trust.med1)
r2beta(trust.med1, method = 'nsj', partial = T)

with(agg_all_data, ZYmediate(congruent, trust, runcorrz,nboot = 10000))
fixef(trust.med)

med.trust <- mediate(pain.med1_1 , trust.med1, treat = "congruent", mediator = "runcorrz", sims = 5000,dropobs=T )
summary(med.trust)


library(mediation)

#----Plots----------------

tiff(filename="/home/pavel/Projects/MEA_doc_pat/output/runningcorr_082818.tiff",  width=1000, height=800, res=200)

p <- ggplot(agg_all_data[!is.na(agg_all_data["runcorrz"]),], aes(factor(congruent), runcorrz))
p<-p+geom_violin(fill = "pink", colour = "#3366FF")
p<-p+geom_point(size=2,colour="#000000",position=position_jitter(width=0.025))
p<-p+geom_line(aes(group=pat_id))
p<-p+ scale_x_discrete(breaks=c("0", "1"),labels=c("Discordant", "Concordant"))
p<-p+xlab("")
p<-p+ylab("Movement synchrony")
p<-p+theme_bw()
p<-p+theme(axis.text.x = element_text(colour="black",size=14))
p<-p+theme(axis.title.y = element_text(colour="black",size=14))
p
dev.off()


#tiff("/home/pavel/Projects/MEA_doc_pat/output/sync_pain_trust_180618.tiff", width = 10, height = 10, units = 'in', res = 300)

png("/home/pavel/Projects/MEA_doc_pat/output/sync_pain_trust_060918.png", width = 10, height = 10, units = 'in', res = 300)

p1<-ggplot(data = agg_all_data, aes(x = runcorrz , y = CP_PEAK)) +
  geom_point( size = 2) +
  geom_smooth(method="rlm")+theme(text = element_text(size=18),legend.title = element_blank())+ylab("Pain rating") +xlab("Movement synchrony")
p1 + annotate("text", x = -0.1, y = 80, label = "Some text")

p2<-ggplot(data = agg_all_data, aes(x = runcorrz , y = trust)) +
  geom_point( size = 2) +
  geom_smooth(method="rlm")+theme(text = element_text(size=18),legend.title = element_blank())+ylab("Trust in clinician") +xlab("Movement synchrony")

S
#names(agg_all_data)
dev.off()

plot(agg_all_data$runcorrz,agg_all_data$trust)
names(agg_all_data)
