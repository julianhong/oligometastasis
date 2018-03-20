#RPA.r
#imports oligomet database and reformats variables appropriately
#computes RPA and prints out outcomes data

#clear all the variables except for the data
rm(list = ls())

#import data
DB.FOR.ANALYSIS <- read.csv("DB20170418.csv", header = TRUE, check.names = TRUE)

#copy over as oligo
oligo<-DB.FOR.ANALYSIS

#appropriately format variables
oligo$OS..M. <- as.numeric(oligo$OS..M.)
oligo$DF.time <- as.numeric(oligo$DF.time)
oligo$LF.time.p.RT <- as.numeric(oligo$LF.time.p.RT)
oligo$LF <- as.numeric(oligo$LF)
oligo$Dead. <- as.numeric(oligo$Dead.)
oligo$Total.Dose <- as.numeric(oligo$Total.Dose)
oligo$X..lesions.treated <- as.numeric(oligo$X..lesions.treated)
oligo$DFS.time <- as.numeric(oligo$DFS.time)
oligo$M1.chem <- as.factor(oligo$M1.chem)
oligo$Gender <- as.factor(oligo$Gender)
oligo$LF.group <- as.factor(oligo$LF.group)
oligo$prior.curartive.local.tx <- as.factor(oligo$prior.curartive.local.tx)
oligo$BED75up <- as.factor(oligo$BED75up)
oligo$initial.lung <- as.factor(oligo$initial.lung)
oligo$initial.hilum..mediast <- as.factor(oligo$initial.hilum..mediast)
oligo$initial.liver <- as.factor(oligo$initial.liver)
oligo$initial.bone <- as.factor(oligo$initial.bone)
oligo$initial.adrenal <- as.factor(oligo$initial.adrenal)
oligo$abd.or.pelvic.LNs <- as.factor(oligo$abd.or.pelvic.LNs)
oligo$JH.primary.chemo <- as.factor(oligo$JH.primary.chemo)

#RPA w OS
library(rpart)
library(survival)
set.seed(220)
#include all variables to start
osrpa <- rpart(Surv(OS..M., Dead.) ~site + age + Gender + Dn..M1..M.+ organs.initial + 
                 prior.curartive.local.tx + M1.chem + JH.primary.chemo + 
                 X..lesions.treated + initial.lung + BED75up +
                 initial.hilum..mediast + initial.liver + initial.adrenal + initial.bone + 
                 abd.or.pelvic.LNs, data = oligo)
library(rpart.plot)
prp(osrpa) #print pre-parsed
 
#check to prune the tree
plotcp(osrpa) #print out plot to determine the appropriate complexity parameter
prunedos <- prune(osrpa, 0.018)
prp(prunedos)

oligo$osclass[1:nrow(oligo)] <- prunedos$where[1:nrow(oligo)]
oskm <-survfit(Surv(OS..M., Dead.) ~ osclass, oligo)
library(survminer)
ggsurvplot(oskm, xlab="Time (months)", ylab="Overall survival", pval = TRUE, 
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           break.time.by = 12, risk.table = TRUE, legend = "right", 
           palette = c("#00427C","#E9B200", "firebrick3", "forestgreen", "darkorchid"),
           legend.labs = c("Class 1", "Class 2", "Class 3", "Class 4", "Class 5"), 
           linetype = "strata", xlim = c(0,70))
summary(oskm,times=c(12*3))
#log-rank
survdiff(Surv(OS..M., Dead.) ~oligo$osclass, oligo, rho = 0)

#cox to reassess the #s
oligo$osclass <- as.factor(oligo$osclass)
osclasscox <- coxph(Surv(OS..M., Dead.)~osclass, data = oligo)
summary(osclasscox)

########
#RPA with DFS
library(rpart)
library(survival)
set.seed(220)
pfsrpa <- rpart(Surv(DFS.time, DFS) ~site + age + Gender + Dn..M1..M.+ organs.initial + 
                  prior.curartive.local.tx + M1.chem + JH.primary.chemo + 
                  X..lesions.treated + initial.lung + BED75up +
                  initial.hilum..mediast + initial.liver + initial.adrenal + initial.bone + 
                  abd.or.pelvic.LNs, data = oligo)

library(rpart.plot)
prp(pfsrpa) #print pre-parsed
plotcp(pfsrpa)
prunedpfs <- prune(pfsrpa, 0.038)
prp(prunedpfs)

#exclude patients with missing PFS data in generating curves
cleanpfs <- oligo[c("DFS.time", "DFS")]
cleanpfs <- na.omit(cleanpfs)

cleanpfs$pfsclass[1:nrow(cleanpfs)] <- prunedpfs$where[1:nrow(cleanpfs)]
pfskm <-survfit(Surv(DFS.time, DFS) ~ pfsclass, cleanpfs)
library(survminer)
ggsurvplot(pfskm, xlab="Time (months)", ylab="Progression-free survival", break.time.by = 12, 
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           risk.table = TRUE, legend = "right", palette = c("#00427C","#E9B200"), 
           pval = TRUE, xlim = c(0,70), linetype = "strata", 
           legend.labs = c("Class 1", "Class 2"))
summary(pfskm,times=c(12*3))
#log-rank
survdiff(Surv(DFS.time, DFS) ~cleanpfs$pfsclass, cleanpfs, rho = 0)

#cox to reassess the #s
cleanpfs$pfsclass <- as.factor(cleanpfs$pfsclass)
pfsclasscox <- coxph(Surv(DFS.time, DFS)~pfsclass, data = cleanpfs)
summary(pfsclasscox)

################
#RPA w TMC
library(rpart)
library(survival)
set.seed(220)
tmcrpa <- rpart(Surv(LF.time.p.RT, LF) ~site + age + Gender + histo + Dn..M1..M.+ organs.initial + 
                  prior.curartive.local.tx + M1.chem + JH.primary.chemo + 
                  X..lesions.treated + BED75up + initial.CNS + initial.lung + 
                  initial.hilum..mediast + initial.liver + initial.adrenal + initial.bone + 
                  abd.or.pelvic.LNs, data = oligo)

library(rpart.plot)
prp(tmcrpa) #print pre-parsed

#check to prune the tree
plotcp(tmcrpa) #no local minimum

##########
#Kaplan-Meier analyses
#whole cohort OS
allos.surv <-survfit(Surv(oligo$OS..M., oligo$Dead.) ~1)
ggsurvplot(allos.surv, xlab="Time (months)", ylab="Overall survival", 
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           break.time.by = 12, risk.table = TRUE, color = c("#00427C"),
           xlim = c(0,70))
summary(allos.surv,times=c(12*3))
allos.surv

#whole cohort PFS
allpfs.surv <-survfit(Surv(oligo$DFS.time, oligo$DFS) ~1)
ggsurvplot(allpfs.surv, xlab="Time (months)", ylab="Progression-free survival", 
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           break.time.by = 12, risk.table = TRUE, color = c("#00427C"),
           xlim = c(0,70))
summary(allpfs.surv,times=c(12*3))
summary(allpfs.surv,times=c(12*5))
allpfs.surv

#whole cohort TMC
alltmc.surv <-survfit(Surv(oligo$LF.time.p.RT, oligo$LF) ~1)
ggsurvplot(alltmc.surv, xlab="Time (months)", ylab="Control of all treated metastases", 
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           break.time.by = 12, risk.table = TRUE, color = c("#00427C"),
           xlim = c(0,70))
summary(alltmc.surv,times=c(12*3))
summary(alltmc.surv,times=c(12*5))
alltmc.surv

#whole cohort OS by BED75
osbed <-survfit(Surv(OS..M., Dead.) ~ BED75up, oligo)
ggsurvplot(osbed, xlab="Time (months)", ylab="Overall survival", pval = TRUE,
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           break.time.by = 12, risk.table = TRUE, legend = "right", 
           palette = c("#00427C","#E9B200"),
           legend.labs = c("<75 BED", 
                           ">=75BED"),
           linetype = "strata", xlim = c(0,70))
survdiff(Surv(OS..M., Dead.) ~ BED75up, oligo, rho = 0)
summary(osbed,times=c(12*3))
summary(osbed,times=c(12*5))
osbed

#whole cohort pfs by BED75
pfsbed <-survfit(Surv(DFS.time, DFS) ~ BED75up, oligo)
ggsurvplot(pfsbed, xlab="Time (months)", ylab="Progression-free survival", pval = TRUE,
           font.x = c(24, "plain", "black"), font.y = c(24, "plain", "black"),
           font.tickslab = c(20, "plain", "black"), font.legend = c(20,"plain", "black"),
           risk.table.fontsize = 6, font.main = c(20, "plain", "black"),
           break.time.by = 12, risk.table = TRUE, legend = "right", 
           palette = c("#00427C","#E9B200"),
           legend.labs = c("<75 BED", 
                           ">=75BED"),
           linetype = "strata", xlim = c(0,70))
survdiff(Surv(DFS.time, DFS) ~ BED75up, oligo, rho = 0)
summary(pfsbed,times=c(12*3))
summary(pfsbed,times=c(12*5))
pfsbed