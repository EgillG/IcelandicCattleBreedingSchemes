# Script to analyse results of simulations
library(tidyverse)
rm(list=ls())
# set the proper working directory
setwd("C:/Users/au589863/OneDrive - Aarhus universitet/Documents/phdProject4_Simulations/MoBPS/Results/46Test/")
matrices = c("M1","M2","M2R","M2D","M1D","M1R","M1_05","Ped", "M1D_Est")

df = NULL
for (scenario in matrices) {
  for (rep in c(1:10)) {
    nData=read.table(paste(scenario,"_",rep,"_info.txt", 
                           sep = ""))
    nData=round(nData,4)
    nData$Gen <- as.numeric(row.names(nData))
    nData$Scen = scenario
    nData$Rep = rep
    sd <- nData[nData$Gen==1,]$GeneticVarA**0.5
    # nData$StBV <- (nData$BV-mean(nData$BV))/sd(nData$BV)
    # nData$StBV <- (nData$BV-min(nData$BV))/(max(nData$BV)-min(nData$BV))*100+100
    nData$stBV <- (nData$BV-nData[nData$Gen==1,]$BV)/sd
    #    nData$stBV <- (nData$BV-min(nData$BV,na.rm = T))/sd(nData$BV,na.rm = T)
    #nData$stBV <- 100+(nData$BV-min(nData$BV,na.rm = T))/sd(nData$BV,na.rm = T)*10
    nData$stGenicVarA <- (nData$varA-max(nData$varA,na.rm = T))/sd(nData$varA,na.rm = T)
    nData$stGenicVarA <- 100+(nData$varA-max(nData$varA,na.rm = T))/(sd(nData$varA,na.rm = T))*10
    nData$stGenicVarA <- 1+(nData$varA-nData[nData$Gen==6,]$varA)/nData[nData$Gen==6,]$varA
    
    nData$stGeneticVarA <- (nData$GeneticVarA-max(nData$GeneticVarA,na.rm = T))/sd(nData$GeneticVarA,na.rm = T)
    nData$stGeneticVarA <- 100+(nData$GeneticVarA-max(nData$GeneticVarA,na.rm = T))/(sd(nData$GeneticVarA,na.rm = T))*10
    nData$stGeneticVarA <- 1+(nData$GeneticVarA-nData[nData$Gen==6,]$GeneticVarA)/nData[nData$Gen==6,]$GeneticVarA
    df=rbind(df,nData)
  }
}

df <- df[,c(16,17,18,1,19,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,21)]
#sd <-sd(df[df$Gen==1,]$BV)
head(df)
tail(df)
colnames(df)
df$Scen <- as.factor(df$Scen)
df$SelfRel <- (df$SelfRel-0.5)*2
#Subset data for regression:
dfRegs <- df[df$Gen>5,]
lm1 <- lm(stBV~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm2 <- lm(-log(1-Kinship)~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm3 <- lm(stGeneticVarA~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm4 <- lm(stGenicVarA~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm5 <- lm(-log(1-F_drift_n)~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm6 <- lm(-log(1-F_hom_n)~as.factor(Rep)+Gen:Scen, data = dfRegs)

coefficients <- cbind(lm1$coefficients,lm2$coefficients,lm3$coefficients,lm4$coefficients,lm5$coefficients,lm6$coefficients)
colnames(coefficients) <- c("stBV","kinship","stVarA","stGenicVarA","F_drift_n","F_hom_n")
coefficients <- as.data.frame(coefficients)
coefficients
cotemp <- round(tail(coefficients,length(matrices)),6)[c(3,4,7,1,6,5,8,2,9),c(2,1,3,4)]
cotemp <- cbind(cotemp, c(cotemp$kinship/cotemp$kinship[9]), cotemp$stBV/cotemp$stBV[9], cotemp$stVarA/cotemp$stVarA[9], cotemp$stGenicVarA/cotemp$stGenicVarA[9])
cotemp$kinship <- cotemp$kinship*100
cotemp$stVarA <- cotemp$stVarA*100
cotemp$stGenicVarA <- cotemp$stGenicVarA*100
colnames(cotemp) <- c("Kinship", "BV","GeneticVar","GenicVar","%Kinship","%BV","%GeneticVar","%GenicVar")
write.table(tail(cotemp,length(matrices)), file = "Table1.txt", sep = "\t")

res_maf <- read.table("ResultsMaf/homAll0.001.txt", header = T)
dfRegs <- res_maf[res_maf$Gen>5,]
lm1 <- lm(-log(1-F_hom_qtl)~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm2 <- lm(-log(1-F_hom_n)~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm3 <- lm(-log(1-F_drift_qtl)~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm4 <- lm(-log(1-F_drift_n)~as.factor(Rep)+Gen:Scen, data = dfRegs)
coefficients <- cbind(lm1$coefficients,lm2$coefficients,lm3$coefficients,lm4$coefficients)
cotemp <- round(tail(coefficients,length(matrices)),6)[c(3,4,7,1,6,5,8,2,9),]
cotemp <- cotemp*100
colnames(cotemp) <- c("F_hom_qtl","F_hom_n","F_drift_qtl","F_drift_n")
write.table(cotemp,"F_hom0.001.txt")

res_maf <- data.frame(res_maf)
agg <- aggregate(res_maf[res_maf$Gen==max(res_maf$Gen),c(4:7)], by = list(res_maf[res_maf$Gen==max(res_maf$Gen),]$Scen), FUN = "mean")
write.table(agg, "LastGen0.001.txt", sep = "\t")
df <- cbind(df[order(df$Scen),c(1,2,3,4,5,6,7,8,9,10,11,18,19,20,21)],res_maf[order(res_maf$Scen),4:7])

# Now go for the MAF filtered F_hom and F_qtl for markers.

res_maf <- read.table("ResultsMaf/homMarkers0.001.txt", header = T)
dfRegs <- res_maf[res_maf$Gen>5,]
lm1 <- lm(-log(1-F_hom_m)~as.factor(Rep)+Gen:Scen, data = dfRegs)
lm2 <- lm(-log(1-F_drift_m)~as.factor(Rep)+Gen:Scen, data = dfRegs)
coefficients <- cbind(lm1$coefficients,lm2$coefficients)
colnames(coefficients) <- c("F_hom_markers","F_drift_markers")
coefficients <- tail(coefficients,9)[c(3,4,7,1,6,5,8,2,9),]
coefficients <- coefficients * 100
write.table(coefficients,"F_homMarkers0.001.txt")
res_maf <- data.frame(res_maf)
agg <- aggregate(res_maf[res_maf$Gen==max(res_maf$Gen),c(4:5)], by = list(res_maf[res_maf$Gen==max(res_maf$Gen),]$Scen), FUN = "mean")
agg <- agg[c(3,4,7,1,6,5,8,2,9),]
write.table(agg, "LastGenMarkers0.001.txt", sep = "\t")
df <- cbind(df[order(df$Scen),],res_maf[order(res_maf$Scen),4:5])
##################################################

library(dplyr)
Old <- c("M1","M1D","M1D_Est","M1R","Ped", "M2","M2D", "M2R", "M1_05")
New <- c("VR1 All","VR1 Base","VR1 Old","VR1 Current","Pedigree","VR2 All","VR2 Base","VR2 Current", "VR1 0.5")

df$Scen <- recode(df$Scen, !!!setNames(New,Old))

# Contrasts for P-values
# df$Scen <- relevel(df$Scen, ref = "VR1 Base")
# df$Scen <- relevel(df$Scen, ref = "VR2 Base")

mat <- matrix(NA, nrow = 9, ncol = 9)
colnames(mat) <-   unique(as.character(df$Scen))
rownames(mat) <-  unique(as.character(df$Scen))
df <- df[order(as.character(df$Scen)),]
lm1 <- lm(stBV~as.factor(Rep)+Gen:Scen+Gen, data = df[df$Gen>5,])
summary(lm1)
lm1 <- lm(-log(1-F_drift_m)~as.factor(Rep)+Gen:Scen, data = df[df$Gen>5,])
summary(lm1)

# Remember to regress on -log(1-Kinship also)
for (i in 1:9) {
  nafn <- paste(unique(as.character(df$Scen))[i])
  print(nafn)
  df$Scen <- relevel(df$Scen, ref = nafn)
  lm1 <- lm(SegAlleles*3000~as.factor(Rep)+Gen:Scen+Gen, data = df[df$Gen>5,])
  summary(lm1)
  # Estimate
  # print(cbind(as.character(sort(unique(as.character(df$Scen))))[-i],lm1$coefficients[(length(lm1$coefficients)-length(unique(df$Scen))+2):length(lm1$coefficients)]))
  # P value
  print(as.matrix(summary(lm1)$coefficients[(length(lm1$coefficients)-length(unique(df$Scen))+2):length(lm1$coefficients),4]))
#  mat[i,c(1:9)[-i]] <-  tail(summary(lm1)$coefficients,8)[,4]
}

df1 <- df %>% group_by(Gen, Scen) %>% summarise(
  BV = mean(BV),
  stBV = mean(stBV),
  Kinship = mean(Kinship),
  Homo = mean(1-Heteroz),
  SegMarkers = mean(SegAlleles),
  SegQTL = mean(SegQTL),
  SegNeutral = mean(SegNeutral),
  F_drift_neutral = mean(F_drift_n),
  F_hom_neutral = mean(F_hom_n),
  F_drift_qtl = mean(F_drift_qtl),
  F_hom_qtl = mean(F_hom_qtl),
   F_drift_markers = mean(F_drift_m),
   F_hom_markers = mean(F_hom_m),
  VarA = mean(varA),
  stVarA = mean(stGenicVarA),
  GenetVarA = mean(GeneticVarA),
  stGenetVarA = mean(stGeneticVarA),
  Inbreeding = mean(SelfRel),
  GeneticVarA = mean(GeneticVarA)
)

df1 <- df1[df1$Gen>5,]
df1$Gen <- df1$Gen - 1

pl <- ggplot(data = df1, aes(x = -log(1-Kinship), y = stBV, colour = Scen, linetype = Scen)) +
  geom_line() +
  theme_minimal() +
  ylab("Breeding value") +
  theme(legend.position="top") +
  labs(color = "Scenario", linetype = "Scenario")
pl
ggsave("BV_Kinship.png", height = 3.5, width = 4.5, dpi = 1800) 
ggsave("BV_Kinship.png",height = 6, width = 7, dpi = 1800) 

longdf <- gather(df1, method, value, F_hom_neutral, Kinship,F_drift_neutral)

pl <- ggplot(data = longdf, aes(x = -log(1-value), y = stBV, colour = Scen, linetype = Scen)) +
  geom_line() +
  facet_wrap(~method, scales = "free_x", ncol=3) +
  theme_minimal() +
  ylab("Genetic gain") +
  xlab("-log(1-F)") +
  theme(legend.position="top") +
  labs(color = "Scenario", linetype = "Scenario")
pl
ggsave("BV_Kinship_facets.png",height = 4, width = 7, dpi = 600) 
