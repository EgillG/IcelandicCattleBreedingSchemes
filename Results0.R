source("Functions_MoBPS_simulations_20_01_2022.R")
library("tidyverse")
hom <- matrix(nrow=19,ncol = 7)
colnames(hom) <- c("Rep","Scen","Gen","F_hom_qtl", "F_drift_qtl","F_hom_n", "F_drift_n")
#matrices <- c("M2D", "M1D")
matrices = c("M1","M2","M2R","M2D","M1D","M1R","M1_05","Ped")

#for (MAF in c(0, 0.001,0.01,0.05)) {
for (MAF in "zeroCheck") {
homAll = NULL
setwd("Results")
for (scenario in matrices) {
  for (rep in 1:8) {
    p_qtl <- read.table(paste(scenario,"_",rep,"_p_qtl.txt", 
                              sep = ""))
    bla <- apply(X = t(as.numeric(p_qtl[1,])), MARGIN = 2, FUN = funcMaf )
    mafCut <- p_qtl
    p_neu <- read.table(paste(scenario,"_",rep,"_p_neutral.txt", 
                              sep = ""))
    bla <- apply(X = t(as.numeric(p_neu[1,])), MARGIN = 2, FUN = funcMaf )
    mafCutNeu <- p_neu
    for (i in 1:19){
      hom[i,] <- c(rep,scenario,i,F_hom(mafCut,i,1), F_drift(mafCut,i,1), F_hom(mafCutNeu,i,1), F_drift(mafCutNeu,i,1))
    }
    homAll <- rbind(homAll, hom)
    
  }
}
setwd("../ResultsMaf")
write.table(homAll,  paste("homAll",MAF,".txt",sep=""),quote = F, sep = "\t", row.names = F)

homAll
homAll <- as.data.frame(homAll)
homAll$Gen<-as.numeric(homAll$Gen)
homAll$F_hom_qtl <- as.numeric(homAll$F_hom_qtl)
homAll$F_hom_n <- as.numeric(homAll$F_hom_n)
homAll$F_drift_qtl <- as.numeric(homAll$F_drift_qtl)
homAll$F_drift_n <- as.numeric(homAll$F_drift_n)
homAll2 <-  homAll %>% group_by(Gen, Scen) %>% summarise(
  F_hom_qtl = mean(F_hom_qtl),
  F_drift_qtl = mean(F_drift_qtl),
  F_hom_n = mean(F_hom_n),
  F_drift_n = mean(F_drift_n),
)
  
ggplot(homAll2, mapping = aes(Gen, F_drift_qtl, colour = Scen)) +
  geom_line()
ggsave(paste("F_hom_qtl_",MAF,".png",sep=""))

lm1 <- lm(-log(1-F_hom_qtl)~as.factor(Rep)+Gen:Scen, data = homAll)
lm2 <- lm(-log(1-F_drift_qtl)~as.factor(Rep)+Gen:Scen, data = homAll)
lm3 <- lm(-log(1-F_hom_n)~as.factor(Rep)+Gen:Scen, data = homAll)
lm4 <- lm(-log(1-F_drift_n)~as.factor(Rep)+Gen:Scen, data = homAll)

coefficients <- cbind(lm1$coefficients,lm2$coefficients,lm3$coefficients,lm4$coefficients)
colnames(coefficients) <- c("F_hom_qtl", "F_drift_qtl", "F_hom_n", "F_drift_n")
write.table(coefficients,paste("F_Maf",MAF,".txt",sep=""))
setwd("../")
}
