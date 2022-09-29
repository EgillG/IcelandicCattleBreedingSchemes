source("Functions_MoBPS_simulations_20_01_2022.R")
library("tidyverse")
# Generation to cut MAF
GenMaf=6
# Number of generations
Gens = 21
# Number of replicates
NumRep = c(1:10)
#matrices <- c("M2D", "M1D")
matrices = c("M1D_Est", "M1","M2","M2R","M2D","M1D","M1R","M1_05","Ped")

hom <- matrix(nrow=Gens,ncol = 7)
colnames(hom) <- c("Rep","Scen","Gen","F_hom_qtl", "F_drift_qtl","F_hom_n", "F_drift_n")

for (MAF in c(0.001,0.01,0,0.05)) {
homAll = NULL
setwd("Results")
for (scenario in matrices) {
  for (rep in NumRep) {
    p_qtl <- read.table(paste(scenario,"_",rep,"_p_qtl.txt", 
                              sep = ""))
    bla <- apply(X = t(as.numeric(p_qtl[GenMaf,])), MARGIN = 2, FUN = funcMaf )
    mafCut <- p_qtl[,bla>MAF]
    p_neu <- read.table(paste(scenario,"_",rep,"_p_neutral.txt", 
                              sep = ""))
    bla <- apply(X = t(as.numeric(p_neu[GenMaf,])), MARGIN = 2, FUN = funcMaf )
    mafCutNeu <- p_neu[,bla>MAF]
    for (i in 1:Gens){
      hom[i,] <- c(rep,scenario,i,F_hom(mafCut,i,6), F_drift(mafCut,i,6), F_hom(mafCutNeu,i,6), F_drift(mafCutNeu,i,6))
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

setwd("../")
}
