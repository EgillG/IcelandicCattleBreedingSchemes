source("Functions_MoBPS_simulations_20_01_2022.R")
library("tidyverse")
# Generation to cut MAF
GenMaf=6
# Number of generations
Gens = 21
# Number of replicates
NumRep = 1:10
matrices = c("M1D_Est", "M1","M2","M2R","M2D","M1D","M1R","M1_05","Ped")
hom <- matrix(nrow=Gens,ncol = 5)
colnames(hom) <- c("Rep","Scen","Gen","F_hom_m", "F_drift_m")

for (MAF in c(0.001,0.01,0.05)) {
#for (MAF in c(0.05)) {
homAll = NULL
setwd("Results")
for (scenario in matrices) {
  for (rep in NumRep) {
    p_qtl <- read.table(paste(scenario,"_",rep,"_p_markers.txt",sep = ""))
    bla <- apply(X = t(as.numeric(p_qtl[GenMaf,])), MARGIN = 2, FUN = funcMaf )
    mafCut <- p_qtl[,bla>MAF]
    for (i in 1:Gens){
      hom[i,] <- c(rep,scenario,i,F_hom(mafCut,i,6), F_drift(mafCut,i,6))
    }
    homAll <- rbind(homAll, hom)
    
  }
}
setwd("../ResultsMaf")
write.table(homAll,  paste("homMarkers",MAF,".txt",sep=""),quote = F, sep = "\t", row.names = F)

homAll
homAll <- as.data.frame(homAll)
homAll$Gen<-as.numeric(homAll$Gen)
homAll$F_drift_qtl <- as.numeric(homAll$F_drift_m)
homAll$F_drift_n <- as.numeric(homAll$F_drift_m)
homAll2 <-  homAll %>% group_by(Gen, Scen) %>% summarise(
  F_hom_m = mean(F_hom_m),
  F_drift_m = mean(F_drift_m),
)
  
ggplot(homAll2, mapping = aes(Gen, F_drift_m, colour = Scen)) +
  geom_line()
ggsave(paste("F_hom_m_",MAF,".png",sep=""))

setwd("../")
}
