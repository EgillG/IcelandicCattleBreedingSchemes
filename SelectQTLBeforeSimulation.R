#  
library(MoBPS)
library(stringr)
library(miraculix)
library(RandomFieldsUtils)
library(plyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
# args are: no. of QTL, neutral markers,replicate number
print(args)

# load QMSim data. 
#Number of QTL
qtl = as.integer(args[1])
NeutralMarkers = as.integer(args[2])
NoRep = args[3]

#setwd("C:/Users/au589863/OneDrive - Aarhus universitet/Documents/phdProject4_Simulations/MoBPS/")
mapfile = paste("data_",args[3],".map",sep="")
pedfile = paste("data_",args[3],".ped",sep="")
freqfile = paste("plink_",args[3],".frq",sep="")
map <- as.matrix(read.table(mapfile))
ped <- as.matrix(read.table(pedfile))

# read frq file
frq <- read.table(freqfile)
# write SNP and chromosome number to a dataframe QTLsnp.
QTLsnp = frq[,c(1,2,3,4,5)]
# Make some filtering on allele frequencies for QTL

# Seems to convert centimorgans into morgans.
map[,3]=as.numeric(map[,3])*0.01

# Convert PED-file into haplotype dataset (one haplotype per colum)
# This is according to page 26 in the MoBPS manual.
nsnp <- (ncol(ped)-6)/2
haplo1 <- ped[,1:nsnp*2+6-1]
haplo2 <- ped[,1:nsnp*2+6]
haplo <- t(rbind(haplo1, haplo2)[c(0,nrow(haplo1)) + sort(rep(1:nrow(haplo1),2)),])
# I set dataset = haplo-1. It seems that dataset=haplo gives wrong results.
# I think that bpcm.conversion does not make a difference.
population <- creating.diploid(dataset = haplo-1,
                               map = map,
                               bpcm.conversion = 0,
                               share.genotyped = 0)
summary(population)
rm(haplo1,haplo2,ped,haplo,map)

# get.map exports a matrix with map information:
# chromosome, snp-name, bp-position, morgan position
# my map does not have position because QMSim only outputs Morgans, not physical distances.
snpmap = get.map(population, use.snp.nr = TRUE)[,c(1,2)]
# snpnr is the running number of the SNP
snpnr = get.map(population, use.snp.nr = FALSE)[,c(2)]
# bind these two together to have an object with snp nr. and snp chromosome number
snpmap = cbind(snpmap,snpnr)

# I randomly sample markers to assign effects.
# head(QTLsnp)
# head(snpmap)
# Remove the header line from QTLsnp
# merge qtlsnps with the snpmap, 
snps = merge(snpmap,QTLsnp, by.x = "snpnr", by.y = "V1")
rm(QTLsnp)
# remove unnecessary columns
snps = snps[,c(2,3,5,6,7)]
# set column names
colnames(snps) = c("CHR","SNP","A1","A2","MAF")
# observe the snps data frame
head(snps)
tail(snps)
# filter the markers
# Draw alleles with MAF < 0.1 and MAF > 0.
# Make boolean list of neutral markers.
NeuList = rep(FALSE, length(snpmap[,1]))
rm(snpmap)

sam = sample(rownames(snps[snps$MAF>0,]),NeutralMarkers)
for (i in sam) {
  NeuList[as.numeric(i)] =TRUE
}
snps=snps[!NeuList,]
snps = snps[snps$MAF<0.5 & snps$MAF > 0,]
# randomly draw markers to be QTL, qtl object is the number of QTL
snps = snps[sample(nrow(snps),qtl),]
# dim(snps)
# order the dataframe according to chromosome and SNP within chromosome
snps = snps[order(snps$CHR,snps$SNP),]
# head(snps)

# Have to specify (draw) markers to be effect markers.
# First simulate the distribution, using random gamma
effects = rgamma(qtl, shape = 0.4, scale = 1.66)
# plot(effects)
# hist(effects)
# Make matrix for effects with five columns
effx=matrix(nrow = qtl, ncol = 5)
# The probability in rbinom has to be changed from 0.5 if the distribution of allele effects if effects are supposed to be assigned according to frequency distribution.

for (bla in 1:qtl) {
  # If A1 is the minor allele, there is a x% probability to assign negative effect
  if (snps[bla,]$A1==1) {
    if (rbinom(1,1,0.5)==1) {
    effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), -effects[bla],0, effects[bla])
    }
    else {
      effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), effects[bla],0, -effects[bla])
    }
    #    print(paste("YES A1:",snps[bla,]$A1, snps[bla,]$MAF, -effects[bla]))
    }
  # else if A2 is the minor allele, there is a 1-x probability to assign negative effect
  else {
    if (rbinom(1,1,0.5)==1) {
      effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), effects[bla],0, -effects[bla])
    }
    else {
      effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), -effects[bla],0, effects[bla])
    }
   }
}
population = creating.trait(population, real.bv.add = effx)
rm(effx)

population <- breeding.diploid(population)

# Add a genotyping array
# Do not use QTL, but assign markers evenly spaced over the genome.
qtlList = as.integer(rownames(snps))
markerIncluded = rep(TRUE,length(NeuList))
for (i in qtlList) {
  markerIncluded[i] = FALSE
}
for (i in 1:length(NeuList)) {
  if (NeuList[i] == TRUE) {
  markerIncluded[i] = FALSE
  }
}

# Filter out some of the very low MAF variants.
colnames(frq) <- c("SNP","CHR","A1","A2","MAF","NoObs")
for (i in 1:length(markerIncluded)) {
  if  (markerIncluded[i]==FALSE) {
    next
  } else if (frq[i,]$MAF==0.0) {
     markerIncluded[i] = FALSE
  } else if (frq[i,]$MAF<0.01) {
    if (rbinom(1,1,0.5)==1){
      markerIncluded[i] = FALSE
    } 
#  } else if (frq[i,]$MAF>=0.01 & frq[i,]$MAF<0.02) {
#    if (rbinom(1,1,0.2)==1){
#      markerIncluded[i] = FALSE
#    }
  }
} 

# Print some marker statistics:
sum(markerIncluded)
png(paste("Figures/markerHistogram_",args[3],".png",sep=""))
hist(frq$MAF[markerIncluded], nclass=50)
dev.off()
png(paste("Figures/markerHistogram_",args[3],"_01.png",sep=""))
hist(frq[frq$MAF>0.01 & markerIncluded,]$MAF)
dev.off()
png(paste("Figures/markerHistogram_",args[3],"_001.png",sep=""))
hist(frq[frq$MAF>0.001 & markerIncluded,]$MAF)
dev.off()
print(paste("Number of markers above MAF = 0.01: ", length(frq[frq$MAF>0.01 & markerIncluded,]$MAF)))
print(paste("Number of markers below MAF = 0.01: ", length(frq[frq$MAF<0.01 & markerIncluded,]$MAF)))
rm(frq)

population <- add.array(population,
          marker.included = markerIncluded,
          array.name = "YarraY")
print("Saving R image.")
RDataFile <- paste("data_",args[3],".RData",sep="")
rm(args)
save.image(file=RDataFile)

