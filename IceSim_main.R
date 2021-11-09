# Updated 11.08.2021 with numerous fixes
# Can use marker-based matrices and pedigree
# Fixed the funcSegQTL function
# Uses bulls equally in each generation.
# Writes relevant parameters (some of them) to a matrix, info.
# SegAlleles currently does not work
# 
library(MoBPS)
library(stringr)
library(miraculix)
library(RandomFieldsUtils)
library(plyr)
library(dplyr)


# library(tidyverse)
# library(plyr)
# Set working drive
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
# args are: no. of repetition, method, number of selected bulls/generation
print(args)
# load QMSim data. 
# 
# Set Minor allele frequency threshold for G2 and G3
G2Maf=0.5
G3Maf=0.1
# Set a frequency for markers to be qtl
QTLMaf <- 0.5
#Number of QTL
qtl = 1000
# Number of offspring each generation
breedSize = 1000
# Number of males and females selected each generation
SelMales=as.integer(args[3])
SelFemales=breedSize/2
#heritability
h2 = 0.4
# Number of markers on each chromosome has to be set:
numMarkers=1000
# Modulus for Neutral marker alleles to track.
NeutralMarkers = 1500
# Variable of maximum number of matings per bull (for OCS).
#MaxMate=2000
MaxMate=breedSize/SelMales
# number of generations for PBLUP
# Pblupgen needs to be 5. Do not change.
Pblupgen = 5
#ssGBLUP
ssGBLUPgen = 4


# Matrix for storing results of genetic gain, heterozygosity and inbreeding.
info = matrix(nrow = (Pblupgen + ssGBLUPgen+1), ncol  = 8)
info=data.frame(info)
colnames(info) = c("BV", "Coancestry","Heteroz", "SegAlleles", "SegQTL", "SegNeutral","F_drift","F_hom")
source("../Functions_MoBPS_simulations_19_10.R")

mapfile = paste("../../QMSim/data_",args[1],".map",sep="")
pedfile = paste("../../QMSim/data_",args[1],".ped",sep="")
freqfile = paste("../../QMSim/plink_",args[1],".frq",sep="")
map <- as.matrix(read.table(mapfile))
ped <- as.matrix(read.table(pedfile))

# read .frq file:
frq <- read.table(freqfile)
# write SNP and chromosome number to a dataframe QTLsnp.
QTLsnp = frq[,c(1,2,3,4,5)]
# Make some filtering on allele frequencies for QTL
# QTLsnp <- frq[frq$V5<QTLMaf & frq$V5>0,c(1,2,3,4,5)]

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

# get.map exports a matrix with map information:
# chromosome, snp-name, bp-position, morgan position
# my map does not have position because QMSim only outputs Morgans, not physical distances.
snpmap = get.map(population, use.snp.nr = TRUE)[,c(1,2)]
# snpnr is the running number of the SNP
snpnr = get.map(population, use.snp.nr = FALSE)[,c(2)]
# bind these two together to have an object with snp nr. and snp chromosome number
snpmap = cbind(snpmap,snpnr)

# I randomly sample markers to assign effects.
# Here I have to sample markers with allele frequency below 10%.
# head(QTLsnp)
# head(snpmap)
# Remove the header line from QTLsnp
QTLsnp = QTLsnp[-1,]
# merge qtlsnps with the snpmap, 
snps = merge(snpmap,QTLsnp, by.x = "snpnr", by.y = "V2")
# remove unnecessary columns
snps = snps[,c(2,3,5,6,7)]
# set column names
colnames(snps) = c("CHR","SNP","A1","A2","MAF")
# observe the snps data frame
head(snps)
tail(snps)
# filter the markers
# Draw alleles with MAF < 0.1 and MAF > 0.
# Make boolean list of IBD markers.
IBDList = rep(FALSE, length(snpmap[,1]))
#IBDList[snpmap[,2] %% NeutralMod == 0] = TRUE
sam = sample(rownames(snps[snps$MAF>0,]),NeutralMarkers)
for (i in sam) {
  IBDList[as.numeric(i)] =TRUE
}
snps=snps[!IBDList,]
snps = snps[snps$MAF<QTLMaf & snps$MAF > 0,]
# randomly draw markers to be QTL, qtl object is the number of QTL
snps = snps[sample(nrow(snps),qtl),]
# dim(snps)
# order the dataframe according to chromosome and SNP within chromosome
snps = snps[order(snps$CHR,snps$SNP),]
# head(snps)

# Have to specify (draw) markers to be effect markers.
# First simulate the distribution, using random gamma
effects = rgamma(qtl, 1, 1)
# plot(effects)
# hist(effects)
# Make matrix for effects with five columns
effx=matrix(nrow = qtl, ncol = 5)
# Try and assign the minor allele always to be negative.
for (bla in 1:qtl) {
  # If A1 is the minor allele, there is a x% probability to assign negative effect
  if (snps[bla,]$A1==1) {
    if (rbinom(1,1,0.8)==1) {
    effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), -effects[bla],0, effects[bla])
    }
    else {
      effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), effects[bla],0, -effects[bla])
    }
    #    print(paste("YES A1:",snps[bla,]$A1, snps[bla,]$MAF, -effects[bla]))
    }
  # else if A2 is the minor allele, there is a 1-x probability to assign negative effect
  else {
    if (rbinom(1,1,0.8)==1) {
      effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), effects[bla],0, -effects[bla])
    }
    else {
      effx[bla,] = c(as.numeric(snps[bla,2]),as.numeric(snps[bla,1]), -effects[bla],0, effects[bla])
    }
   }
}
# head(effx)
# tail(effx)
# dim(effx)
# View(effx)
population = creating.trait(population, real.bv.add = effx)


cohorts = c(length = (Pblupgen+ssGBLUPgen)*2+2)
cohorts[1] = "Cohort_1_M"
cohorts[2] = "Cohort_1_F"
# Set index for cohorts matrix
cohIndex = 3


population <- breeding.diploid(population)

# Add a genotyping array
# Do not use QTL, but assign markers evenly spaced over the genome.
qtlList = get.qtl(population)
markerIncluded = rep(TRUE,population$info$snp[1]*population$info$chromosome)
for (i in qtlList) {
  markerIncluded[i] = FALSE
}
for (i in 1:length(IBDList)) {
  if (IBDList[i] == TRUE) {
  markerIncluded[i] = FALSE
  }
}

# Matrix for storing allele frequency changes 
p_qtl = matrix(nrow = (Pblupgen + ssGBLUPgen+1), ncol  = qtl)
p_marker = matrix(nrow = (Pblupgen + ssGBLUPgen+1), ncol  = sum(markerIncluded))
p_neutral = matrix(nrow = (Pblupgen + ssGBLUPgen+1), ncol  = sum(IBDList))

population <- add.array(population,
          marker.included = markerIncluded,
          array.name = "YarraY")
population <- bv.standardization(population, gen = 1)

for(gen in 1:(Pblupgen+1)){
    print(paste("PBLUP Generation: ", gen))
    population <- breeding.diploid(population, 
                                   heritability = 0.4, 
                                   phenotyping.cohorts = c(cohorts[cohIndex-1]), 
                                   bve = TRUE, # estimate breeding values
                                   relationship.matrix="kinship", # use pedigree relationships
                                   bve.cohorts = cohorts, # Generations of individuals to consider in bve
                                   share.genotyped = 0)

   print(paste("phenotypes for cohort", cohorts[cohIndex-1],"were generated for generation ",gen))
   # Generate a copy of the selected bulls.
   population <- breeding.diploid(population, 
                                  selection.size = c(SelMales,SelFemales),
                                  copy.individual.m = TRUE,
                                  share.genotyped = 0,                                  #add.gen = gen-1,
                                  selection.criteria = "bve",
                                  selection.f.cohorts = cohorts[cohIndex-1],
                                  selection.m.cohorts = cohorts[cohIndex-2],
                                  name.cohort = paste("SelectedBulls",gen, sep=""))
   # Genotype the copies.
   population <- breeding.diploid(population,
                                  genotyped.cohorts = paste("SelectedBulls",gen, sep=""),
                                  genotyped.array = 2
                                  )  
   print(paste(cohorts[cohIndex-2]," were used for generating copies to genotype: ","SelectedBulls",gen,sep=""))
   
   # Apply selection, here I have to set the number of matings so that females
   # have 1 offspring each
   population <- breeding.diploid(population, 
                                 breeding.size = breedSize,
                                 name.cohort = paste("PBLUP",gen+1, sep=""),
                                 selection.size = c(SelMales,SelFemales),
                                 heritability = 0.4,
                                 selection.f.cohorts = cohorts[cohIndex-1],
                                 selection.m.cohorts = paste("SelectedBulls",gen, sep=""),
                                 selection.criteria = "bve",
                                 share.genotyped = 0,
                                 max.mating.pair = 2,
#                                 max.offspring = c(breedSize/SelMales,2),
                                 store.effect.freq = TRUE)

   mrg <- getHet(population,c(cohorts[cohIndex-1], cohorts[cohIndex-2]),numMarkers)
   p_qtl[gen,] <- mrg$freq0
   info[gen,1:6] <- getAll(mrg, gen, c(cohorts[cohIndex-1], cohorts[cohIndex-2]), population, funcSegQTL, func)
   p_marker[gen,] <- get_p(population, c(cohorts[cohIndex-1], cohorts[cohIndex-2]), markerIncluded) 
   p_neutral[gen,] <- get_p(population, c(cohorts[cohIndex-1], cohorts[cohIndex-2]), IBDList) 
   info[gen,7] = sum((p_neutral[gen,] - p_neutral[gen-1,])**2/(p_neutral[1,]*(1-p_neutral[1,])))/NeutralMarkers
   info[gen,8] = 1-sum(2*(p_neutral[gen,]*(1-p_neutral[gen,]))/(2*(p_neutral[1,]*(1-p_neutral[1,]))))/NeutralMarkers 

  # Replace previous male generation with selected bulls by updating the cohort vector
  cohorts[cohIndex-2] =  paste("SelectedBulls",gen, sep="")
  # Assign the males to cohIndex.
  cohorts[cohIndex] =  paste("PBLUP",gen+1,"_M", sep="")
  # Assign females to cohIndex + 1
  cohorts[cohIndex + 1] =  paste("PBLUP",gen+1,"_F", sep="")
  cohIndex =+ cohIndex+2
  print(cohorts)
  print(paste("generation", gen))
            }
summary(population)
cohorts
# Construct two vectors to replace MoBPS generation values with actual generation values
# Repeat the mobps generations, they are Pblup + ssGBLUP + 1 (founder)
# Maybe this can all be replaced by add.gen
MoBPSGen = c(1:(2*(Pblupgen+ssGBLUPgen+1)))
# make a vector of number of times to repat the actual generation number.
# 1 is only once, subsequent Pblup generations are repeated twice,
# ssGBLUP generations only
sbla = seq(2,2,length.out = (length(MoBPSGen)%/%2))
# build the vector using rep.
ve = rep(1:(length(MoBPSGen)%/%2), times = sbla)
gen=7
print("Info matrix looks like this before ssGBLUP stage:\n")
print(info)
print("starting genomic selection")
RDataFile <- paste(args[2],"_",args[1],gen,"BeforeGenomicSelection.RData",sep="")
save.image(file=RDataFile)

# 
###############################################################################
# Start of genomic selection
for (gen in (Pblupgen+2):(ssGBLUPgen+Pblupgen+1)){
print(paste("list of cohorts when genomic selection starts in generation ",gen, ":", toString(cohorts), sep=""))
# Here select 2000 bulls based on parent average to be genotyped.
# Check whether the parent mean is computed correctly (for the selection candidates)
population <- breeding.diploid(population,
                               bve.parent.mean = TRUE,
                               bve.insert.cohorts = cohorts[cohIndex-2]
                               )
# Now use the parent mean to preselect bulls for genotyping
population <- breeding.diploid(population,
                               selection.size = breedSize/2*0.5,
                               copy.individual.m = TRUE,
                               selection.criteria = "bve",
                               selection.m.cohorts = cohorts[cohIndex-2],
#                               add.gen = Pblupgen+gen-2, #check if this is correct.
                               name.cohort = paste("GenoBulls",gen, sep=""))
cohorts[cohIndex-2] <- paste("GenoBulls",gen, sep="")
print(paste("her er eg a undan ebv, kynslod: ",gen))

# Female parents of selection candidates are phenotyped, but not males.
# All females in the current generation are genotyped at this stage
population <- breeding.diploid(population,
                               phenotyping.cohorts = cohorts[cohIndex-3], #phenotype previous generation
                               genotyped.cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]), #genotype all females, and 2000 bulls
                               bve = FALSE)
print(paste("her er eg aftur eftir ebv, kynslod: ",gen))

# Now use DMU for single step breeding value estimation.
## First write necessary files.
ped <- data.frame(get.pedigree(population, gen = 1:length(population$breeding)))
ped2 = data.frame(get.pedigree(population, gen = 1:length(population$breeding), id = T))
matc=data.frame(ped$offspring, ped2$offspring)
colnames(matc) = c("ID","RecodeID")

pedig <- data.frame(get.pedigree(population, gen=1:length(population$breeding)))
# Put in generation number (not exactly correct but good enough):
pedig$Generation <- substr(str_extract(pedig$offspring,"_.$"),2,2)
pedig[is.na(pedig$Generation),]$Generation <- substr(str_extract(pedig[is.na(pedig$Generation),]$offspring,"_..$"),2,3)
pedig$Generation <- mapvalues(pedig$Generation, MoBPSGen, ve, F)

pedig$offspring <- mapvalues(pedig$offspring, matc$ID, matc$RecodeID, warn_missing = F)
pedig$father <- mapvalues(pedig$father, matc$ID, matc$RecodeID, warn_missing = F)
pedig$mother <- mapvalues(pedig$mother, matc$ID, matc$RecodeID, warn_missing = F)
pedig <- unique(pedig)
pedig <- pedig[pedig$father!=pedig$offspring,]
write.table(pedig, 
            "pedigree.txt",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE, 
            sep = " ")

pheno = data.frame(t(get.pheno(population, cohorts = cohorts)))
pheno$ID <- rownames(pheno)
pheno$ID <- mapvalues(pheno$ID, matc$ID, matc$RecodeID, F)
pheno$mean = 1
pheno[is.na(pheno)] <- -99
pheno$Trait.1 <- round(pheno$Trait.1,2)
pheno = pheno[,c(2,3,1)]
write.table(pheno, "phenotypes.txt",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE, 
            sep = " ")
genos = data.frame(t(get.geno(population, 
                   cohorts = cohorts[-seq(2, to = Pblupgen*2+2, by = 2)],
                   non.genotyped.as.missing = TRUE)))
genos[is.na(genos)] <- 9
genos$ID <- row.names(genos)
genos$RecodeID <- mapvalues(genos$ID, matc$ID, matc$RecodeID, F)
genos <- genos %>% relocate(RecodeID,ID, .before = X1)
write.table(genos,
            "Gmatrix/gmat.dat", 
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE, 
            sep = " ")
    # make mapfile for GMATRIX
    # mapfile is simply the number of snps and a 0 or 1, depending on whether
    # the marker is genotyped or not
      mapfile = cbind(seq(1,population$info$snp[1]*population$info$chromosome),
                      rep(0,population$info$snp[1]*population$info$chromosome))
      # This loop checks whether the marker is genotyped
      # (I know that there is a more efficient way)
      for (i in 1:dim(mapfile)[1]) {
        if (markerIncluded[i]) {
          mapfile[i,2] <- 1
        }
      }
    write.table(cbind(mapfile,rep(1,dim(mapfile)[1])), 
                "Gmatrix/gmat.map",
                row.names = FALSE, 
                col.names = FALSE,
                quote = F,
                sep = " ")
    # Make id file for GMATRIX
    idfile <- cbind(genos$RecodeID,genos$RecodeID,rep(1,length(genos$RecodeID)))
    # Should this code not be run? I don't remember what it does
    # genotyped = get.genotyped(population, gen=1:length(population$breeding))
    # for (i in 1:dim(idfile)[1]) {
    #   if (genotyped[i]) {
    #     idfile[i,3] <- 1
    #   }
    # }
    write.table(idfile, 
                "Gmatrix/gmat.id",
                row.names = FALSE, 
                col.names = FALSE,
                quote = F,
                sep = " ")
# write the necessary variances (I use the true values)
varA = sum(get.qtl.variance(population = population, cohorts=cohorts)[[1]][,3])
varE = ((1-h2)/h2)*varA
write.table(rbind(c(1,1,1, varA),c(2,1,1,varE)),
                "variances.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = F,
                sep = " ")

system(paste("python3 gmatPar.py ", getwd(), " M1", sep = ""))
system(paste("bash DMU.sh ", args[2], sep=""))

#DMU ebvs are inserted back into MoBPS and extracted again for ocs computations (silly indeed).
dmuBVE = data.frame(read.table("dmuSS.SOL", skip = 2, stringsAsFactors = F)[,c(5,8)])
dmuBVE$V5 <- mapvalues(dmuBVE$V5, matc$RecodeID, matc$ID, F)
population <- insert.bve(population, dmuBVE)
###############################################################################
# Try OCS
###############################################################################
  ped <- data.frame(get.pedigree(population, gen = 1:length(population$breeding)))
  ped2 = data.frame(get.pedigree(population, gen = 1:length(population$breeding), id = T))
  matc=data.frame(ped$offspring, ped2$offspring)
  colnames(matc) = c("ID","RecodeID")
  ped <- data.frame(get.pedigree(population, cohorts = cohorts[c(cohIndex-1,cohIndex-2)]))
    #ped <- data.frame(get.pedigree(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2])))
  # Assign sex
  # First assign all animals as 1 (males)
  ped$sex = "1"
  # Assign 2 to animals with F as their first letter
  ped[substr(ped$offspring,1,1) == "F",]$sex = "2"
  bve <- get.bve(population, cohorts = cohorts[c(cohIndex-1,cohIndex-2)])
  # The input to EVA, evaIn contains ped and ebvs.
  # Use merge to combine ped and bve.
  evaIn <- merge(ped, data.frame(t(bve)), by.x = "offspring", by.y = "row.names", sort = F, all.x = TRUE)
  evaIn$ID = mapvalues(evaIn$offspring,matc$ID,matc$RecodeID, warn_missing = FALSE)
  
  #make dataframe of individual MoBPS IDs and IDs for EVA.
  evaIn$mother <- 0
  evaIn$father <- 0
  evaIn$generation = 1
  evaIn$maxmatings = 1
  evaIn[evaIn$sex==2,]$maxmatings = args[3]
  evaIn = evaIn[c(6,2,3,4,7,8,5,1)]
  colnames(evaIn)=c("ID", "sire", "dam", "sex" ,"generation", "maxmatings",
                    "ebv","MoBPSID")
  evaIn$ebv = round(evaIn$ebv, 2)
  meanEBV=mean(evaIn[evaIn$sex== 2 & evaIn$maxmatings>0,]$ebv)
  # Write selection candidates before adding pseudofemale
  write.table(evaIn[,c(1,4)], "SelCands",
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  evaIn=rbind(evaIn[evaIn$sex==1,],c("00001","0","0",2,1,args[3],round(meanEBV,2),"PseudoFemale"))
  write.table(evaIn, "evaIn.txt",
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  if (args[2]=="Ped") {
  bla = kinship.exp.store(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]), depth.pedigree = Pblupgen+ssGBLUPgen)
  newFile=matrix(0,ncol=3, nrow=dim(bla)[1]*(dim(bla)[1]-1)/2+dim(bla)[1])
  count=1
  for (i in 1:dim(bla)[1]) {
    for (j in i:dim(bla)[1]) {
      #print(paste(i,j))
      newFile[count,1] = rownames(bla)[i]
      newFile[count,2] = colnames(bla)[j]
      newFile[count,3] = bla[i,j]
      count = count + 1
    }
  }
  
  write.table(newFile,"ReducedMatrix", quote = F, row.names = F, col.names = F)
  }
# Write genotypes if marker-based GMATRIX is used
  if (substr(args[2],1,1) == "M") { 
    # Here I retrieve the genotypes of desired cohorts
print(paste("these cohorts genotypes were written in generation ",gen, " ", toString(cohorts[-seq(2, to = Pblupgen*2+2, by = 2)]), sep=""))

    genos = t(get.geno(population, 
                       cohorts = cohorts[-seq(2, to = Pblupgen*2+2, by = 2)],
                       non.genotyped.as.missing = TRUE))
    genos = data.frame(genos)
    genos$ID <- row.names(genos)
    # Assign correct recoded IDs
    genos$RecodeID <- mapvalues(genos$ID,
                                matc$ID,
                                matc$RecodeID,
                                warn_missing = FALSE)
    # Move ID columns to the front.
    genos <- genos %>% relocate(RecodeID,ID, .before = X1)
    # Set missing genotypes to 9 (these are effect markers)
    genos[is.na(genos)] <- 9
    write.table(genos,"gmat.dat", 
                row.names = FALSE, 
                col.names = FALSE,
                quote = FALSE, 
                sep = " ")
    # make mapfile for GMATRIX
    # mapfile is simply the number of snps and a 0 or 1, depending on whether
    # the marker is genotyped or not
      mapfile = cbind(seq(1,population$info$snp[1]*population$info$chromosome),
                      rep(0,population$info$snp[1]*population$info$chromosome))
      # This loop checks whether the marker is genotyped
      # (I know that there is a more efficient way)
      for (i in 1:dim(mapfile)[1]) {
        if (markerIncluded[i]) {
          mapfile[i,2] <- 1
        }
      }
    write.table(cbind(mapfile,rep(1,dim(mapfile)[1])), 
                "gmat.map",
                row.names = FALSE, 
                col.names = FALSE,
                quote = F,
                sep = " ")
    # Make id file for GMATRIX
    idfile <- cbind(genos$RecodeID,genos$RecodeID,rep(1,length(genos$RecodeID)))

    write.table(idfile, 
                "gmat.id",
                row.names = FALSE, 
                col.names = FALSE,
                quote = F,
                sep = " ")
    if (args[2] == "M1D" || args[2] == "M2D") {
      # find the frequency of alleles in the base generation.
      freqfile <-  get.geno(population, gen = 1, non.genotyped.as.missing = FALSE)
      # rowMeans can be used to compute allele frequencies (for the larger allele)
      freq <- rowMeans(freqfile)/2
      write.table(freq, 
                  "mapbase.dat",
                  row.names = TRUE, 
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }
    if (args[2] == "M1D_Est" || args[2] == "M2D_Est") {
      # Now use the genotypes of older bulls. This uses the first three
      # generations, should perhaps be changed.
      freqfile <-  get.geno(population, 
                            cohorts = c(cohorts[c(1,3,5)]),
                            non.genotyped.as.missing = FALSE)
      # rowMeans can be used to compute allele frequencies (for the larger allele)
      freq <- rowMeans(freqfile)/2
      write.table(freq, 
                  "mapbase.dat",
                  row.names = TRUE, 
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }
    
    if (args[2] == "M1R" || args[2] == "M2R") {
      freqfile <-  get.geno(population, 
                            cohorts = c(cohorts[cohIndex-1],cohorts[cohIndex-2]),
                            non.genotyped.as.missing = FALSE)
      # rowMeans can be used to compute allele frequencies (for the larger allele)
      freq <- rowMeans(freqfile)/2
      write.table(freq, 
                  "mapbase.dat",
                  row.names = TRUE, 
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }    
  }

# Write the haplotypes of genotyped animals if method is haplotype-based.
# Here I have to take care of cohorts.
if (substr(args[2],1,1) == "H") { 
haplo=get.haplo(population, cohorts = cohorts[-c(2,4,6,8,10)],
                non.genotyped.as.missing = TRUE)
haplo[is.na(haplo)] <- 9
print("Writing haplotypes")
write.table(haplo,"haplo.hap", row.names = FALSE, col.names = FALSE,
             quote = FALSE, sep = "\t")
# Use stringr to replace _set1 and _set2 with ""
id=colnames(haplo)
id=str_replace(id,"_set[1-2]","")
id <- mapvalues(id,matc$ID,matc$RecodeID, F)
# Write sample IDs
write.table(unique(id),"haplo.sample", row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")
# Make and write genetic map
ma=get.map(population)
ma=data.frame(ma[,1],ma[,2], population$info$snp.position,"A","C")
write.table(ma, "haplo.map", row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")
rm(ma)
rm(haplo)
}
system(paste("bash EVA.sh",args[1],args[2],args[3],sep = " "))

EVAfilename = paste("evaSim/Candidates.txt",sep="")
# Read the Candidates.txt file from EVA
df <- read.table(EVAfilename, skip = 6, header = T, nrows = breedSize)

# Select males
df<- df[df$Sex == 1,]
# sort the list descending according to number of matings allocated and select males
# according to SelMales
df <- head(df[order(-df$N.matings),], SelMales)
print(paste("her er eg fyrir OCS urval, kynslod:", gen))

# Assign huge breeding value estimate to selected bulls.
dfBVE = cbind(df$txt, 1e9)
population <- insert.bve(population, bves = dfBVE)

##############################################################
# Make R Data file. This is temporary
#RDataFile <- paste(args[2],"/gen",gen,args[2],"_",args[1],"TempStuff.RData",sep="")
#save.image(file=RDataFile)
##############################################################

population <- breeding.diploid(population, 
                               breeding.size = breedSize,
                               name.cohort = paste("ssGBLUP",gen, sep=""),
                               selection.size = c(SelMales,SelFemales),
                               heritability = 0.4,
                               selection.f.cohorts = cohorts[cohIndex-1],
                               selection.m.cohorts = cohorts[cohIndex-2],
                               selection.criteria = "bve",
                               max.mating.pair = 2,
                               max.offspring = c(breedSize/SelMales,2))
# Here I can decide whether to make new cohort of selected bulls,
# to replace the cohort of genotyped bulls.
a <- get.pedigree(population, gen = length(population$breeding), raw=TRUE)
table(a[,6])	
   print(paste("write genomic information for generation",gen, sep = " "))
   mrg <- getHet(population,c(cohorts[cohIndex-1], cohorts[cohIndex-2]),numMarkers)
   p_qtl[gen,] <- mrg$freq0
   info[gen,1:6] <- getAll(mrg, gen, c(cohorts[cohIndex-1], cohorts[cohIndex-2]), population, funcSegQTL, func)
   p_marker[gen,] <- get_p(population, c(cohorts[cohIndex-1], cohorts[cohIndex-2]), markerIncluded)
   p_neutral[gen,] <- get_p(population, c(cohorts[cohIndex-1], cohorts[cohIndex-2]), IBDList)
   info[gen,7] = sum((p_neutral[gen,] - p_neutral[gen-1,])**2/(p_neutral[1,]*(1-p_neutral[1,])))/NeutralMarkers
   info[gen,8] = 1-sum(2*(p_neutral[gen,]*(1-p_neutral[gen,]))/(2*(p_neutral[1,]*(1-p_neutral[1,]))))/NeutralMarkers


# Update cohortlist
cohorts[cohIndex] =  paste("ssGBLUP",gen,"_M", sep="")
cohorts[cohIndex+1] =  paste("ssGBLUP",gen,"_F", sep="")
cohIndex =+ cohIndex+2
print(paste("end of generation",gen, sep = " "))      
RDataFile <- paste(args[2],"_",args[1],gen,".RData",sep="")
save.image(file=RDataFile)
}
summary(population)
RDataFile <- paste(args[2],"_",args[1],".RData",sep="")
save.image(file=RDataFile)


qt <- get.qtl.variance(population, gen = 1:length(population$breeding))
qt=data.frame(qt)
tail(qt)
head(qt)
write.table(qt, paste(args[2],"_",args[1],"_qt.txt",sep=""), quote = F, sep = "\t")
write.table(info, file = paste(args[2],"_",args[1],"_info.txt",sep=""))

write.table(p_qtl, file = paste(args[2],"_",args[1],"_p_qtl.txt",sep=""))
write.table(p_marker, file = paste(args[2],"_",args[1],"_p_markers.txt",sep=""))
write.table(p_neutral, file = paste(args[2],"_",args[1],"_p_neutral.txt",sep=""))
