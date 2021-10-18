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
breedSize = 10000
# Number of males and females selected each generation
SelMales=as.integer(args[3])
SelFemales=breedSize/2
#heritability
h2 = 0.4
# Number of markers on each chromosome has to be set:
numMarkers=1000
# Number of Neutral marker alleles to track. Not currently in use.
NeutralMarkers = 500
# Variable of maximum number of matings per bull (for OCS).
#MaxMate=2000
MaxMate=breedSize/SelMales
# number of generations for PBLUP
# Pblupgen needs to be 5. Do not change.
Pblupgen = 5
#ssGBLUP
ssGBLUPgen = 15

# Matrix for storing results of genetic gain, heterozygosity and inbreeding.
info = matrix(nrow = (Pblupgen + ssGBLUPgen+1), ncol  = 6)
info=data.frame(info)
colnames(info) = c("BV", "Coancestry","Heteroz", "SegAlleles", "SegQTL", "DriftVar")
# Matrix for storing allele frequency changes for computing drift variance
# (variance of allele frequency changes)
p = matrix(nrow = (Pblupgen + ssGBLUPgen+1), ncol  = qtl)

#####################
#define functions
#####################
# function to check whether alleles are segregating or fixed.
# Might be better to return 0 or 1 depending on fixation, to reduce memory use.
func <- function(x) {
  if (length(table(x)) == 1){
    # if length of table is 1, one allele or the other is fixed.
    # 0 is assigned to a fixed allele, 1 to a segregating allele.
    return(0)
  }
  else {
    return(1)
  }
}
# This function does the same as before, but on different data, so it counts whether
# number of animals homozygous for one QTL allele is equal to the number of loci
funcSegQTL <- function(x){
  # x is the input frequencies (homo1,hetero,homo2)
  # If homo1 and hetero is 0, return 0, or if 
  # hetero and homo2 is 0, return 0
  # which means the QTL is fixed.
  if (x[1] == 0 & x[2] == 0) {
    return(0)
  }
  else if (x[2] == 0 & x[3] == 0) {
    return(0)
  }
  else {
    return(1)
  }
}


# return minor allele frequency from vector of frequencies.
funcMaf <- function(x){
  # x is the input
  #   If x is larger than 0.5, return 1-x
  if (x > 0.5) {
    return(1-x)
  }  else {
    return(x)
  }

}
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

# Set values for cohort one into info matrix
# Get kinship from IBD. 
info[1,2] <- kinship.emp.fast(population = population, gen = 1,ibd.obs = breedSize/10, hbd.obs = breedSize/10)[1]
population <- breeding.diploid(population, 
                               heritability = 0.4, 
                               phenotyping = "non_obs_f",
                               # phenotyping = "non_obs_f", # phenotype non-phenotyped females
                               bve = TRUE, # estimate breeding values
                               relationship.matrix="kinship", # use pedigree relationships
                               bve.cohorts = cohorts, # Generations of individuals to consider in bve
                               share.genotyped = 0,
                               store.effect.freq = TRUE)
bv <- get.bv(population, gen = 1)
info[1,1] <- mean(bv)
df=data.frame(get.effect.freq(population, gen = 1))
# info$real.bv.add : Lists with an overview of all single marker QTLs for each trait
head(population$info$real.bv.add[[1]][,c(1,3,4,5)])
# take out the effects
effects = data.frame(population$info$real.bv.add[[1]][,c(1,2,3,4,5)])
head(effects)
# change to character
effects$X1=effects$X1 + (effects$X2-1)*numMarkers
effects$X1=as.character(effects$X1)
effects = effects[,-2]
# find the frequencies.
head(df)
df$X1 = row.names(df)
# merge 
mrg = merge(effects, df, by = "X1")
mrg$freq1 = (mrg$Homo0*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
mrg$freq2 = (mrg$Homo1*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
colnames(mrg) = c("SNP","Effect0", "EffectHet","Effect1", "Homo0","Hetero","Homo1","freq0","freq1")
p[1,] <- mrg$freq0
# pring mean heterozygosity and store in info matrix.
print(c("mean heterozygosity of QTL, gen:", paste(1), paste(mean(sum(mrg$Hetero)/(sum(mrg$Homo0+mrg$Hetero+mrg$Homo1))))))
info[1, 3] = mean(sum(mrg$Hetero)/(sum(mrg$Homo0+mrg$Hetero+mrg$Homo1)))
freq1=data.frame(mrg[,c(2,8)])
freq2=data.frame(mrg[,c(4,9)])
afi <- apply(mrg[,c(5,6,7)], MARGIN = 1, FUN = funcSegQTL)
sum(afi)
info[1, 5] <- sum(afi)/qtl
# compute total heterozygosity of markers and number of fixed alleles
# Get genotypes: 
bl=get.geno(population, gen = 1)
# Assign
fixSeg <- apply(bl, MARGIN = 1, FUN = func)
# Sum of fixSeg is the number of segregating alleles
sum(fixSeg)
# Put in percentage of alleles segregating.
info[1,4] <- sum(fixSeg)/length(fixSeg)
colnames(freq1) = c("effect","frequency")
colnames(freq2) = c("effect","frequency")
dat = rbind(freq1,freq2)
#plot(dat$frequency,dat$effect)
gen=2

# Add a genotyping array
# Do not use QTL, but assign markers evenly spaced over the genome.
qtlList = get.qtl(population)

markerIncluded = rep(TRUE,population$info$snp[1]*population$info$chromosome)
for (i in qtlList) {
  markerIncluded[i] = FALSE
}
population <- add.array(population,
          marker.included = markerIncluded,
          array.name = "YarraY")


for(gen in 1:Pblupgen+1){
    print(paste("PBLUP Generation: ", gen))
    population <- breeding.diploid(population, 
                                   heritability = 0.4, 
                                   phenotyping.cohorts = c(cohorts[cohIndex-3]), 
                                   bve = TRUE, # estimate breeding values
                                   relationship.matrix="kinship", # use pedigree relationships
                                   bve.cohorts = cohorts, # Generations of individuals to consider in bve
                                   share.genotyped = 0)

   print(paste("phenotypes for cohort", cohorts[cohIndex-3],"were generated for generation ",gen))
   # Generate a copy of the selected bulls.
   population <- breeding.diploid(population, 
                                  selection.size = c(SelMales,SelFemales),
                                  copy.individual.m = TRUE,
                                  share.genotyped = 0,                                  #add.gen = gen-1,
                                  selection.criteria = "bve",
                                  selection.f.cohorts = cohorts[cohIndex-1],
                                  selection.m.cohorts = cohorts[cohIndex-2],
                                  name.cohort = paste("SelectedBulls",gen-1, sep=""))
   # Genotype the copies.
   population <- breeding.diploid(population,
                                  genotyped.cohorts = paste("SelectedBulls",gen-1, sep=""),
                                  genotyped.array = 2
                                  )  
   print(paste(cohorts[cohIndex-2]," were used for generating copies to genotype: ","SelectedBulls",gen-1,sep=""))
   
   # Apply selection, here I have to set the number of matings so that females
   # have 1 offspring each
   population <- breeding.diploid(population, 
                                 breeding.size = breedSize,
                                 name.cohort = paste("PBLUP",gen, sep=""),
                                 selection.size = c(SelMales,SelFemales),
                                 heritability = 0.4,
                                 selection.f.cohorts = cohorts[cohIndex-1],
                                 selection.m.cohorts = paste("SelectedBulls",gen-1, sep=""),
                                 selection.criteria = "bve",
                                 share.genotyped = 0,
                                 max.mating.pair = 2,
#                                 max.offspring = c(breedSize/SelMales,2),
                                 store.effect.freq = TRUE)

        # To compute the heterozygosity
        df=data.frame(get.effect.freq(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2])))
        # info$real.bv.add : Lists with an overview of all single marker QTLs for each trait
        head(population$info$real.bv.add[[1]][,c(1,3,4,5)])
        # take out the effects
        effects = data.frame(population$info$real.bv.add[[1]][,c(1,2,3,4,5)])
        head(effects)
        # change to character
        effects$X1=effects$X1 + (effects$X2-1)*numMarkers
        effects$X1=as.character(effects$X1)
        effects = effects[,-2]
        # find the frequencies.
        head(df)
        df$X1 = row.names(df)
        # merge 
        mrg = merge(effects, df, by = "X1")
        mrg$freq1 = (mrg$Homo0*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
        mrg$freq2 = (mrg$Homo1*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
        colnames(mrg) = c("SNP","Effect0", "EffectHet","Effect1", "Homo0","Hetero","Homo1","freq0","freq1")
        p[gen,] <- mrg$freq0
        info[gen,6] <- var(p[gen,]-p[gen-1,])
        # print mean heterozygosity and store in info matrix.
        print(c("mean heterozygosity of QTL, gen:", paste(gen), paste(mean(sum(mrg$Hetero)/(sum(mrg$Homo0+mrg$Hetero+mrg$Homo1))))))
        info[gen, 3] = mean(sum(mrg$Hetero)/(sum(mrg$Homo0+mrg$Hetero+mrg$Homo1)))
        freq1=data.frame(mrg[,c(2,8)])
        freq2=data.frame(mrg[,c(4,9)])
        afi <- apply(mrg[,c(5,6,7)], MARGIN = 1, FUN = funcSegQTL)
        sum(afi)
        info[gen, 5] <- sum(afi)/qtl
        # compute total heterozygosity of markers and number of fixed alleles
        # Get genotypes: 
        bl=get.geno(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]))
        # Assign
        fixSeg <- apply(bl, MARGIN = 1, FUN = func)
        # Sum of fixSeg is the number of segregating alleles
        sum(fixSeg)
        # Put in percentage of alleles segregating.
        info[gen,4] <- sum(fixSeg)/length(fixSeg)
        colnames(freq1) = c("effect","frequency")
        colnames(freq2) = c("effect","frequency")
        dat = rbind(freq1,freq2)
        #plot(dat$frequency,dat$effect)
        # hist(dat$frequency, breaks = 100)
        # Get kinship from IBD. 
        info[gen,2] <- kinship.emp.fast(population = population,
                                              cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]),
                                              ibd.obs = breedSize/10,
                                              hbd.obs = breedSize/10)[1]
        bv <- get.bv(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]))
        info[gen,1] <- mean(bv)
  # Replace previous male generation with selected bulls by updating the cohort vector
  cohorts[cohIndex-2] =  paste("SelectedBulls",gen-1, sep="")
  # Assign the males to cohIndex.
  cohorts[cohIndex] =  paste("PBLUP",gen,"_M", sep="")
  # Assign females to cohIndex + 1
  cohorts[cohIndex + 1] =  paste("PBLUP",gen,"_F", sep="")
  cohIndex =+ cohIndex+2
  print(paste(cohorts,"generation", gen))
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
# 
###############################################################################
# Start of genomic selection
for (gen in (Pblupgen+2):(ssGBLUPgen+Pblupgen+1)){
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
                               name.cohort = paste("GenoBulls",gen-1, sep=""))
cohorts[cohIndex-2] <- paste("GenoBulls",gen-1, sep="")
print(paste("her er eg a undan ebv, kynslod: ",gen))
# Female parents of selection candidates are phenotyped, but not males.
# All females in the current generation are genotyped at this stage
population <- breeding.diploid(population,
                               phenotyping.cohorts = cohorts[cohIndex-3], #phenotype previous generation
                               genotyped.cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]), #genotype all females, and 2000 bulls
                               heritability = 0.4,
                               bve = TRUE,
                               relationship.matrix="vanRaden",
                               remove.effect.position = TRUE,
                               singlestep.active = TRUE,
                               bve.cohorts = cohorts,
                               genotyped.array = 2)
#                                 bve.gen=1:(index*2+1))# Has to be 10 or 11 because the copied individuals
#                                 above (for genotyping of bulls) get their own generation.
print(paste("her er eg aftur eftir ebv, kynslod: ",gen))

###############################################################################
# Try OCS
###############################################################################
# Is the pedigree containing all the relevant animals? 
# Maybe produce whole pedigree rather than specify cohorts.
# Start here to construct the pedigree for input to EVA.
# First take out the whole pedigree, to make lists of IDs to convert.
ped <- data.frame(get.pedigree(population, gen = 1:length(population$breeding)))

matc=data.frame(ped$offspring, row.names(ped))
colnames(matc) = c("ID","RecodeID")
ped <- get.pedigree(population, cohorts = cohorts)
ped <- data.frame(ped)
# Assign sex
# First assign all animals as 1 (males)
ped$sex = "1"
# Assign 2 to animals with F as their last letter
ped[substr(ped$offspring,1,1) == "F",]$sex = "2"
# Use regular expressions to add a column with generation number
# Start by generations 1-9, look for rows with _ followed by any character.
# Print the character after _.
ped$Generation <- substr(str_extract(ped$offspring,"_.$"),2,2)
# Next find row with _ followed by two characters. Print two characters after _
ped[is.na(ped$Generation),]$Generation <- substr(str_extract(ped[is.na(ped$Generation),]$offspring,"_..$"),2,3)
# Make column with maximum matings for input to EVA
ped$maxmatings = 0
# Modify generation values using mapvalues from the plyr package
# Males with _3 are generation 2.
# Females with _4 are 3 and males with _5 are also 3 etc.
ped$Generation = mapvalues(ped$Generation, MoBPSGen, ve)
# change to integer
ped$Generation = as.integer(ped$Generation)
# Set maxmatings to value MaxMate for selection candidates. Set to 1 for all females in
# generation gen, and to 0 for bulls with ebv lower than mean ebv
# Make a column maxmatings. Start by assigning all animals in the current 
# generation to MaxMate
ped[ped$Generation == gen-1,]$maxmatings = MaxMate
# Assign females in current generation to 1 mating
ped[ped$Generation == gen-1 & ped$sex == "2",]$maxmatings = 1
# make matrix of ebvs.
bve <- get.bve(population, cohorts = cohorts)
# The input to EVA, evaIn contains ped and ebvs.
# Use merge to combine ped and bve.
evaIn <- merge(ped, data.frame(t(bve)), by.x = "offspring", by.y = "row.names", sort = F, all.x = TRUE)
evaIn$maxmatings=as.character(evaIn$maxmatings)
# Apply truncation selection here on the GEBVs, only 20% best bulls are selecte
# This simply takes the last breedsize/4 number of bulls in the evaIn file (which should be the genotyped bulls)
evaIn[(dim(evaIn)[1]-breedSize/4):(dim(evaIn)[1]),]$maxmatings = MaxMate

#Big <- merge(ped, data.frame(t(pheno)), by.x = "offspring", by.y = "row.names", sort = F)
tail(evaIn)
#evaIn=evaIn[order(evaIn$Generation),]
evaIn$ID = mapvalues(evaIn$offspring,matc$ID,matc$RecodeID, warn_missing = FALSE)
evaIn=evaIn[,c(8,2,3,4,5,6,7,1)]

#make dataframe of individual MoBPS IDs and IDs for EVA.
evaIn$mother <- mapvalues(evaIn$mother,matc$ID,as.integer(matc$RecodeID), warn_missing = FALSE)
evaIn$father <- mapvalues(evaIn$father,matc$ID,matc$RecodeID, warn_missing = FALSE)

colnames(evaIn)=c("ID", "sire", "dam", "sex" ,"generation", "maxmatings", 
                    "ebv","MoBPSID")
evaIn$ebv = round(evaIn$ebv, 2)
#evaIn[is.na(evaIn$ID),]$ID = 0
evaIn$ID=as.integer(evaIn$ID)
evaIn=evaIn[order(as.integer(evaIn$generation),as.integer(evaIn$ID)),]
evaIn[evaIn$maxmatings != 0 & evaIn$sex==1,]$maxmatings=1
evaIn2=evaIn[evaIn$maxmatings != 0 & evaIn$sex==1,]
meanEBV=mean(evaIn[evaIn$sex==2& evaIn$maxmatings>0,]$ebv)
evaIn2=rbind(evaIn2,c("00001","0","0",2,gen-1,args[3],meanEBV,"PseudoFemale"))
write.table(evaIn2, "evaIn.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(evaIn[,c(1,4,6)], "SelCands",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

# Write genotypes if marker-based GMATRIX is used
  if (substr(args[2],1,1) == "M") { 
    # Here I retrieve the genotypes of desired cohorts
    genos = t(get.geno(population, 
                       cohorts = cohorts[-c(2,4,6,8,10)],
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
    # Should this code not be run? I don't remember what it does
    # genotyped = get.genotyped(population, gen=1:length(population$breeding))
    # for (i in 1:dim(idfile)[1]) {
    #   if (genotyped[i]) {
    #     idfile[i,3] <- 1
    #   }
    # }
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
      write.table(freq[freq!=0], 
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
      write.table(freq[freq!=0], 
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
      write.table(freq[freq!=0], 
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
a <- get.pedigree(population, cohorts = get.cohorts(population)[c(length(get.cohorts(population))-1,length(get.cohorts(population))-2)], raw=TRUE)
png(paste(args[2],"_",args[2],"_",args[1],gen,".png",sep=""))

hist(a[,6], nclass=500, xlab="sire nr.", ylab="times used", main="Frequency of use for each sire")
dev.off()
	print(paste("write genomic information for generation",gen, sep = " "))
        df=data.frame(get.effect.freq(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2])))
        # info$real.bv.add : Lists with an overview of all single marker QTLs for each trait
        head(population$info$real.bv.add[[1]][,c(1,3,4,5)])
        # take out the effects
        effects = data.frame(population$info$real.bv.add[[1]][,c(1,2,3,4,5)])
        head(effects)
        # change to character
        effects$X1=effects$X1 + (effects$X2-1)*numMarkers
        effects$X1=as.character(effects$X1)
        effects = effects[,-2]
        # find the frequencies.
        head(df)
        df$X1 = row.names(df)
        # merge 
        mrg = merge(effects, df, by = "X1")
        mrg$freq1 = (mrg$Homo0*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
        mrg$freq2 = (mrg$Homo1*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
        colnames(mrg) = c("SNP","Effect0", "EffectHet","Effect1", "Homo0","Hetero","Homo1","freq0","freq1")
        p[gen,] <- mrg$freq0
        info[gen,6] <- var(p[gen,]-p[gen-1,])
        # print mean heterozygosity and store in info matrix.
        print(c("mean heterozygosity of QTL, gen:", paste(gen), paste(mean(sum(mrg$Hetero)/(sum(mrg$Homo0+mrg$Hetero+mrg$Homo1))))))
        info[gen, 3] = mean(sum(mrg$Hetero)/(sum(mrg$Homo0+mrg$Hetero+mrg$Homo1)))
        freq1=data.frame(mrg[,c(2,8)])
        freq2=data.frame(mrg[,c(4,9)])
        afi <- apply(mrg[,c(5,6,7)], MARGIN = 1, FUN = funcSegQTL)
        sum(afi)
        info[gen, 5] <- sum(afi)/qtl
        # compute total heterozygosity of markers and number of fixed alleles
        # Get genotypes: 
        bl=get.geno(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]))
        # Assign
        fixSeg <- apply(bl, MARGIN = 1, FUN = func)
        # Sum of fixSeg is the number of segregating alleles
        sum(fixSeg)
        # Put in percentage of alleles segregating.
        info[gen,4] <- sum(fixSeg)/length(fixSeg)
        colnames(freq1) = c("effect","frequency")
        colnames(freq2) = c("effect","frequency")
        dat = rbind(freq1,freq2)
        #plot(dat$frequency,dat$effect)
        # hist(dat$frequency, breaks = 100)
        # Get kinship from IBD. 
        info[gen,2] <- kinship.emp.fast(population = population,
                                              cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]),
                                              ibd.obs = breedSize/10,
                                              hbd.obs = breedSize/10)[1]
        bv <- get.bv(population, cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]))
        info[gen,1] <- mean(bv)


# Update cohortlist
cohorts[cohIndex] =  paste("ssGBLUP",gen,"_M", sep="")
cohorts[cohIndex+1] =  paste("ssGBLUP",gen,"_F", sep="")
cohIndex =+ cohIndex+2
print(paste("end of generation",gen, sep = " "))      
}
summary(population)
RDataFile <- paste(args[2],"_",args[1],".RData",sep="")
save.image(file=RDataFile)

png(paste(args[2],"_",args[1],"_BVdev.png",sep=""))
#bv.development(population,cohort=get.cohorts(population),display.cohort.name = TRUE)
bv.development(population,gen=1:length(population$breeding),display.cohort.name = TRUE)
dev.off()
png(paste(args[2],"_",args[1],"_kinshipdev.png",sep=""))
kinship.development(population, gen=1:length(population$breeding),display.cohort.name = TRUE, ibd.obs = 10000)
dev.off()

png(paste(args[2],"_",args[1],"_BV.png",sep=""))
plot(info$BV)
dev.off()

png(paste(args[2],"_",args[1],"_Coancestry.png",sep=""))
plot(info$Coancestry)
dev.off()

png(paste(args[2],"_",args[1],"_Heteroz.png",sep=""))
plot(info$Heteroz)
dev.off()

png(paste(args[2],"_",args[1],"_SegAlleles.png",sep=""))
plot(info$SegQTL)
dev.off()

png(paste(args[2],"_",args[1],"_DriftVariance.png",sep=""))
plot(info$DriftVar)
dev.off()


qt <- get.qtl.variance(population, gen = 1:length(population$breeding))
qt=data.frame(qt)
tail(qt)
head(qt)
write.table(qt, paste(args[2],"_",args[1],"_qt.txt",sep=""), quote = F, sep = "\t")
write.table(info, file = paste(args[2],"_",args[1],"_info.txt",sep=""))

write.table(p, file = paste(args[2],"_",args[1],"_p.txt",sep=""))
