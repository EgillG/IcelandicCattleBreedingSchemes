# Updated 11.08.2021 with numerous fixes
# Can use marker-based matrices and pedigree
# Fixed the funcSegQTL function
# Uses bulls equally in each generation.
# Writes relevant parameters (some of them) to a matrix, info.
# SegAlleles currently does not work
# Updated 21.04.2022
# Removed useless statistics from the info matrix.
# Added variance of genetic values (genetic variance) to the info matrix
# Fixed the genic variance in PBLUP generations
# Now BLUP uses genetic variance as variance component, instead of genic.
# M1D and M2D now uses generation 2 as base, to reduce the effect of sampling from generation 1 to 2.
# F_hom and Fdrift computed with generation 2 as base in PBLUP and Pblupgen as base in genomic BLUP.

# Updated 30th May 2022.
# Fixed the genetic variance for simulation of phenotypes and breeding value estimation.
# Now, the genetic variance of the base population is used for both.

# Update 1st June 2022
# Fixed an issue with writing out all genotypes for OCS. 
# ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
# The program does not work if MoBPS generations are more than 89. Need to put a test to check that gens < 89

library(MoBPS)
library(stringr)
library(miraculix)
library(RandomFieldsUtils)
library(plyr)
library(dplyr)

options(scipen = 100)
RFoptions(warn_address=FALSE)

# library(tidyverse)
# library(plyr)
# Set working drive
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
# args are: no. of repetition, method, number of selected bulls/generation
print(args)
# Number of single step GBLUP evaluations
numSS = 3
# Truncation of generations for GBLUP
GBLUPgens = 9
# Proportion of males genotyped
# Number of offspring each generation
breedSize = 12000

GenoSize = 2000
#GenoSize = breedSize/2*0.25
NeuMarkMAF <- 0.001
# load QMSim data. 
#Number of QTL
qtl = 3000
# Number of males and females selected each generation
SelMales=as.integer(args[3])
SelFemales=breedSize/2
#heritability
h2 = 0.4
#Number of neutral markers
NeutralMarkers = 3000
# Variable of maximum number of matings per bull (for OCS).
#MaxMate=2000
MaxMate=breedSize/SelMales
# number of generations for PBLUP
Pblupgen = 5
# Generations for ssGBLUP and GBLUP
ssGBLUPgen = 16

print(paste("Running simulation with ", Pblupgen, " generations of PBLUP and ", ssGBLUPgen," generations of genomic selection.",sep=""))
print(paste("The number of animals simulated each generation is ", breedSize, " and the number of genotyped males is ", GenoSize , 
    ". The number of generations to truncate for GBLUP is ", GBLUPgens, ". The number of generations using single step evaluation is ", numSS,sep=""))
# Matrix for storing results of genetic gain, heterozygosity and inbreeding.
info = matrix(nrow = (Pblupgen + ssGBLUPgen), ncol  = 15)
info=data.frame(info)
colnames(info) = c("BV", "Kinship","SelfRel","Heteroz", "SegAlleles", "SegQTL", "SegNeutral","F_drift_n","F_hom_n","F_drift_qtl","F_hom_qtl","F_drift_m","F_hom_m","varA","GeneticVarA")
source("../Functions_MoBPS_simulations_20_01_2022.R")

RDatafile = paste("../../QMSim_11_02_2022/data_",args[1],".RData",sep="")
load(RDatafile)

cohorts = c(length = (Pblupgen+ssGBLUPgen)*2+2)
cohorts[1] = "Cohort_1_M"
cohorts[2] = "Cohort_1_F"
# Set index for cohorts matrix
cohIndex = 3


# Matrix for storing allele frequency changes 
p_qtl = matrix(nrow = (Pblupgen + ssGBLUPgen), ncol  = qtl)
p_marker = matrix(nrow = (Pblupgen + ssGBLUPgen), ncol  = sum(markerIncluded))
p_neutral = matrix(nrow = (Pblupgen + ssGBLUPgen), ncol  = sum(NeuList))

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

for(gen in 1:(Pblupgen)){
  #RDataFile <- paste(args[2],"_",args[1],gen,"PBLUP.RData",sep="")
  #save.image(file=RDataFile)
    print(paste("PBLUP Generation: ", gen))
    print("phenotyping cows")
if (gen==1) {
varABase = sd(get.bv(population = population, cohorts=c("Cohort_1_M", "Cohort_1_F")))**2
varEBase = ((1-h2)/h2)*varABase
}
print(paste("Using a residual variance of:", round(varEBase,3), sep = " "))
print(paste("Using a genetic variance of:", round(varABase,3),"for a heritability of", h2, sep=" "))
    population <- breeding.diploid(population,
                                   sigma.e = varEBase,
                                   phenotyping.cohorts = c(cohorts[cohIndex-1]),
                                   bve = FALSE, # estimate breeding values
                                   mutation.rate = 0,
                                   remutation.rate = 0,
                                   share.genotyped = 0)
   print(paste("phenotypes for cohort", cohorts[cohIndex-1],"were generated for generation ",gen))

    print("Running Pedigree based evaluation in DMU")
ped <- data.frame(get.pedigree(population, cohorts = cohorts))
ped2 = data.frame(get.pedigree(population, cohorts = cohorts, id = T))
matc=data.frame(ped$offspring, ped2$offspring)
colnames(matc) = c("ID","RecodeID")
matc$RecodeID <- as.integer(matc$RecodeID)
rm(ped,ped2)
pedig <- data.frame(get.pedigree(population, cohorts = cohorts))
# Put in generation number (not exactly correct but good enough):
pedig$Generation <- substr(str_extract(pedig$offspring,"_.$"),2,2)
pedig[is.na(pedig$Generation),]$Generation <- substr(str_extract(pedig[is.na(pedig$Generation),]$offspring,"_..$"),2,3)
pedig$Generation <- mapvalues(pedig$Generation, MoBPSGen, ve, F)

pedig$offspring <- mapvalues(pedig$offspring, matc$ID, matc$RecodeID, warn_missing = F)
pedig$father <- mapvalues(pedig$father, matc$ID, matc$RecodeID, warn_missing = F)
pedig$mother <- mapvalues(pedig$mother, matc$ID, matc$RecodeID, warn_missing = F)
pedig <- unique(pedig)
pedig <- pedig[pedig$father!=pedig$offspring,]
if (gen > 1) {
	pedig[substr(pedig$father,1,1)=="M",]$father <- 0
	pedig[substr(pedig$mother,1,1)=="M",]$mother <- 0
             }
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
pheno[is.na(pheno)] <- -9999
pheno$Trait.1 <- round(pheno$Trait.1,4)
pheno = pheno[,c(2,3,1)]
write.table(pheno, "phenotypes.txt",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE, 
            sep = " ")
rm(pheno)
# write the necessary variances (I use the true values)
   cohortsAll <- get.cohorts(population)
   lenCo <- length(get.cohorts(population))
# Use genetic variance of the base for breeding value estimation.
varABase = sd(get.bv(population = population, cohorts=c("Cohort_1_M", "Cohort_1_F")))**2
varEBase = ((1-h2)/h2)*varABase
write.table(rbind(c(1,1,1, varABase),c(2,1,1,varEBase)),
                "variances.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = F,
                sep = " ")

system("bash r_dmu5 dmuPBLUP")
system("awk '{print $5,$8}' dmuPBLUP.SOL > temp")
dmuBVE = data.frame(read.table("temp", skip = 2, stringsAsFactors = F))
dmuBVE$V1 <- mapvalues(dmuBVE$V1, matc$RecodeID, matc$ID, F)
population <- insert.bve(population, dmuBVE)
print("generating a copy of selected bulls")
   # Generate a copy of the selected bulls.
   population <- breeding.diploid(population, 
                                  selection.size = c(SelMales*3,SelFemales),
                                  copy.individual.m = TRUE,
                                  share.genotyped = 0,                                  #add.gen = gen-1,
                                  selection.criteria = "bve",
                                  selection.f.cohorts = cohorts[cohIndex-1],
                                  selection.m.cohorts = cohorts[cohIndex-2],
                                  mutation.rate = 0,
                                  remutation.rate = 0,
                                  name.cohort = paste("SelectedBulls",gen, sep=""))
   # Genotype the copies.
   population <- breeding.diploid(population,
                                  genotyped.cohorts = paste("SelectedBulls",gen, sep=""),
                                  genotyped.array = 2,
                                  mutation.rate = 0,
                                  remutation.rate = 0
                                  )  
   print(paste(cohorts[cohIndex-2]," were used for generating copies to genotype: ","SelectedBulls",gen,sep=""))
   
   # Apply selection, here I have to set the number of matings so that females
   # have 1 offspring each
   population <- breeding.diploid(population, 
                                 breeding.size = breedSize,
                                 name.cohort = paste("PBLUP",gen+1, sep=""),
                                 selection.size = c(3*SelMales,SelFemales),
                                 #heritability = 0.4,
                                 miraculix.cores = 4,
                                 selection.f.cohorts = cohorts[cohIndex-1],
                                 selection.m.cohorts = paste("SelectedBulls",gen, sep=""),
                                 selection.criteria = "bve",
                                 share.genotyped = 0,
                                 max.mating.pair = 2,
                                 max.offspring = c(breedSize/SelMales,2),
                                 mutation.rate = 0,
                                 remutation.rate = 0,
                                 store.effect.freq = TRUE)
population <- clean.up(population)
#Genotype females of the last PBLUP generation
if (gen == Pblupgen) {
population <- breeding.diploid(population,
                               genotyped.cohorts = cohorts[cohIndex-1], #genotype all females of last generation
                               mutation.rate = 0,
                               genotyped.array = 2,
                               remutation.rate = 0,
                               bve = FALSE)
                               }
# This is to genotype generation = gen (all selection candidates)
   cohortsAll <- get.cohorts(population)
   lenCo <- length(get.cohorts(population))
# Start obtaining genomic information
   mrg <- getHet(population, c(cohortsAll[lenCo-3], cohortsAll[lenCo-4]),qtl)
   p_qtl[gen,] <- mrg$freq0
   info[gen,1:7] <- getAll(mrg, gen, c(cohortsAll[lenCo-3], cohortsAll[lenCo-4]), population, funcSegQTL, func)
   p_marker[gen,] <- get_p(population, c(cohortsAll[lenCo-3], cohortsAll[lenCo-4]), markerIncluded) 
   p_neutral[gen,] <- get_p(population, c(cohortsAll[lenCo-3], cohortsAll[lenCo-4]), NeuList) 

   info[gen,8] = F_drift(p_neutral, gen, 2)
   info[gen,9] = F_hom(p_neutral, gen, 2)
   info[gen,10] = F_drift(p_qtl, gen, 2)
   info[gen,11] = F_hom(p_qtl, gen, 2)
   info[gen,12] = F_drift(p_marker, gen, 2)
   info[gen,13] = F_hom(p_marker, gen, 2)

   info[gen,14] = sum(get.qtl.variance(population = population, cohorts=c(cohortsAll[lenCo-3], cohortsAll[lenCo-4]))[[1]][,3])
   info[gen,15] = sd(get.bv(population = population, cohorts=c(cohortsAll[lenCo-3], cohortsAll[lenCo-4])))**2

   print(round(info,4))

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
cohorts
gen=Pblupgen+1
print("Info matrix looks like this before ssGBLUP stage:\n")
print(info)
print("starting genomic selection")

# 
###############################################################################
# Start of genomic selection
for (gen in (Pblupgen+1):(ssGBLUPgen+Pblupgen)){
ptmLong <- proc.time()
#RDataFile <- paste(args[2],"_",args[1],gen,"GenomicSel.RData",sep="")
#save.image(file=RDataFile)

print(paste("list of cohorts when genomic selection starts in generation ",gen, ":", toString(cohorts), sep=""))
# Here select 2000 bulls based on parent average to be genotyped.
# Check whether the parent mean is computed correctly (for the selection candidates)
population <- breeding.diploid(population,
                               bve.parent.mean = TRUE,
                               bve.insert.cohorts = cohorts[cohIndex-2],
                               mutation.rate = 0,
                               remutation.rate = 0
                                  )
# Now use the parent mean to preselect bulls for genotyping
population <- breeding.diploid(population,
                               selection.size = GenoSize,
                               copy.individual.m = TRUE,
                               selection.criteria = "bve",
                               selection.m.cohorts = cohorts[cohIndex-2],
                               mutation.rate = 0,
                               genotyped.array = 2,
                               remutation.rate = 0,
#                              add.gen = Pblupgen+gen-2, #check if this is correct.
                               name.cohort = paste("GenoBulls",gen, sep=""))
cohorts[cohIndex-2] <- paste("GenoBulls",gen, sep="")
print(paste("her er eg a undan ebv, kynslod: ",gen))

# Female parents of selection candidates are phenotyped, but not males.
# All females in the current generation are genotyped at this stage
print(paste("Using a residual variance of:", varEBase,"to simulate phenotypes", sep = " "))
population <- breeding.diploid(population,
                               phenotyping.cohorts = cohorts[cohIndex-3], #phenotype previous generation
                               genotyped.cohorts = c(cohorts[cohIndex-1], cohorts[cohIndex-2]), #genotype all females, and 2000 bulls
                               genotyped.array = 2,
                               sigma.e = varEBase,
                               mutation.rate = 0,
                               remutation.rate = 0,
                               bve = FALSE)
print(paste("her er eg aftur eftir ebv, kynslod: ",gen))

# Now use DMU for single step breeding value estimation.
## First write necessary files.
ped <- data.frame(get.pedigree(population, cohorts = get.cohorts(population)))
ped2 = data.frame(get.pedigree(population, cohorts = get.cohorts(population), id = T))
matc=data.frame(ped$offspring, ped2$offspring)
colnames(matc) = c("ID","RecodeID")
matc$RecodeID <- as.integer(matc$RecodeID)
rm(ped,ped2)

pedig <- data.frame(get.pedigree(population, cohorts = cohorts))
# Put in generation number (not exactly correct but good enough):
pedig$Generation <- substr(str_extract(pedig$offspring,"_.$"),2,2)
pedig[is.na(pedig$Generation),]$Generation <- substr(str_extract(pedig[is.na(pedig$Generation),]$offspring,"_..$"),2,3)
pedig$Generation <- mapvalues(pedig$Generation, MoBPSGen, ve, F)

pedig$offspring <- mapvalues(pedig$offspring, matc$ID, matc$RecodeID, warn_missing = F)
pedig$father <- mapvalues(pedig$father, matc$ID, matc$RecodeID, warn_missing = F)
pedig$mother <- mapvalues(pedig$mother, matc$ID, matc$RecodeID, warn_missing = F)
pedig <- unique(pedig)
#pedig <- pedig[pedig$father!=pedig$offspring,]
pedig[pedig$father==pedig$offspring,]$father <- 0
pedig[pedig$mother==pedig$offspring,]$mother <- 0
write.table(pedig, 
            "pedigree.txt",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE, 
            sep = " ")


# Variable for cohorts to use. Uses truncation after GBLUPgens is passed

if ((gen-Pblupgen-1) < GBLUPgens) {
     cohorts_to_use <- cohorts[-seq(2, to = (Pblupgen-1)*2, by = 2)]
} else {
     print(paste("Start truncation of genomic information in generation: ", gen,". Using ",GBLUPgens," generations to form G-matrix."))
     cohorts_to_use <- cohorts[(cohIndex-2*GBLUPgens):(cohIndex-1)]
}
print("Writing phenotypes.")
pheno = data.frame(t(get.pheno(population, cohorts = cohorts_to_use)))
pheno$ID <- rownames(pheno)
pheno$ID <- mapvalues(pheno$ID, matc$ID, matc$RecodeID, F)
pheno$mean = 1
pheno$Trait.1 <- round(pheno$Trait.1,4)
pheno = pheno[,c(2,3,1)]
write.table(na.omit(pheno), "phenotypes.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = " ")
rm(pheno)

print("Writing and transposing genotype data.")
ptm <- proc.time()
#system("rm Gmatrix/gmat.dat Gmatrix/gmat.id")
if (gen==Pblupgen+1) {
    coSel <- 1:length(cohorts_to_use)
    system("rm Gmatrix/gmat.dat Gmatrix/gmat.id")
    } else if ((gen-Pblupgen-1) >= GBLUPgens) {
    # use tail to select the individuals included in cohorts_to_use, minus the current generation.

    argument <- paste("head -n -", breedSize," Gmatrix/gmat.dat > Gmatrix/temp", sep="")
    system(argument)
    system("mv Gmatrix/temp Gmatrix/gmat.dat")
    argument <- paste("head -n -", breedSize," Gmatrix/gmat.id > Gmatrix/temp", sep="")
    system(argument)
    system("mv Gmatrix/temp Gmatrix/gmat.id")

    argument <- paste("tail -n ", dim(get.bv(population, cohorts = cohorts_to_use))[2]-GenoSize-breedSize/2-as.numeric(args[3])-breedSize/2," Gmatrix/gmat.dat > Gmatrix/temp", sep="")
    system(argument)
    system("mv Gmatrix/temp Gmatrix/gmat.dat")
    argument <- paste("tail -n ", dim(get.bv(population, cohorts = cohorts_to_use))[2]-GenoSize-breedSize/2-as.numeric(args[3])-breedSize/2," Gmatrix/gmat.id > Gmatrix/temp", sep="")
    system(argument)
    system("mv Gmatrix/temp Gmatrix/gmat.id")
    system("wc -l Gmatrix/gmat.id | echo")
    print("Number of lines in Gmatrix/gmat.id:")
    print(system("wc -l Gmatrix/gmat.id | cat"))

    coSel <- c(length(cohorts_to_use)-3,length(cohorts_to_use)-2,length(cohorts_to_use)-1,length(cohorts_to_use))
    } else {
    coSel <- coSel <- c(length(cohorts_to_use)-3,length(cohorts_to_use)-2,length(cohorts_to_use)-1,length(cohorts_to_use))
    argument <- paste("head -n -", breedSize," Gmatrix/gmat.dat > Gmatrix/temp", sep="")
    system(argument)
    system("mv Gmatrix/temp Gmatrix/gmat.dat")
    argument <- paste("head -n -", breedSize," Gmatrix/gmat.id > Gmatrix/temp", sep="")
    system(argument)
    system("mv Gmatrix/temp Gmatrix/gmat.id")
    }

for (co in cohorts_to_use[coSel]) {
    print(paste("Writing genotypes for cohort",co))
    genos = data.frame(t(get.geno(population, 
                   cohorts = co,
                   non.genotyped.as.missing = TRUE)))
    genos <- genos[,markerIncluded]
    genos$ID <- row.names(genos)
    genos$RecodeID <- mapvalues(genos$ID, matc$ID, matc$RecodeID, F)
    genos <- genos %>% relocate(RecodeID,ID)
    write.table(genos,
            "Gmatrix/gmat.dat", 
            append = TRUE,
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE, 
            sep = " ")
    # make mapfile for GMATRIX
    # mapfile is simply the number of snps and a 0 or 1, depending on whether
    # the marker is genotyped or not
      mapfile = cbind(seq(1,sum(markerIncluded)),
                      rep(1,sum(markerIncluded)))
    # Make id file for GMATRIX
    idfile <- cbind(genos$RecodeID,genos$RecodeID,rep(1,length(genos$RecodeID)))
    rm(genos)
    write.table(idfile, 
                "Gmatrix/gmat.id",
                append = TRUE,
     	           row.names = FALSE, 
                col.names = FALSE,
                quote = F,
                sep = " ")
    }          
print(paste("Writing of genotype data finished in ", (proc.time()-ptm)[1], " seconds, in generation ", gen, sep=""))

print("Memory use of R for objects:")
print(sort( sapply(ls(),function(x){object.size(get(x))}))[(length(sapply(ls(),function(x){object.size(get(x))}))-5):length(sapply(ls(),function(x){object.size(get(x))}))])	
write.table(cbind(mapfile,rep(1,dim(mapfile)[1])),
                "Gmatrix/gmat.map",
                row.names = FALSE,
                col.names = FALSE,
                quote = F,
                sep = " ")

# write the necessary variances (I use the true values)
   cohortsAll <- get.cohorts(population)
   lenCo <- length(get.cohorts(population))

varABase = sd(get.bv(population = population, cohorts=c("Cohort_1_M", "Cohort_1_F")))**2
varEBase = ((1-h2)/h2)*varABase
write.table(rbind(c(1,1,1, varABase),c(2,1,1,varEBase)),
                "variances.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = F,
                sep = " ")

singleS = "No"
if (gen < (Pblupgen+1+numSS)) {
      singleS = "Yes"
      }
system(paste("python3 gmatPar.py ", getwd(), " M1 ", singleS,  sep = ""))
ptm <- proc.time()

system(paste("bash DMU.sh", args[2], singleS, args[1], sep=" "))
print(paste("Setting up of GRM and genetic evaluation completed in ", (proc.time()-ptm)[3]/60, " minutes, in generation ", gen, sep=""))

#DMU ebvs are inserted back into MoBPS and extracted again for ocs computations (silly indeed).
dmuBVE = data.frame(read.table("dmuSS.SOL", skip = 2, stringsAsFactors = F))
#dmuBVE$V2 <- as.integer(dmuBVE$V2)

# I have to remove some individuals from matc, to insert correct IDs.
if (((Pblupgen + 1) * 3 + 1) < length(get.cohorts(population))) {
  matcNotUse <- rownames(data.frame(get.id(population, cohorts = get.cohorts(population)[ c(seq(1,(Pblupgen+1)*3, by = 3),seq((Pblupgen+1)*3+1,length(get.cohorts(population)), by = 4))])))

} else {
matcNotUse <- rownames(data.frame(get.id(population, cohorts = get.cohorts(population)[seq(1,(Pblupgen+1)*3, by = 3)])))
}

matc_reduced <- matc[!(matc$ID %in% matcNotUse),]
dmuBVE$V1 <- mapvalues(dmuBVE$V1, matc_reduced$RecodeID, matc_reduced$ID, F)

population <- insert.bve(population, dmuBVE)
#rm(dmuBVE)
###############################################################################
# Try OCS
###############################################################################
if (gen==Pblupgen+1) {
  tempco <- c(paste("PBLUP",gen,"_M",sep=""),paste("PBLUP",gen,"_F",sep=""))
  } else {
  tempco <- c(paste("ssGBLUP",gen,"_M",sep=""),paste("ssGBLUP",gen,"_F",sep=""))
  }

  ped <- data.frame(get.pedigree(population, cohorts = tempco))
  # Assign sex
  # First assign all animals as 1 (males)
  ped$sex = "1"
  # Assign 2 to animals with F as their first letter
  ped[substr(ped$offspring,1,1) == "F",]$sex = "2"
  bve <- get.bve(population, cohorts = cohorts[c(cohIndex-1,cohIndex-2)])
  rusl <- merge(matc, data.frame(t(bve)), by.x = "ID", by.y = "row.names", sort = F, all.x = TRUE)

  ped$offspring <- mapvalues(ped$offspring, matc$ID, matc$RecodeID, warn_missing = F)
  ped$mother <- mapvalues(ped$mother, matc$ID, matc$RecodeID, warn_missing = F)
  ped$father <- mapvalues(ped$father, matc$ID, matc$RecodeID, warn_missing = F)

  # The input to EVA, evaIn contains ped and ebvs.
  # Use merge to combine ped and bve.
  evaIn <- merge(ped, rusl, by.x = "offspring", by.y = "RecodeID", sort = F, all.x = TRUE)
  rm(ped,bve,rusl)

  #make dataframe of individual MoBPS IDs and IDs for EVA.
  evaIn$mother <- 0
  evaIn$father <- 0
  evaIn$generation = 1
  evaIn$maxmatings = 1
  evaIn[is.na(evaIn$Trait.1),]$maxmatings=0
  evaIn[is.na(evaIn$Trait.1),]$Trait.1=-900
  evaIn[evaIn$Trait.1<median(evaIn[evaIn$sex==1 & evaIn$maxmatings==1,]$Trait.1) | evaIn$sex==2,]$maxmatings <- 0

  #Only set bulls with EBV above the median as selection candidates.
  evaIn[evaIn$sex==2,]$maxmatings = args[3]
  evaIn = evaIn[c(1,2,3,4,7,8,6,5)]
  colnames(evaIn)=c("ID", "sire", "dam", "sex" ,"generation", "maxmatings", "ebv","MoBPSID")
  evaIn$ebv = round(evaIn$ebv, 4)
  meanEBV=mean(evaIn[evaIn$sex== 2 & evaIn$maxmatings>0,]$ebv)  
  # Remove the duplicated IDs of the genotyped animals. Make table with the duplicate IDs and put into variable.
  tavle <- table(evaIn$ID)
  Genotyped <- names(tavle[tavle==2])
  if (length(population$breeding)> 9) {
        higherGen <- max(substr(str_extract(evaIn$MoBPSID,"_..$"),2,3),na.rm=T)
        ru <- (!(evaIn$ID %in% Genotyped& substr(str_extract(evaIn$MoBPSID,"_..$"),2,3)<higherGen))
        ru[is.na(ru)] <- FALSE 
        evaIn <- evaIn[ru,]
      } else {
        higherGen <- max(substr(str_extract(evaIn$MoBPSID,"_.$"),2,2),na.rm=T)
        evaIn <- evaIn[!(evaIn$ID %in% Genotyped & substr(str_extract(evaIn$MoBPSID,"_.$"),2,2)<higherGen),]
      }
  # Write selection candidates before adding pseudofemale
  write.table(evaIn[,c(1,4)], "SelCands",
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  evaIn=rbind(evaIn[evaIn$sex==1,],c("00001","0","0",2,1,args[3],round(meanEBV,4),"PseudoFemale"))
  write.table(evaIn, "evaIn.txt",
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  rm(evaIn)
  if (args[2]=="Ped") {
  print("Writing Pedigree information")
# Select the cohorts for which to write pedigree information
if (gen==Pblupgen+1) {
  tempco <- c(paste("PBLUP",gen,"_M",sep=""),paste("PBLUP",gen,"_F",sep=""))
  } else{
  tempco <- c(paste("ssGBLUP",gen,"_M",sep=""),paste("ssGBLUP",gen,"_F",sep=""))
  }
  bla = kinship.exp.store(population, cohorts = tempco, depth.pedigree = GBLUPgens)
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
  rm(bla)
  write.table(newFile,"ReducedMatrix", quote = F, row.names = F, col.names = F)
  rm(newFile)
  }
# Write genotypes if marker-based GMATRIX is used
  if (substr(args[2],1,1) == "M" | args[2]=="Ped") { 
    # Here I retrieve the genotypes of desired cohorts
# I have to remove the genotyped bulls + females from gmat.datm and replace with all bulls and females
# Use head to remove the bulls and females
RemoveBulls = breedSize/2+GenoSize
argument <- paste("head -n -", RemoveBulls," Gmatrix/gmat.dat > Gmatrix/temp", sep="")
system(argument)
system("mv Gmatrix/temp Gmatrix/gmat.dat")
argument <- paste("head -n -", RemoveBulls," Gmatrix/gmat.id > Gmatrix/temp", sep="")
system(argument)
system("mv Gmatrix/temp Gmatrix/gmat.id")

# Now append the current generation 6000 males and 6000 females to the genotype data file
if (gen==Pblupgen+1) {
tempco <- c(paste("PBLUP",gen,"_M",sep=""),paste("PBLUP",gen,"_F",sep=""))
} else{
tempco <- c(paste("ssGBLUP",gen,"_M",sep=""),paste("ssGBLUP",gen,"_F",sep=""))
}
for (co in tempco) {
print(paste("Writing genotypes for cohort",co))
genos = data.frame(t(get.geno(population,
                   cohorts = co,
                   non.genotyped.as.missing = FALSE)))
genos <- genos[,markerIncluded]
genos$ID <- row.names(genos)
genos$RecodeID <- mapvalues(genos$ID, matc$ID, matc$RecodeID, F)
genos <- genos %>% relocate(RecodeID,ID)
write.table(genos,
            "Gmatrix/gmat.dat",
            append = TRUE,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = " ")
    # Make id file for GMATRIX
    idfile <- cbind(genos$RecodeID,genos$RecodeID,rep(1,length(genos$RecodeID)))
    rm(genos)
    write.table(idfile,
                "Gmatrix/gmat.id",
                append = TRUE,
                row.names = FALSE,
                col.names = FALSE,
                quote = F,
                sep = " ")
}

print("Memory use of R for objects:")
print(sort( sapply(ls(),function(x){object.size(get(x))}))[(length(sapply(ls(),function(x){object.size(get(x))}))-5):length(sapply(ls(),function(x){object.size(get(x))}))])	
rm(genos)

    if (args[2] == "M1D" || args[2] == "M2D") {
      # find the frequency of alleles in the base generation.
      freqfile <-  get.geno(population, cohorts = c("PBLUP2_M","PBLUP2_F"), non.genotyped.as.missing = FALSE)
      freqfile <- freqfile[markerIncluded,]
    # rowMeans can be used to compute allele frequencies (for the larger allele)
      freqfile <- rowMeans(freqfile)/2
      freqfile <- cbind(seq(1,length(freqfile)),freqfile)
      write.table(freqfile, 
                  "mapbase.dat",
                  row.names = FALSE, 
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }
    if (args[2] == "M1D_Est" || args[2] == "M2D_Est") {
      # Now use the genotypes of older bulls. This uses the first generation of selected bulls.
      freqfile <-  get.geno(population, 
                            cohorts = c(cohorts[c(1,3,5,7)]),
                            non.genotyped.as.missing = FALSE)
      freqfile <- freqfile[markerIncluded,]
      # rowMeans can be used to compute allele frequencies (for the larger allele)
      freqfile <- rowMeans(freqfile)/2
      freqfile <- cbind(seq(1,length(freqfile)),freqfile)
      write.table(freqfile, 
                  "mapbase.dat",
                  row.names = FALSE, 
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }
    
    if (args[2] == "M1R" || args[2] == "M2R") {
      freqfile <-  get.geno(population, 
                            cohorts = c(cohorts[cohIndex-1],cohorts[cohIndex-2]),
                            non.genotyped.as.missing = FALSE)
      freqfile <- freqfile[markerIncluded,]
      # rowMeans can be used to compute allele frequencies (for the larger allele)
      freqfile <- rowMeans(freqfile)/2
      freqfile <- cbind(seq(1,length(freqfile)),freqfile)
      write.table(freqfile, 
                  "mapbase.dat",
                  row.names = FALSE, 
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }
    if (args[2] == "M1_05" || args[2] == "M2_05") {
      freqfile <- rep(0.5,sum(markerIncluded))
      freqfile <- cbind(seq(1,sum(markerIncluded)),freqfile)
      write.table(freqfile,
                  "mapbase.dat",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = F,
                  sep = " ")
    }
    
  }
rm(freqfile)

ptm <- proc.time()
system(paste("bash EVA.sh",args[1],args[2],args[3],sep = " "))
print(paste("Setting up of GRM and optimum contribution computations completed in ", (proc.time()-ptm)[3]/60, " minutes, in generation ", gen, sep=""))

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
rm(df)
##############################################################
# Make R Data file. This is temporary
#RDataFile <- paste(args[2],"/gen",gen,args[2],"_",args[1],"TempStuff.RData",sep="")
#save.image(file=RDataFile)
##############################################################

population <- breeding.diploid(population, 
                               breeding.size = breedSize,
                               name.cohort = paste("ssGBLUP",gen+1, sep=""),
                               selection.size = c(SelMales,SelFemales),
                               #heritability = 0.4,
                               selection.f.cohorts = cohorts[cohIndex-1],
                               selection.m.cohorts = cohorts[cohIndex-2],
                               selection.criteria = "bve",
                               max.mating.pair = 2,
                               mutation.rate = 0,
                               remutation.rate = 0,
                               share.genotyped = 0,
                               miraculix.cores = 4,
                               max.offspring = c(breedSize/SelMales,2))
population <- insert.bve(population, dmuBVE)

# Generate copies of selected bulls
population <- breeding.diploid(population,
                               selection.size = c(SelMales,SelFemales),
                               copy.individual.m = TRUE,
                               genotyped.array = 2,                                  
                               selection.criteria = "bve",
                               selection.f.cohorts = cohorts[cohIndex-1],
                               selection.m.cohorts = cohorts[cohIndex-2],
                               mutation.rate = 0,
                               add.gen = length(population$breeding)-1,
#                               delete.individual = length(population$breeding)-1,
#                               delete.sex=1,
                               remutation.rate = 0,
                               name.cohort = paste("SelectedBulls",gen, sep=""))

population <- clean.up(population)

# Here I can decide whether to make new cohort of selected bulls,
# to replace the cohort of genotyped bulls.
#a <- get.pedigree(population, gen = length(population$breeding), raw=TRUE)
#table(a[,6])	
# Replave the ebvs with dmu ones.
#rm(dmuBVE)
   print(paste("write genomic information for generation",gen, sep = " "))
   cohortsAll <- get.cohorts(population)
   lenCo <- length(get.cohorts(population))
if (gen == (Pblupgen+1)) {
   lenCo <- lenCo+1
}
   mrg <- getHet(population,c(cohortsAll[lenCo-5], cohortsAll[lenCo-6]),qtl)
   p_qtl[gen,] <- mrg$freq0
   info[gen,1:7] <- getAll(mrg, gen, c(cohortsAll[lenCo-5], cohortsAll[lenCo-6]), population, funcSegQTL, func)
   p_marker[gen,] <- get_p(population, c(cohortsAll[lenCo-5], cohortsAll[lenCo-6]), markerIncluded) 
   p_neutral[gen,] <- get_p(population, c(cohortsAll[lenCo-5], cohortsAll[lenCo-6]), NeuList) 

   info[gen,8] = F_drift(p_neutral, gen, Pblupgen)
   info[gen,9] = F_hom(p_neutral, gen, Pblupgen)
   info[gen,10] = F_drift(p_qtl, gen, Pblupgen)
   info[gen,11] = F_hom(p_qtl, gen, Pblupgen)
   info[gen,12] = F_drift(p_marker, gen, Pblupgen)
   info[gen,13] = F_hom(p_marker, gen, Pblupgen)

   info[gen,14] = sum(get.qtl.variance(population = population, cohorts=c(cohortsAll[lenCo-5], cohortsAll[lenCo-6]))[[1]][,3])
   info[gen,15] = sd(get.bv(population = population, cohorts=c(cohortsAll[lenCo-5], cohortsAll[lenCo-6])))**2
	rm(mrg)
print(info)
# Update cohortlist
cohorts[cohIndex] =  paste("ssGBLUP",gen+1,"_M", sep="")
cohorts[cohIndex+1] =  paste("ssGBLUP",gen+1,"_F", sep="")
# Replace the previous cohort of genotyped bulls with the selected bulls.
cohorts[cohIndex-2] =  paste("SelectedBulls",gen, sep="")
		
cohIndex =+ cohIndex+2
print(paste("end of generation",gen, sep = " "))      
#RDataFile <- paste(args[2],"_",args[1],gen,"PostGenomSel.RData",sep="")
#save.image(file=RDataFile)
write.table(round(info,4), file = paste(args[2],"_",args[1],"_info.txt",sep=""))
print(paste("Genetic evaluation and selection completed in ", (proc.time()-ptmLong)[3]/3600, " hours, in generation ", gen, sep=""))
}
summary(population)
RDataFile <- paste(args[2],"_",args[1],".RData",sep="")
save.image(file=RDataFile)

write.table(round(info,4), file = paste(args[2],"_",args[1],"_info.txt",sep=""))
print("Writing QTL allele frequencies")
write.table(round(p_qtl,4), file = paste(args[2],"_",args[1],"_p_qtl.txt",sep=""))
print("Writing marker frequencies")
write.table(round(p_marker,4), file = paste(args[2],"_",args[1],"_p_markers.txt",sep=""))
print("Writing neutral loci frequencies")
write.table(round(p_neutral,4), file = paste(args[2],"_",args[1],"_p_neutral.txt",sep=""))
print("Writing QTL effects")
#write.table(get.qtl.effects(population=population)[1], file = paste(args[2],"_",args[1],"_qtl_effects.txt",sep=""), row.names = F)

#qt <- get.qtl.variance(population, gen = 1:length(population$breeding))
#qt=data.frame(qt)
#tail(qt)
#head(qt)
#write.table(qt, paste(args[2],"_",args[1],"_qt.txt",sep=""), quote = F, sep = "\t")
system("rm Gmatrix/gmat.dat Gmatrix/gmat Gmatrix/invgmat Gmatrix.gmat ReducedMatrix")
