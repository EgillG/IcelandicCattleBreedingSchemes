# This is a function to check whether alleles are segregating or fixed.
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

getHet<- function(x,y,z) {
  #  This function extracts the heterozygosity of markers and returns a data frame
  df = data.frame(get.effect.freq(population=population, cohorts = y, sort = T))
  effects = data.frame(x$info$real.bv.add[[1]][,c(1,2,3,4,5)])
  mrg <- cbind(df,effects[order(effects$X2,effects$X1),],snps)[,c(10,6,7,8,1,2,3)]
  mrg$freq1 = (mrg$Homo0*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
  mrg$freq2 = (mrg$Homo1*2 + mrg$Hetero)/((mrg$Homo0+mrg$Hetero+mrg$Homo1)*2)
  colnames(mrg) = c("SNP","Effect0", "EffectHet","Effect1", "Homo0","Hetero","Homo1","freq0","freq1")
  return(mrg)
}

getAll<- function(x, y, w, p, func1, func2) {
  #  This function extracts:
  # 1 mean bv of cohorts
  # 2 mean coancestry of cohorts
  # 3 mean heterozygosity of QTL
  # 4 Percentage of segregating QTL
  # 5 Percentage of segregating genotyped markers
  # 6 Percentage of segregating neutral markers
  # 
  # 
  out=rep(0,7)
  bv <- get.bv(population, cohorts = w)
  out[1] <- mean(bv)
  # Get kinship from IBD. 
  kins <- kinship.emp.fast(population = p,
                             cohorts = w,
                             ibd.obs = (length(bv)*(length(bv)-1)/400),
                             hbd.obs = length(bv))
  out[2] <- kins[1]
  out[3] <- kins[2]
  out[4] = mean(sum(x$Hetero)/(sum(x$Homo0+x$Hetero+x$Homo1)))
  freq1=data.frame(x[,c(2,8)])
  freq2=data.frame(x[,c(4,9)])
  afi <- apply(x[,c(5,6,7)], MARGIN = 1, FUN = func1)
  out[6] <- sum(afi)/qtl
  # compute total heterozygosity of markers and number of fixed alleles
  # Get genotypes: 
  # Assign
  geno <- get.geno(p, cohorts = w)
  fixSeg <- apply(geno, MARGIN = 1, FUN = func2)
  rm(geno)
  # Sum of fixSeg is the number of segregating alleles
  # Put in percentage of alleles segregating.
  out[5] <- sum(fixSeg[markerIncluded])/length(fixSeg[markerIncluded])
  out[7] <- sum(fixSeg[NeuList])/length(fixSeg[NeuList])
  colnames(freq1) = c("effect","frequency")
  colnames(freq2) = c("effect","frequency")
  dat = rbind(freq1,freq2)
  return(out)
}

get_p <- function(p, w, List) {
  geno1 = get.geno(p, cohorts = w)
  p = (rowSums(geno1)/(2*dim(geno1)[2]))[List]
  return(p)
}

F_drift <- function(freq, gen, refGen) {
   freq <- freq[,freq[refGen,] != 0 & freq[refGen,] != 1]
   nLoci <- dim(freq)[2]
   delta <- (freq[gen,] - freq[refGen,])**2
   pq <- (freq[refGen,]*(1-freq[refGen,]))
   result <- sum(delta/pq, na.rm=T)/nLoci
   return(result)
}

F_hom <- function(freq, gen, refGen) {
  freq <- freq[,freq[refGen,] != 0 & freq[refGen,] != 1]
  nLoci <- dim(freq)[2]
  pq_t <- 2*(freq[gen,]*(1-freq[gen,]))
  pq_0 <- 2*(freq[refGen,]*(1-freq[refGen,]))
  sum <- sum(pq_t/pq_0, na.rm=T)
  result <- sum/nLoci
  return(1-result)
}
