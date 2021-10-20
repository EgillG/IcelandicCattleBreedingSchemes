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

getHet<- function(x,y,z) {
  #  This function extracts the heterozygosity of markers and returns a data frame
  df = data.frame(get.effect.freq(population=x, cohorts = y))
  effects = data.frame(x$info$real.bv.add[[1]][,c(1,2,3,4,5)])
  effects$X1=effects$X1 + (effects$X2-1)*numMarkers
  effects$X1=as.character(effects$X1)
  effects = effects[,-2]
  # find the frequencies.
  df$X1 = row.names(df)
  # merge 
  mrg = merge(effects, df, by = "X1")
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
  out=rep(0,6)
  bv <- get.bv(population, cohorts = w)
  out[1] <- mean(bv)
  # Get kinship from IBD. 
  out[2] <- kinship.emp.fast(population = p,
                             cohorts = w,
                             ibd.obs = length(bv)/10,
                             hbd.obs = length(bv)/10)[1]
  
  out[3] = mean(sum(x$Hetero)/(sum(x$Homo0+x$Hetero+x$Homo1)))
  freq1=data.frame(x[,c(2,8)])
  freq2=data.frame(x[,c(4,9)])
  afi <- apply(x[,c(5,6,7)], MARGIN = 1, FUN = func1)
  out[5] <- sum(afi)/qtl
  # compute total heterozygosity of markers and number of fixed alleles
  # Get genotypes: 
  # Assign
  fixSeg <- apply(get.geno(p, cohorts = w), MARGIN = 1, FUN = func2)
  # Sum of fixSeg is the number of segregating alleles
  # Put in percentage of alleles segregating.
  out[4] <- sum(fixSeg)/length(fixSeg)
  out[6] <- sum(fixSeg[IBDList])/length(fixSeg[IBDList])
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
