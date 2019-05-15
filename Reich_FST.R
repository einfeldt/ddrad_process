## Calculate pairwise population FST on a genind object using Reich's unbiased estimate (from: https://www.nature.com/articles/nature08365)
# Advantage over W&C estimate: no upwards bias with inbreeding, unequal sample sizes, small sample sizes (n < 10)

# code modified from: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091237#s6, supplement: https://doi.org/10.1371/journal.pone.0091237.s002)

## Function for calculating Reich FST between two population variant count dataframes
# Note: If no loci removed by call.rate being too low, returns error
Reich.Fst <- function(pop1,pop2, call.rate=0.75, top.number=10){
  # remove the SNPs that are not in common between the 2 populations
  snp.to.keep <- intersect(row.names(pop1),row.names(pop2))
  if (length(snp.to.keep) == 0){print("Error: no SNP in common");return(NULL)}
  pop1.k <- pop1[snp.to.keep,]
  pop2.k <- pop2[snp.to.keep,]
  # Commented section: unnecessary if reference allele and variant built in same order, left for modification if building population tables using a different method than below
  # change the reference allele if is not concordant between the 2 populations
  # if (sum(pop1.k$A1 == pop2.k$A1) != length(snp.to.keep)){
  #   idx <- which(pop1.k$A1 != pop2.k$A1)
  #   idx.rev <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 == pop2.k$A2)
  #   idx.rm  <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 != pop2.k$A2)
  #   if(length(idx.rev) > 0){
  #     provv <- pop1.k$A1[idx.rev]
  #     pop1.k$A1[idx.rev] <- pop1.k$A2[idx.rev]
  #     pop1.k$A2[idx.rev] <- provv
  #     provv <- pop1.k$x0[idx.rev]
  #     pop1.k$x0[idx.rev] <- pop1.k$x2[idx.rev]
  #     pop1.k$x2[idx.rev] <- provv}
  #   if(length(idx.rm) > 0){      
  #     pop1.k <- pop1.k[-idx.rm,]
  #     pop2.k <- pop2.k[-idx.rm,]}}
  # remove SNPs with low call rate in one or both populations
  x0 <- pop1.k$x0
  x1 <- pop1.k$x1
  x2 <- pop1.k$x2
  s <- x0 + x1 + x2
  y0 <- pop2.k$x0
  y1 <- pop2.k$x1
  y2 <- pop2.k$x2
  t <- y0 + y1 + y2
  idx.rm.pop1 <- which(s < max(s)*call.rate)
  idx.rm.pop2 <- which(t < max(t)*call.rate)
  idx.rm.all <- union(idx.rm.pop1,idx.rm.pop2)
  x0 <- x0[-idx.rm.all]
  x1 <- x1[-idx.rm.all]
  x2 <- x2[-idx.rm.all]
  s  <- s[-idx.rm.all]
  y0 <- y0[-idx.rm.all]
  y1 <- y1[-idx.rm.all]
  y2 <- y2[-idx.rm.all]
  t  <- t[-idx.rm.all]
  # compute SNP_Fst and global Fst estimators in presence of inbreeding   
  e.x <- ((x1 + 2*x2)/(2*s) - (y1 + 2*y2)/(2*t))^2 + x1/(4*s*s) + y1/(4*t*t)
  e.h1 <- (x0*x2 + (x0 + x2)*x1/2 + x1*(x1-1)/4)/(s*(s-1))
  e.h2 <- (y0*y2 + (y0 + y2)*y1/2 + y1*(y1-1)/4)/(t*(t-1))
  N <- e.x - e.h1/s - e.h2/t
  D <- N + e.h1 + e.h2
  Fst.v <- N/D
  names(Fst.v) <- row.names(pop1.k[-idx.rm.all,])
  Fst.o <- Fst.v[order(Fst.v,decreasing=TRUE)]
  F.global <- sum(N)/sum(D)
  se1 <- sd(N)/sqrt(length(N))
  se2 <- sd(D)/sqrt(length(N))
  se.F <- sqrt(se1*se1 + se2*se2)
  F_L95 <- F.global - 1.96*se.F
  F_U95 <- F.global + 1.96*se.F
  Z <- F.global/se.F
  p <- 2*(1 - pnorm(Z))
  if(p < 2e-16){p <- "less than 2e-16"}
  output <- list()
  output[[1]] <- c(F.global,F_L95,F_U95,p)
  names(output[[1]]) <- c("Reich.Fst","L.95%.CI","U.95%.CI","p.val")
  output[[2]] <- data.frame(Fst.o[1:top.number])
  names(output[[2]]) <- c("Reich.Fst")
  return(output)
}

## Generate list of population dataframes from a genind object
DATA_list <- seppop(DATA_genind)
allPopTables <- lapply(DATA_list, makePopVarTable)
# Vectors of which allele (1 first or 2 second) is variant (minor) and common (major)
minorAllele <- as.integer(apply(DATA_genind@tab, 2, function(X) sum(X, na.rm=T))[seq(1,ncol(DATA_genind@tab)-1,2)] > apply(DATA_genind@tab, 2, function(X) sum(X, na.rm=T))[seq(2,ncol(DATA_genind@tab),2)])+1
majorAllele <- as.integer(apply(DATA_genind@tab, 2, function(X) sum(X, na.rm=T))[seq(1,ncol(DATA_genind@tab)-1,2)] < apply(DATA_genind@tab, 2, function(X) sum(X, na.rm=T))[seq(2,ncol(DATA_genind@tab),2)])+1
# Function for generating table
makePopVarTable <- function(gi_obj){
  # x0: 0 copies of variant allele
  x0in <- apply(gi_obj@tab[,seq(1,ncol(gi_obj@tab)-1,2)+(minorAllele-1)], 2, function(X) sum(X==0,na.rm=T))
  # x1: 1 copy of variant allele and 1 copy of common allele (if major is 1, minor will be 1, so can count those equal to 1)
  x1in <- apply(gi_obj@tab[,seq(1,ncol(gi_obj@tab)-1,2)+(minorAllele-1)], 2, function(X) sum(X==1,na.rm=T))
  # x2: 2 copies of variant allele
  x2in <- apply(gi_obj@tab[,seq(1,ncol(gi_obj@tab)-1,2)+(minorAllele-1)], 2, function(X) sum(X==2,na.rm=T))
  # Define dataframe
  popTable <- data.frame(row.names=locNames(gi_obj), x0=x0in, x1=x1in, x2=x2in, A1=as.character(majorAllele), A2=as.character(minorAllele))
  return(popTable)
}

# Function for calculating pairwise population comparisons
pairwise.reich.fst <- function(popTab){
  # m <- list(...)
  n_mod <- length(popTab)
  # All possible pairwise comparisons
  combs <- t(combn(x=names(popTab), m=2))
  # Calculate Reich FST
  # comp_value <- apply(X=combs, MARGIN=1, function(ind) pchisq(2 * (logLik(m[[ind[2]]]) - logLik(m[[ind[1]]])), df = abs(m[[ind[1]]]$df.residual - m[[ind[2]]]$df.residual), lower.tail = FALSE))
  comp_value <- apply(X=combs, MARGIN=1, function(popComp) Reich.Fst(popTab[[popComp[1]]],popTab[[popComp[2]]]))
  FST_val <- sapply(comp_value, "[[", 1)[1,]
  FST_95low <- sapply(comp_value, "[[", 1)[2,]
  FST_95upp <- sapply(comp_value, "[[", 1)[3,]
  FST_pval <- sapply(comp_value, "[[", 1)[4,]
  df_out <- data.frame(combs, FST_val, FST_95low, FST_95upp, FST_pval)
  names(df_out) <- c("pop_1","pop_2","FST_Reich", "FST_95low", "FST_95upp", "FST_pval")
  return(df_out)
}

ReichFST_all <- pairwise.reich.fst(allPopTables)

