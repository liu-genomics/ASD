#################################################################
# TADA-Denovo: analysis of de novo data 
#################################################################

# Bayes factor of de novo counts of a gene 
# x: the de novo count
# N: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events)
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
bayes.factor.denovo <- function(x, N, mu, gamma.mean, beta) {
  marg.lik0 <- dpois(x, 2*N*mu)
  marg.lik1 <- dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu))
  BF <- marg.lik1/marg.lik0
  
  return (BF)
}

# Genome-wide application of TADA for denovo data
# Input: counts, N, mu, gamma.mean, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: m x K matrix, where m and K have the same meanings. mu[i,j] is the mutation rate of the j-th category of the i-th gene
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA.denovo <- function(counts, N, mu, gamma.mean, beta) {
  m <- dim(counts)[1]
  K <- dim(counts)[2]
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  # Compute BFs for each of the j-th category of mutations
  for (i in 1:m) {
    for (j in 1:K)  BF[i,j] <- bayes.factor.denovo(counts[i,j], N, mu[i,j], gamma.mean[j], beta[j])
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene: de novo only
# mu.gene: the mutation rates of a gene (K-dim. vector)
# N: the sample size
# l: the number of permutations
# gamma.mean, beta: RR of de novo mutations (vectors)
# Output: BF - l BFs from permutation; sample - permutate data
permute.gene.denovo <- function(mu.gene, N, l, gamma.mean, beta, dn.max=5) {
  K <- length(mu.gene)
  BF.gene <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  
  # permutation of l times
  count.gene <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N, mu.gene[j], gamma.mean[j], beta[j])
    }
    
    # permutation
    count.gene[,j] <- rpois(l, 2*N*mu.gene[j])
    for (i in 1:l) {
      x <- count.gene[i,j]
      cond.range.dn <- (x <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.gene[i,j] <- table.dn[x+1]
      } else {
        BF.gene[i,j] <- bayes.factor.denovo(x, N, mu.gene[j], gamma.mean[j], beta[j]) 
      }
    }
  }
  
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample=count.gene))
}

# Genome-wide application of TADA for denovo data: the difference with TADA.denovo is the report of p-values. 
# Input: counts, N, mu, mu.frac, gamma.mena, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: m x K matrix, where m and K have the same meanings. mu[i,j] is the mutation rate of the j-th category of the i-th gene
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# l: the number of permutations to obtain the null distribution
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs. 
TADAp.denovo <- function(counts, N, mu, gamma.mean, beta, l=100, dn.max=5) {
  m <- dim(counts)[1]
  K <- dim(counts)[2]
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs
  rs <- TADA.denovo(counts, N, mu, gamma.mean, beta)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene.denovo(mu[i,], N, l, gamma.mean, beta, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# MOM estimation of hyperprior parameters from de novo data
#################################################################

# Prob. of having d or more de novo mutations under H1 
# Use simulation, but could also use analytic form 
multihit.prob <- function(N, mu, gamma.mean, beta, d=2, S=100) {
  p <- numeric(S)
  gamma <- rgamma(S, gamma.mean*beta, rate=beta)
  for (i in 1:S) {
    p[i] <- 1 - ppois(d-1, 2*N*mu*gamma[i])
  }
  return (mean(p))
}

# Estimate the number of multihit genes in a genome. 
# d: the parameter of the multiplicity test. 
# Returns: M0 - the number of multihit genes from non-risk genes; M1 - the number from risk genes. 
count.multihit <- function(N, mu, pi, gamma.mean, beta, d=c(2,3), S=2) {
  m <- length(mu)
  M0 <- numeric(length(d))
  M1 <- numeric(length(d))
  
  # M1: the number of causal genes having d or more de novo mutations
  p.alt <- array(0, dim=c(m, length(d)))  # p1[i,j] = P(X_i >= d_j|H1)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.alt[i,j] <- multihit.prob(N, mu[i], gamma.mean, beta, d=d[j], S=S)
    }
  }
  for (j in 1:length(d)) { 
    M1[j] <- m * pi  * mean(p.alt[,j]) 
  }
  
  # M0: the number of non-causal genes having d or more de novo mutations
  p.null <- array(0, dim=c(m, length(d)))  # p0[i,j] = P(X_i >= d_j|H0)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.null[i,j] <- 1 - ppois(d[j] - 1, 2*N*mu[i])
    }
  }
  for (j in 1:length(d)) { 
    M0[j] <- m * (1-pi) * mean(p.null[,j]) 
  }
  
  result <- data.frame(d=d, M0=M0, M1=M1)
  return (result)
}

# Estimating relative risk and the number of multiple hits from de novo data
# Input: sample size (N), mutation rates of all genes (mu), observed number of de novo events (C), beta (parameter of the prior distribution of gamma), k (number of disease genes)
# Output: the average relative risk (gamma.mean), the expected number of multi-hit genes (M)
denovo.MOM <- function(N, mu, C, beta, k) {
  m <- length(mu) # number of genes
  
  # enrichment of de novo events
  nu <- C / (2 * N * sum(mu))
  
  # MOM estimator of gamma.mean
  gamma.mean <- (nu-1)*m/k +1
  
  # expected M (choose d = 2)
  rs <- count.multihit(N, mu, k/m, gamma.mean, beta, d=2)    
  M <- sum(rs$M1) + sum(rs$M0)
  
  return (list(gamma.mean=gamma.mean, M=M))
}

#################################################################
# Useful functions
#################################################################

# Bayesian FDR control (PMID:19822692, Section2.3)
# BF: a sorted vector of BFs (in decreasing order)
# pi0: the prior probability that the null model is true
# alpha: the FDR target
# Return: the q-value of each BF, and the number of findings with q below alpha. 
Bayesian.FDR <- function(BF, pi0, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}