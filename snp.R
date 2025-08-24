library(rstan)
library(coda)
library(bayestestR)
library(parallel)
source('snp_joint.R')
#alleles are always numbered 1..num.alleles
#alleles.controls (alleles.cases) are matrices of dimension num.sites over num.controls(num.cases) , 
#i.e for each site (row) specify the allele, for each patient (column)
#the number of alleles can be more than 4, for example if we consider homo- vs heterozygotes
#The row names of allele.controls must be the names of the genomic sites
snp=function(alleles.controls,alleles.cases,num.alleles,mcmc.warmup=3000,mcmc.iter=4000,allele.level=F,n.cores=1) 
{
  theta=snp.joint(alleles.controls = alleles.controls ,alleles.cases = alleles.cases ,num.alleles = num.alleles,opt.iter = 1000)
  
  functions.string <- "\n  \n  functions{\n  
                  real dirichlet_lpdf_m(vector y,real theta){\n
                  real s=dirichlet_lpdf(y|rep_vector(50.0,num_elements(y)));\n  
                  real l=dirichlet_lpdf(y|rep_vector(1.0,num_elements(y)));\n  
                  if (l>s)\n  \n  return log(theta)+l;\n  
                  return log(1-theta)+s;\n  \n  
                    }\n  \n
                  }\n"
  
  modelString = paste0(functions.string, "data {\n   
                         int<lower=1> Nalleles;
                         int<lower=1> Ncondition1;\n  
                         int<lower=1> Ncondition2;\n  
                         int alleles_controls[Ncondition1];
                         int alleles_cases[Ncondition2];
                         \n\n} \n  parameters {
                         simplex[Nalleles] prob_controls;\n 
                         simplex[Nalleles] alpha;\n 
                         
                         }")
  
  modelString = paste(modelString, "\n  model {\n 
                        \n\n  
                          ")
  
  modelString = paste(modelString, "\n  
                         target+=dirichlet_lpdf_m(alpha,",theta,");  \n
                         prob_controls ~  dirichlet(rep_vector(1.0,Nalleles)); \n
                         for (j in 1:Ncondition1)
                            alleles_controls[j]~categorical(prob_controls);
                         for (j in 1:Ncondition2)
                            alleles_cases[j]~categorical((prob_controls.*alpha)/dot_product(prob_controls,alpha));
                         
                        }")
  
  stan.mod = stan_model(model_code = modelString)
  
  results.tab = do.call(rbind, mclapply(1:nrow(alleles.controls), function(site.number) {
  
   prob.controls=unlist(lapply(1:num.alleles,function(x)(0.001+sum(alleles.controls[site.number,]==x))))
   prob.controls=prob.controls/sum(prob.controls)
    
   initf = function() {
      list(theta = 0.1, alpha=rep(1/num.alleles,num.alleles),
           prob_controls=prob.controls)
      
    }
    res = matrix(ncol = 4, nrow = 0)
    colnames(res) = c("Site", "Allele", "FC", "P")
    dataList = list(Ncondition1 = ncol(alleles.controls),
                    Ncondition2 = ncol(alleles.cases),
                    alleles_controls=alleles.controls[site.number,],
                    alleles_cases=alleles.cases[site.number,],
                    Nalleles=num.alleles
    )
    stanFit = sampling(
      object = stan.mod,
      data = dataList,
      cores = 1,
      init = initf,
      chains = 1,
      refresh = 0,
      iter = mcmc.warmup + mcmc.iter,
      warmup = mcmc.warmup,
      thin = 1
    )
    mcmcCoda = mcmc.list(lapply(1:ncol(stanFit), function(x) {
      mcmc(as.array(stanFit)[, x, ])
    }))
       
    site.name = rownames(alleles.controls)[site.number]
    prob.controls = rep(0, num.alleles)
    alpha = rep(0, num.alleles)
    alpha.p = rep(0, num.alleles)
    for (i in (1:num.alleles)) {
      next.frac.col = which(colnames(mcmcCoda[[1]]) ==
                              paste0("prob_controls[", i, "]"))
      prob.controls[i] = mean(mcmcCoda[[1]][, next.frac.col])
      next.alpha.col = which(colnames(mcmcCoda[[1]]) ==
                               paste0("alpha[", i, "]"))
      alpha[i] = mean(mcmcCoda[[1]][, next.alpha.col])
    }
    if (allele.level) {
      uni.val = sum(prob.controls * alpha)
    }
    else {
      uni.val = 1 / num.alleles
    }
    for (i in (1:num.alleles)) {
      next.alpha.col = which(colnames(mcmcCoda[[1]]) ==
                               paste0("alpha[", i, "]"))
      alpha.p[i] = p_rope(as.numeric(mcmcCoda[[1]][, next.alpha.col]) -
                            uni.val,
                          range = c(-0.2 * uni.val, 0.2 * uni.val))$p_ROPE
    }
     prob.cases = (prob.controls * alpha) /
      sum(prob.controls * alpha)
    fc = prob.cases / prob.controls
      
    if (allele.level) {
      for (i in (1:num.alleles))
        res = rbind(res, c(site.name, i, fc[i], alpha.p[i]))
    }
    else {
      min.p = min(alpha.p)
      res = rbind(res, c(site.name, "mutated", fc[which(alpha.p ==
                                                              min(alpha.p))[1]], min.p))
    }
    return(res)
  }, mc.cores = n.cores))
  
  return(results.tab)
}
