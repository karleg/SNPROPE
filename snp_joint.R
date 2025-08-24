snp.joint=function(alleles.controls,alleles.cases,num.alleles, opt.iter = 1000) 
{
 
  functions.string <- "\n  \n  functions{\n  
                real dirichlet_lpdf_m(vector y,real theta){\n
                real s=dirichlet_lpdf(y|rep_vector(50.0,num_elements(y)));\n  
                real l=dirichlet_lpdf(y|rep_vector(1.0,num_elements(y)));\n  
                if (l>s)\n  \n  return log(theta)+l;\n  
                return log(1-theta)+s;\n  \n  
                  }\n  \n
                }\n"
  
  modelString = paste0(functions.string, "data {\n  int<lower=1> Nsites;\n 
                       int<lower=1> Nalleles;
                       int<lower=1> Ncondition1;\n  
                       int<lower=1> Ncondition2;\n  
                       int alleles_controls[Nsites,Ncondition1];
                       int alleles_cases[Nsites,Ncondition2];
                       \n\n} \n  parameters {
                       simplex[Nalleles] prob_controls[Nsites];\n 
                       simplex[Nalleles] alpha[Nsites];\n 
                       real<lower=0,upper=1> theta;
                    
                       }")

  modelString = paste(modelString, "\n  model {\n\n  \n  \n  
                      \n\n  
                        ")
  
  modelString = paste(modelString, "\n  \n  for (i in 1:Nsites)\n  {\n  
                       target+=dirichlet_lpdf_m(alpha[i],theta);  \n
                       prob_controls[i] ~  dirichlet(rep_vector(1.0,Nalleles)); \n
                       for (j in 1:Ncondition1)
                          alleles_controls[i,j]~categorical(prob_controls[i]);
                       for (j in 1:Ncondition2)
                          alleles_cases[i,j]~categorical((prob_controls[i].*alpha[i])/dot_product(prob_controls[i],alpha[i]));
                       
                      }}")
  
  stan.mod = stan_model(model_code = modelString)
  
  prob.controls=matrix(nrow=nrow(alleles.controls),ncol=num.alleles)

  for (i in 1:nrow(prob.controls))  
  {
    prob.controls[i,]=unlist(lapply(1:num.alleles,function(x)(0.05+sum(alleles.controls[i,]==x))))
    prob.controls[i,]=prob.controls[i,]/sum(prob.controls[i,])
  }
                                                                                         
  init_l = list(theta = 0.1, alpha=matrix(rep(1/num.alleles,nrow(alleles.cases)*num.alleles),ncol=num.alleles),
                prob_controls=prob.controls,ncol=num.alleles)
                
  
  
  dataList = list(Ncondition1 = ncol(alleles.controls),
                  Ncondition2 = ncol(alleles.cases),
                  alleles_controls=alleles.controls,
                  alleles_cases=alleles.cases,
                  Nsites=nrow(alleles.controls),
                  Nalleles=num.alleles
                  )
  
  stanFit = optimizing(object = stan.mod, data = dataList, 
                       init = init_l, iter = opt.iter)
  
  theta = stanFit$par[grepl("theta", names(stanFit$par))]
  return(theta)
}
