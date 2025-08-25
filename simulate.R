set.seed(123)
tp=c()
fp=c()
recall=c()
for (sim.itr in 1:10)
{
  
  num.sites=1000
  num.mutations=10
  num.controls=50
  num.cases=50
  genotype=list(c('A','A'),c('A','T'),c('T','T'))
  num.alleles=length(genotype)
  alleles.controls=matrix(nrow=0,ncol=num.controls)
  alleles.cases=matrix(nrow=0,ncol=num.cases)
  for (i in 1:num.sites)
  {
    basic.probs=gtools::rdirichlet(1,rep(1,num.alleles))
    diff.vec=gtools::rdirichlet(1,rep(1,num.alleles))
    diff.probs=(diff.vec*basic.probs)/sum(diff.vec*basic.probs)
    
    alleles.controls=rbind(alleles.controls,matrix(sample(num.alleles,num.controls,prob = basic.probs,replace=T),nrow=1))
    if (i<=num.mutations)
      alleles.cases=rbind(alleles.cases,matrix(sample(num.alleles,num.cases,prob = diff.probs,replace=T),nrow=1))
    else
      alleles.cases=rbind(alleles.cases,matrix(sample(num.alleles,num.cases,prob = basic.probs,replace=T),nrow=1))
    
  }
  
  rownames(alleles.controls)=paste0('site-',1:nrow(alleles.controls))
  
  res=snp(alleles.controls = alleles.controls,alleles.cases = alleles.cases,num.alleles = num.alleles,
          mcmc.warmup = 20000,mcmc.iter = 200000,allele.level = F,n.cores = 20)
  
  res=res[order(as.numeric(res[,'P'])),]
  
  if (min(as.numeric(res[,'P']))>0.05)
  {
    tp=c(tp,0)
    
    fp=c(fp,0)
    
    recall=c(recall,0)
    
    next
    
  }
  
  thresh.idx=max(which(cumsum(as.numeric(res[,'P']))/(1:nrow(res))<=0.05))
  
  tp=c(tp,sum(res[1:thresh.idx,'Site'] %in% paste0('site-',1:num.mutations))/thresh.idx) #TP
  
  fp=c(fp,sum(!res[1:thresh.idx,'Site'] %in% paste0('site-',1:num.mutations))/thresh.idx) #FP
  
  recall=c(recall,sum(res[1:thresh.idx,'Site'] %in% paste0('site-',1:num.mutations))/num.mutations) #Recall
  
  write.table(res,paste0('ropesnp',sim.itr,'.txt'),sep='\t',col.names = T,row.names = F,quote = F)
  
  print(paste(c('TP:','FP:','Recall'),c(tp[length(tp)],fp[length(fp)],recall[length(recall)])))

}

out.file=cbind(tp,fp,recall)

colnames(out.file)=c('TP','FP','Recall')

write.table(out.file,'sim_ropesnp.txt',sep='\t',row.names = F,col.names = T,quote = F)

