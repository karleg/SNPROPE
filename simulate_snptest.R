set.seed(123)
tp=c()
fp=c()
recall=c()
for (sim.itr in 1:10)
{
  
  num.sites=10000
  num.mutations=100
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
  
  sample.file=matrix(nrow=0,ncol=3)
  
  write.table(x =matrix(c('sample_id','sex','phenotype'),nrow=1) ,file = 'sim.sample',sep = '\t',col.names = F,row.names = F,quote = F)

  write.table(x =matrix(c(0,'D','B'),nrow=1) ,file = 'sim.sample',sep = '\t',append = T,col.names = F,row.names = F,quote = F)
  
  for (i in 1:num.controls)
     sample.file=rbind(sample.file,c(i,'M','control'))
  
  for (i in (num.controls+1):(num.controls+num.cases))
    sample.file=rbind(sample.file,c(i,'M','case'))
  
  write.table(sample.file,'sim.sample',sep='\t',quote = F,row.names = F,col.names = F,append = T)
  
  map.file=matrix(nrow=0,ncol=4)
  
  
  for (i in 1:num.sites)
  {
    chr=sample(1:22,1)
    bp=sample(1:10^6,1)
    map.file=rbind(map.file,c(chr,paste0('snp',i),0,bp))
  }
  write.table(map.file,'sim.map',sep=' ',quote = F,row.names = F,col.names = F)
  
  ped.file=matrix(nrow=0,ncol=6+num.sites)
  for (i in 1:(num.controls+num.cases))
  {
    if (i<=num.controls)
      ped.file=rbind(ped.file,c(1,i,0,0,1,1,unlist(lapply(1:num.sites,function(j)paste(genotype[[alleles.controls[j,i]]],collapse = ' ')))))
    else
      ped.file=rbind(ped.file,c(1,i,0,0,1,2,unlist(lapply(1:num.sites,function(j)paste(genotype[[alleles.cases[j,i-num.controls]]],collapse = ' ')))))
  }
  
  write.table(ped.file,'sim.ped',sep='  ',quote = F,row.names = F,col.names = F)

  system('sh snptest.sh')

  snptest.out=read.table('snptest.out',header = T)
  
  padj=p.adjust(snptest.out$frequentist_add_lrt_pvalue,'BH')
  
  if (min(padj,na.rm = T)>0.05)
  {
    tp=c(tp,0)
    
    fp=c(fp,0)
    
    recall=c(recall,0)
    
    next
    
  }
  
  tp=c(tp,sum(unique(snptest.out$rsid[padj<=0.05 & !is.na(padj)]) %in% paste0('snp',1:num.mutations))/sum(padj<=0.05,na.rm=T)) #TP
  
  fp=c(fp,sum(!unique(snptest.out$rsid[padj<=0.05 & !is.na(padj)]) %in% paste0('snp',1:num.mutations))/sum(padj<=0.05,na.rm=T)) #FP
  
  recall=c(recall,sum(unique(snptest.out$rsid[padj<=0.05 & !is.na(padj)]) %in% paste0('snp',1:num.mutations))/num.mutations) #Recall
}


out.file=cbind(tp,fp,recall)

colnames(out.file)=c('TP','FP','Recall')

write.table(out.file,'sim_snptest.txt',sep='\t',row.names = F,col.names = T,quote = F)
