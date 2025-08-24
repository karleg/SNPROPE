source('snp.R')
library(GenomicRanges)

cellsnp.files=list.files('cellSNP_output.txt/',full.names = T)

snp.tab=read.table(cellsnp.files[1],sep='\t')

colnames(snp.tab)[10]=gsub('.txt.gz','',unlist(strsplit(cellsnp.files[1],'_'))[4])

colnames(snp.tab)[c(1,2,4,5)]=c('chr','bp','allele','alt')

alleles=unlist(lapply(snp.tab[,10],function(x)unlist(strsplit(x,':'))[1]))

snp.tab[,10]=alleles

snp.tab=snp.tab[,(colnames(snp.tab) %in% c('chr','bp')) | grepl('SRR',colnames(snp.tab))]

col.names=colnames(snp.tab)

snp.tab=do.call(rbind,lapply(split(snp.tab,paste(snp.tab$chr,snp.tab$bp)),function(m){
  
  v=unlist(c(m[1,1:2]))
  if (sum(m[,3]=='0/0',na.rm = T)>0 && sum(m[,3]!='0/0',na.rm = T)==0)
     return(c(v,2)) #reference allele
  if (sum(m[,3]!='0/0',na.rm = T)>0)
     return(c(v,1)) #mutation
  return(c(v,3)) #no calls
}))

colnames(snp.tab)=col.names

metadata=read.csv('SraRunTable_Hooks.csv')

control.ids=metadata$Run[metadata$cell_type=='normal']

case.ids=metadata$Run[metadata$cell_type=='tumour']

for (file in cellsnp.files[2:length(cellsnp.files)])
{
  next.tab=read.table(file,sep='\t')
  
  colnames(next.tab)[10]=gsub('.txt.gz','',unlist(strsplit(file,'_'))[4])
  
  colnames(next.tab)[c(1,2,4,5)]=c('chr','bp','allele','alt')
  
  alleles=unlist(lapply(next.tab[,10],function(x)unlist(strsplit(x,':'))[1]))
  
  next.tab[,10]=alleles
  
  next.tab=next.tab[,(colnames(next.tab) %in% c('chr','bp')) | grepl('SRR',colnames(next.tab))]

  col.names=colnames(next.tab)
  
  next.tab=do.call(rbind,lapply(split(next.tab,paste(next.tab$chr,next.tab$bp)),function(m){
    
    v=unlist(c(m[1,1:2]))
    if (sum(m[,3]=='0/0',na.rm = T)>0 && sum(m[,3]!='0/0',na.rm = T)==0)
      return(c(v,2)) #reference allele
    if (sum(m[,3]!='0/0',na.rm = T)>0)
      return(c(v,1)) #mutation
    return(c(v,3)) #no calls
  }))
  
  colnames(next.tab)=col.names
  
  snp.tab=merge(snp.tab,next.tab,by = c('chr','bp'),all=T)
  
  snp.tab=snp.tab[,(colnames(snp.tab) %in% c('chr','bp')) | grepl('SRR',colnames(snp.tab))]
  
  snp.tab=as.data.frame(snp.tab)
  
  print(nrow(snp.tab))
}

rm(next.tab)

pcg=read.table('protein_coding_genes.tsv',header = T,sep='\t')

pcg=pcg[!is.na(pcg$Chromosomes),]

pcg=pcg[!is.na(pcg$Annotation.Genomic.Range.Start),]

pcg=pcg[!is.na(pcg$Annotation.Genomic.Range.Stop),]

pcg=pcg[!duplicated(pcg$NCBI.GeneID),]

snp.coords=GRanges(paste0(gsub('chr','',snp.tab$chr),':',snp.tab$bp,'-',snp.tab$bp))

gene.coords=GRanges(paste0(pcg$Chromosomes,':',pcg$Annotation.Genomic.Range.Start,'-',pcg$Annotation.Genomic.Range.Stop),gene_id=pcg$NCBI.GeneID)

snp.gene.matches=findOverlaps(snp.coords,gene.coords,select='first')

snp.tab=snp.tab[!is.na(snp.gene.matches),]

snp.coords=snp.coords[!is.na(snp.gene.matches)]

snp.gene.matches=snp.gene.matches[!is.na(snp.gene.matches)]

snp.tab=do.call(rbind,lapply(split(snp.tab,snp.gene.matches),function(m){

  unlist(c(m[1,1:2],apply(m[,-c(1,2)],2,function(v)
    {
    if (sum(v==1,na.rm = T)>0)
      return(1)
    if (sum(v==2,na.rm = T)>0)
      return(2)
    return(3)
    })))
}))

metadata=metadata[metadata$Run %in% colnames(snp.tab),]

metadata=metadata[match(colnames(snp.tab)[-c(1:2)],metadata$Run),]

sum(metadata$Run!=colnames(snp.tab)[-c(1:2)])

controls.idx=which(metadata$cell_type=='normal')

cases.idx=which(metadata$cell_type=='tumour')

snp.tab=as.data.frame(snp.tab)

snp.tab=snp.tab[apply(snp.tab[,-c(1,2)],1,function(v)sum(v==3)<=15),] #remove if more than 25% of samples are missing data

alleles.cases=snp.tab[,2+cases.idx]

alleles.controls=snp.tab[,2+controls.idx]

alleles.cases=apply(alleles.cases,2,as.integer)

alleles.controls=apply(alleles.controls,2,as.integer)

rownames(alleles.controls)=rownames(alleles.cases)=paste(snp.tab$chr,snp.tab$bp,sep='-')

res=snp(alleles.controls = alleles.controls,alleles.cases = alleles.cases,num.alleles = 3,mcmc.warmup = 3000,mcmc.iter = 10000,allele.level = T,n.cores = 7)

res=res[order(as.numeric(res[,'P'])),]

gene.ids=c()

for (i in (1:nrow(res)))
{
  chr=gsub('chr','',strsplit(res[i,1],split='-')[[1]][1])
  bp=strsplit(res[i,1],split='-')[[1]][2]
  gene.ids=c(gene.ids,gene.coords[snp.gene.matches[which(snp.coords@seqnames==chr & snp.coords@ranges==paste(bp,bp,sep='-'))]]@elementMetadata[1,1])
}
res=cbind(res,gene.ids)
colnames(res)[5]='GeneId'


write.table(res,'snp_results_Hooks.txt',sep='\t',quote = F,row.names = F,col.names = T)

thresh.idx=max(which(cumsum(as.numeric(res[,'P']))/(1:nrow(res))<=0.05))

selected.genes=unique(res[1:thresh.idx,'GeneId'][res[1:thresh.idx,'Allele']!=3])

#remove alleles if a significant difference between the proportions of missing samples in cases and controls was detected
selected.genes=selected.genes[unlist(lapply(selected.genes,function(x)! ('3' %in% res[res[,'GeneId']==x & (1:nrow(res))<thresh.idx,'Allele'])))]

fc=c()

for (gene in selected.genes)
{
  tab.ctr=table(alleles.controls[rownames(alleles.controls)==res[res[,'GeneId']==gene,'Site']])

  tab.case=table(alleles.cases[rownames(alleles.controls)==res[res[,'GeneId']==gene,'Site']])
  
  fc=c(fc,tab.case['1']/tab.ctr['1'])
   
}

library(ggpubr)

df=data.frame(GeneId=selected.genes,FC=fc)

ggbarplot(df,x='GeneId',y='FC',xlab='Entrez Gene ID',ylab='Fold-change',fill='blue',col='blue',palette = 'npg')

