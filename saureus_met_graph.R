saogenes=keggList("sao")

saogenesid=names(saogenes)
saorn=vector()
for (i in 1:length(saogenesid)){
  geneform=keggGet(saogenesid[i])
  orthid=names(geneform[[1]]$ORTHOLOGY)
  if (length(orthid)>0){
    for (j in 1:length(orthid)){
      orthform=keggGet(orthid[j])
      dblinks=orthform[[1]]$DBLINKS
      rnlinks=grep("RN:",dblinks,fixed=T)
      if (length(rnlinks)>0){
        saorn=c(saorn,dblinks[rnlinks])
      }      
    }
    
  }
   
  print(i)
}


save(saorn,file="saorn.RData")
saorn=load("saorn.RData")

splitrn=strsplit(saorn,split=":")
saorn2=vector()
for (i in 1:length(saorn)){
  allrn=setdiff(unlist(strsplit(splitrn[[i]][2],split=" ")),"")
  saorn2=c(saorn2,allrn)
}

saorn3=unique(saorn2)

all_met=vector()
rncompound=list()
rncompound$react=list()
rncompound$prod=list()
for (i in 1:length(saorn3)){
  eq=keggGet(saorn3[i])[[1]]$EQUATION
  eqsplit=strsplit(eq,"<=>")
  eqsplit=unlist(eqsplit)
  r=unlist(strsplit(eqsplit[1],"+",fixed=T))
  p=unlist(strsplit(eqsplit[2],"+",fixed=T))
  r=gsub(" ", "", r, fixed = TRUE)
  p=gsub(" ", "", p, fixed = TRUE)
  for (j in 1:length(r)){
    r[j]=str_extract(r[j],"C[0-9]{5}")
  }
  for (j in 1:length(p)){
    p[j]=str_extract(p[j],"C[0-9]{5}")
  }
  #Smat[is.element(unimet2,r),i]=-1
  #Smat[is.element(unimet2,p),i]=1
  rncompound$react[[i]]=r
  rncompound$prod[[i]]=p
  all_met=union(all_met,c(r,p))
  print(i)
}


saudirtab=data.frame(from="mm",to="mm",reaction="rr", stringsAsFactors = F)
l=1
for (i in 1:length(saorn3)){
  r=rncompound$react[[i]]
  p=rncompound$prod[[i]]
  if (length(r)*length(p)>0){
    for (f in 1:length(r)){
      for (t in 1:length(p)){
        saudirtab[l,1]=r[f]
        saudirtab[l,2]=p[t]
        saudirtab[l,3]=saorn3[i]
        l=l+1
      }
    }
    
  }
  print(i)
}

save(saudirtab,file="saudirtab2.RData")
load(file="saudirtab2.RData")

all_met=union(saudirtab$from,saudirtab$to)
all_met=sort(all_met)

metname=vector()
metcount=vector()
for (i in 1:length(all_met)){
  print(i)
  metname[i]=keggGet(all_met[i])[[1]]$NAME[1]
  rnpos=union(which(saudirtab$from==all_met[i]),which(saudirtab$to==all_met[i]))
  metcount[i]=length(rnpos)
}

mettab=data.frame(id=all_met,name=metname,count=metcount,path=0,stringsAsFactors = F)
mettab$path[mettab$count<134]=1
mettab[14,4]=0
mettab[7,4]=0
mettab[18,4]=0
mettab[85,4]=0
mettab[28,4]=0
mettab[393,4]=0
mettab[26,4]=0
mettab[160,4]=0
mettab[49,4]=0
mettab[14,4]=0
mettab[88,4]=0
mettab[40,4]=0
mettab[54,4]=0
mettab[14,4]=0
mettab[839,4]=0
mettab[191,4]=0
mettab[210,4]=0
mettab[289,4]=0
mettab[594,4]=0
mettab[848,4]=0
mettab[849,4]=0
mettab[210,4]=0
mettab[16,4]=0
mettab[63,4]=0
mettab[15,4]=0
mettab[193,4]=0
mettab[210,4]=0
mettab[1014,4]=0
mettab[210,4]=0
mettab[1016,4]=0
mettab[25,4]=0
mettab[210,4]=0
mettab[80,4]=0
mettab[210,4]=0
mettab[211,4]=0
mettab[210,4]=0
mettab[413,4]=0
mettab[210,4]=0
mettab[32,4]=0
mettab[30,4]=0
mettab[120,4]=0
mettab[52,4]=0
mettab[444,4]=0
mettab[1015,4]=0
mettab[210,4]=0
mettab[1017,4]=0
mettab[904,4]=0
mettab[108,4]=0
mettab[74,4]=0
mettab[139,4]=0
mettab[150,4]=0
mettab[192,4]=0
mettab[220,4]=0
mettab[1021,4]=0
mettab[505,4]=0
mettab[655,4]=0
mettab[1020,4]=0
mettab[68,4]=0
mettab[109,4]=0
mettab[210,4]=0
mettab[224,4]=0
mettab[210,4]=0
mettab[257,4]=0
mettab[319,4]=0
mettab[412,4]=0
mettab[1019,4]=0
mettab[93,4]=0
mettab[170,4]=0
mettab[427,4]=0
mettab[483,4]=0
mettab[509,4]=0
mettab[708,4]=0
mettab[609,4]=0
mettab[712,4]=0
mettab[718,4]=0
mettab[210,4]=0

write.table(mettab,file="mettab.txt",row.names=F,quote=F,sep="|")
mettab=read.delim(file="mettab_processed.txt",header=T,stringsAsFactors = F)

validmet=mettab$id[mettab$path==1]

saudirtab_v=saudirtab
saudirtab_v=saudirtab_v[is.element(saudirtab_v$from,validmet),]
saudirtab_v=saudirtab_v[is.element(saudirtab_v$to,validmet),]


gsao3=graph_from_data_frame(saudirtab_v[,c(1,2)],directed=F)
gsao3comp=components(gsao3)
gsao3comp$csize

datamet=c("C00058","C00221", "C00267", "C00031", "C00049", "C00073", "C00033", "C00041", "C00186", "C00256", "C00810", "C01769", "C01083", "C00037", "C00022", "C00469", "C05345", "C05378", "C00118", "C00236", "C00197", "C00631", "C00074", "C00024")


vlist3=vertex_attr(gsao3)
mindex3=which(is.element(vlist3[[1]],datamet))
gsao3comp$membership[mindex3]

dmat3=distances(gsao3,v=mindex3,to=mindex3)

metsubnet3=joincompt(mindex3,gsao3)
addednodes3=sort(vlist3[[1]][metsubnet3$uni])

corenet=as_data_frame(simplify(metsubnet3$subg),what="Edges")

corenode=data.frame(id=unique(c(corenet$from,corenet$to)),label="",stringsAsFactors=F)

for (i in 1:length(corenode$id)){
  prelabel=mettab$name[mettab$id==corenode$id[i]]
  #prelabel=strsplit(prelabel,";")
  corenode$label[i]=prelabel[[1]]
}


library(visNetwork)
visNetwork(corenode, corenet) %>% 
  visLayout(randomSeed = 42) 



metdiff=read.table("metdiff2.txt",header=T,stringsAsFactors = F)

metdiff2=metdiff[,-c(13,15)]
inputmet=c("C00058","C00031","C00049","C00073","C00033","C00041","C00186","C00810","C01083","C00037","C01104","C00022","C00543","C00469")
inputmet2=inputmet[-c(11,13)]


gsao4=simplify(gsao3)
vlist4=vertex_attr(gsao4)
inputindex=vector()
for (i in 1:12){
  
  inputindex[i]=which(vlist4[[1]]==inputmet2[i])  
}



vlist4[[1]][inputindex]
inputmet
sourcerwwrmat=data.frame(metabolite=vlist4[[1]])
sinkrwwrmat=data.frame(metabolite=vlist4[[1]])
n=gorder(gsao4)
for (i in 1:9){
  sourcemet=which(metdiff2[i,3:14]<0)
  sinkmet=which(metdiff2[i,3:14]>0)
  sourceid=inputindex[sourcemet]
  sinkid=inputindex[sinkmet]
  sourcew=metdiff2[i,sourcemet+2]
  sinkw=metdiff2[i,sinkmet+2]
  vecsource=rep(0,n)
  for (j in 1:length(sourcemet)){
    vecsource[sourceid[j]]=-sourcew[j]  
  }
  vecsink=rep(0,n)
  for (j in 1:length(sinkmet)){
    vecsink[sinkid[j]]=sinkw[j]  
  }
  sourcerwwrmat[,i+1]=page_rank(gsao4,directed=F,damping=0.85,personalized=vecsource)$vector#rwwr(gsao4,sourceid,-sourcew,0.85)
  sinkrwwrmat[,i+1]=page_rank(gsao4,directed=F,damping=0.85,personalized=vecsink)$vector#rwwr(gsao4,sinkid,sinkw,0.85)
  print(i)
}


s2smat=sourcerwwrmat
s2smat[,2:10]=(sourcerwwrmat[,2:10]*sinkrwwrmat[,2:10])#/(2*abs(sourcerwwrmat[,2:13]-sinkrwwrmat[,2:13]))

reactrwwr=saudirtab_v
reactrwwr$def=""
reactrwwr$t0wt=0
reactrwwr$t0a=0
reactrwwr$t0b=0
reactrwwr$t03wt=0
reactrwwr$t03a=0
reactrwwr$t03b=0
reactrwwr$t36wt=0
reactrwwr$t36a=0
reactrwwr$t36b=0


for (i in 1:length(saudirtab_v$reaction)){
  rowid1=which(sourcerwwrmat$metabolite==saudirtab_v[i,1])
  rowid2=which(sourcerwwrmat$metabolite==saudirtab_v[i,2])
  prodvec=sourcerwwrmat[rowid1,2:10]*sinkrwwrmat[rowid2,2:10]+sourcerwwrmat[rowid2,2:10]*sinkrwwrmat[rowid1,2:10]
  sourcesign=sourcerwwrmat[rowid1,2:10]-sourcerwwrmat[rowid2,2:10]
  sinksign=sinkrwwrmat[rowid1,2:10]-sinkrwwrmat[rowid2,2:10]
  coherentsign=1*(sourcesign*sinksign<0)
  reactdir=sourcesign*coherentsign
  reactrwwr$def[i]=keggGet(saudirtab_v$reaction[i])[[1]]$DEFINITION
  reactrwwr[i,5:13]=prodvec*reactdir
  print(i)
}

reactid=unique(reactrwwr$reaction)

reactrwwru=reactrwwr[1:length(reactid),]

for (i in 1:length(reactid)){
  lines=reactrwwr[reactrwwr$reaction==reactid[i],]
  reactrwwru[i,1:4]=lines[1,1:4]
  reactrwwru[i,5:13]=colSums(lines[,5:13])
}

a=abs(reactrwwru[,5:13])
a[a==0]=1
minrw=min(a)

signa=1*(reactrwwru[,5:13]>=0)-1*(reactrwwru[,5:13]<0)

a=abs(reactrwwru[,5:13])
a[a==0]=minrw*0.5
loga=log2(a)
loga=loga-min(loga)
loga=loga*signa




#reactrwwru[reactrwwru==0]=minrw*0.5
reactrwwru$lrt0a=(-1*(reactrwwru$t0wt>reactrwwru$t0a)+1*(reactrwwru$t0wt<reactrwwru$t0a))*abs(loga$t0wt-loga$t0a)
reactrwwru$lrt03a=(-1*(reactrwwru$t03wt>reactrwwru$t03a)+1*(reactrwwru$t03wt<reactrwwru$t03a))*abs(loga$t03wt-loga$t03a)
reactrwwru$lrt36a=(-1*(reactrwwru$t36wt>reactrwwru$t36a)+1*(reactrwwru$t36wt<reactrwwru$t36a))*abs(loga$t36wt-loga$t36a)
reactrwwru$lrt0b=(-1*(reactrwwru$t0wt>reactrwwru$t0b)+1*(reactrwwru$t0wt<reactrwwru$t0b))*abs(loga$t0wt-loga$t0b)
reactrwwru$lrt03b=(-1*(reactrwwru$t03wt>reactrwwru$t03b)+1*(reactrwwru$t03wt<reactrwwru$t03b))*abs(loga$t03wt-loga$t03b)
reactrwwru$lrt36b=(-1*(reactrwwru$t36wt>reactrwwru$t36b)+1*(reactrwwru$t36wt<reactrwwru$t36b))*abs(loga$t36wt-loga$t36b)

reactrwwru$lpt0a=log2(abs(reactrwwru$t0a*reactrwwru$t0wt))
reactrwwru$lpt03a=log2(abs(reactrwwru$t03a*reactrwwru$t03wt))
reactrwwru$lpt36a=log2(abs(reactrwwru$t36a*reactrwwru$t36wt))
reactrwwru$lpt0b=log2(abs(reactrwwru$t0b*reactrwwru$t0wt))
reactrwwru$lpt03b=log2(abs(reactrwwru$t03b*reactrwwru$t03wt))
reactrwwru$lpt36b=log2(abs(reactrwwru$t36b*reactrwwru$t36wt))

#reactrwwru[reactrwwru==minrw*0.5]=0

tokeep=(rowSums(reactrwwru[,14:19])!=0)

tokeep_a=(rowSums(reactrwwru[,14:16])!=0)

tokeep_b=(rowSums(reactrwwru[,17:19])!=0)

tokeep_a_b=(tokeep_a & tokeep_b)

tokeep_just_a=(tokeep_a & !tokeep_b)

tokeep_just_b=(!tokeep_a & tokeep_b)

reactrwwrs=reactrwwru[tokeep,]

react2path=list()
for (i in 1:nrow(reactrwwru)){
  rform=keggGet(reactrwwru$reaction[i])
  react2path[[i]]=rform[[1]]$PATHWAY
}


pathcountall=as.data.frame(table(unlist(react2path)))

enrichpath=function(react2path,subset){
  pathcountall=as.data.frame(table(unlist(react2path)),stringsAsFactors = F)
  pathsub=react2path[subset]
  pathcountall$subset=pathcount(pathcountall,pathsub)
  pathcountall$fold=rep(0,nrow(pathcountall))
  pathcountall$zscore=rep(0,nrow(pathcountall))
  pathcountall$p=rep(0,nrow(pathcountall))
  
  
  rpathcount=matrix(nrow=nrow(pathcountall),ncol=2500)
  for (i in 1:2500){
    rsubset=sample(length(react2path),length(subset),replace=F)
    rpathsub=react2path[rsubset]
    rpathcount[,i]=pathcount(pathcountall,rpathsub)
  }
  
  for (i in 1:nrow(pathcountall)){
    rmean=mean(rpathcount[i,])
    rsd=sd(rpathcount[i,])
    pathcountall$fold[i]=(pathcountall$subset[i]/length(subset))/(pathcountall$Freq[i]/length(react2path))
    pathcountall$zscore[i]=(pathcountall$subset[i]-rmean)/rsd
    pathcountall$p[i]=sum(rpathcount[i,]>=pathcountall$subset[i])/2500
  }
  pathcountall
  
}

pathcount=function(pathcountall,subsetlist){
  pathvec=unlist(subsetlist)
  countvec=vector()
  for (i in 1:nrow(pathcountall)){
    countvec[i]=sum(pathvec==pathcountall[i,1])
  }
  countvec
}

pcount_a=enrichpath(react2path,which(tokeep_a))
pcount_b=enrichpath(react2path,which(tokeep_b))
pcount_ab=enrichpath(react2path,which(tokeep_a_b))
pcount_just_a=enrichpath(react2path,which(tokeep_just_a))

write.table(pcount_a,file="pcount_a.txt",row.names = F)
write.table(pcount_b,file="pcount_b.txt",row.names = F)


psig_a=pcount_a$Var1[pcount_a$zscore>2 & pcount_a$p<0.05]
psig_b=pcount_b$Var1[pcount_b$zscore>2 & pcount_b$p<0.05]
psig_a_or_b=union(psig_a,psig_b)

pcount_a_or_b=pcount_a
names(pcount_a_or_b)=c("Pathwhay","# reactions", "# changes in -ndh2a", "fold enrichment a", "z-score a", "p a")
pcount_a_or_b$reactions_b=pcount_b$subset
pcount_a_or_b$fold_b=pcount_b$fold
pcount_a_or_b$zscore_b=pcount_b$zscore
pcount_a_or_b$p_b=pcount_b$p
pcount_a_or_b$selected=1*is.element(pcount_a$Var1,psig_a_or_b)
write.table(pcount_a_or_b,file="path_enrich.txt",row.names = F)

pathmatu=matrix(nrow=length(react2path),ncol=length(psig_a_or_b))
for (i in 1:length(psig_a_or_b)){
  for (j in 1:length(react2path)){
    if (is.element(c(psig_a_or_b[i]),react2path[[j]])){
      pathmatu[j,i]=1
    } else {
      pathmatu[j,i]=0
    }
  }
}
pathmatu=as.data.frame(pathmatu,stringsAsFactors=F)
names(pathmatu)=psig_a_or_b
pathmatu=pathmatu[,-c(5,6,10)]




pathmats=pathmatu[tokeep,]

nadvec=vector()
nadpvec=vector()
for (i in 1:nrow(reactrwwru)){
  rlines=which(saudirtab$reaction==reactrwwru$reaction[i])
  mets=union(saudirtab$from[rlines],saudirtab$to[rlines])
  if (is.element(c("C00003"),mets)){
    nadvec[i]=1
  } else {
    nadvec[i]=0
  }
  if (is.element(c("C00005"),mets)){
    nadpvec[i]=1
  } else {
    nadpvec[i]=0
  }
}

row.names(reactrwwrs)=reactrwwrs$def

pdf(file = "sao_path1.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Alanine, aspartate and glutamate metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(7, 40),cexRow = 1.2, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Alanine, aspartate and glutamate metabolism",cex=1.2)
dev.off()

pdf(file = "sao_path2.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Arginine biosynthesis`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(7, 30),cexRow = 1.2, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Arginine biosynthesis",cex=1.2)
dev.off()

pdf(file = "sao_path3.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Biosynthesis of amino acids`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(5, 50),cexRow = 0.7, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Biosynthesis of amino acids",cex=1.2)
dev.off()

pdf(file = "sao_path4.pdf", width = 20,  height = 16) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Biosynthesis of secondary metabolites`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(5, 60),cexRow = 0.7, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Biosynthesis of secondary metabolites",cex=1.2)
dev.off()

pdf(file = "sao_path5.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Cysteine and methionine metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 40),cexRow = 1.2, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Cysteine and methionine metabolism",cex=1.2)
dev.off()

pdf(file = "sao_path6.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Glycerolipid metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 40),cexRow = 1.3, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Glycerolipid metabolism",cex=1.2)
dev.off()

pdf(file = "sao_path7.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Glycine, serine and threonine metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 30),cexRow = 1.2, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Glycine, serine and threonine metabolism",cex=1.2)
dev.off()

pdf(file = "sao_path8.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`One carbon pool by folate`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 50),cexRow = 1, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="One carbon pool by folate",cex=1.2)
dev.off()

pdf(file = "sao_path9.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Pentose phosphate pathway`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 40),cexRow = 1, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Pentose phosphate pathway",cex=1.2)
dev.off()

pdf(file = "sao_path10.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Phenylalanine, tyrosine and tryptophan biosynthesis`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 40),cexRow = 1, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Phenylalanine, tyrosine and tryptophan biosynthesis",cex=1.2)
dev.off()

pdf(file = "sao_path11.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Purine metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 40),cexRow = 0.8, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Purine metabolism",cex=1.2)
dev.off()

pdf(file = "sao_path12.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Amino sugar and nucleotide sugar metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 40),cexRow = 1, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Amino sugar and nucleotide sugar metabolism",cex=1.2)
dev.off()

pdf(file = "sao_path13.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Glycolysis / Gluconeogenesis`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 50),cexRow = 1, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Glycolysis / Gluconeogenesis",cex=1.2)
dev.off()

pdf(file = "sao_path14.pdf", width = 18,  height = 12) # The height of the plot in inches
heatmap3(as.matrix(reactrwwrs[pathmats$`Methane metabolism`==1,14:19]),Colv=NA,balanceColor=T,scale="none",margins=c(6, 30),cexRow = 1, labCol=c("~0h a", "0-3h a", "3-6h a","~0h b", "0-3h b", "3-6h b"),xlab="Methane metabolism",cex=1.2)
dev.off()

logabs=function(mat){
  signmat=mat/abs(mat)
  absmat=abs(mat)
  a=log10(absmat)
  a[a==-Inf]=0
  a[a==0]=min(a)-0.5
  a=a-min(a)+0.5
  a=a/max(a)
  a=a*10
  a=a*signmat
  a
}



write.table(reactrwwru,file="reacttab2.txt",row.names=F,sep="\t",quote=F)

reactrwwru=read.delim("reacttab2.txt",header=T,stringsAsFactors = F)

write.table(reactrwwrs,file="reacttabs2.txt",row.names=F,sep="\t",quote=F)

reactrwwrs=read.delim("reacttabs2.txt",header=T,stringsAsFactors = F)

rmax=vector()
for (i in 1:nrow(reactrwwrs)){
  rmax[i]=max(log10(abs(reactrwwrs[i,14:19])))
}
summary(rmax)
topQ=(rmax>1.9121)
top10=(rmax>-5.35)
top5=(rmax>-4.6)


pcount_Q=enrichpath(react2path,which(topQ))

reactQ=reactrwwrs[topQ,]
react10=reactrwwrs[top10,]
react5=reactrwwrs[top5,]

row.names(reactQ)=reactQ$reaction
row.names(react10)=react10$reaction
row.names(react5)=react5$reaction


write.table(reactQ,file="reactQ_2.txt",row.names=F,sep="\t",quote=F)
write.table(react10,file="react10_2.txt",row.names=F,sep="\t",quote=F)
write.table(react5,file="react_5_2.txt",row.names=F,sep="\t",quote=F)

bigcorenet=corenet

bigcorenet$t0wt=0
bigcorenet$t0a=0
bigcorenet$t0b=0
bigcorenet$t03wt=0
bigcorenet$t03a=0
bigcorenet$t03b=0
bigcorenet$t36wt=0
bigcorenet$t36a=0
bigcorenet$t36b=0

bigcorenet$lr0a=0
bigcorenet$lr03a=0
bigcorenet$lr36a=0
bigcorenet$lr0b=0
bigcorenet$lr03b=0
bigcorenet$lr36b=0

for (i in 1:length(bigcorenet$from)){
  lines_dir=intersect(which(reactrwwrs$from==bigcorenet$from[i]),which(reactrwwrs$to==bigcorenet$to[i]))
  lines_inv=intersect(which(reactrwwrs$from==bigcorenet$to[i]),which(reactrwwrs$to==bigcorenet$from[i]))
  if (length(lines_dir)>0){
    to_sum=colSums(reactrwwrs[lines_dir,5:19])
    bigcorenet[i,3:17]=bigcorenet[i,3:17]+to_sum
  }
  if (length(lines_inv)>0){
    to_sum=colSums(reactrwwrs[lines_inv,5:19])
    bigcorenet[i,3:17]=bigcorenet[i,3:17]-to_sum
  }
}

bigcorenetsign=bigcorenet
bigcorenetsign[,3:17]=bigcorenet[,3:17]/abs(bigcorenet[,3:17])
bigcorenetflux=bigcorenet
bigcorenetflux[,3:17]=abs(bigcorenet[,3:17])

a=log10(bigcorenetflux[,3:11])
a[a==-Inf]=0
a[a==0]=min(a)-0.5
a=a-min(a)+0.5
a=a/max(a)
a=a*10

b=bigcorenetflux[,12:17]
minb=min(b[b!=0])
maxb=max(b)
b=b/maxb
b=b*20+0.5


corenet_t0w=corenet
corenet_t0w$arrows=rep("both",nrow(corenet))
corenet_t0w$arrows[which(bigcorenetsign$t0wt==1)]="last"
corenet_t0w$arrows[which(bigcorenetsign$t0wt==-1)]="first"
corenet_t0w$width=a$t0wt

corenet_t0a=corenet
corenet_t0a$arrows=rep("both",nrow(corenet))
corenet_t0a$arrows[which(bigcorenetsign$t0a==1)]="last"
corenet_t0a$arrows[which(bigcorenetsign$t0a==-1)]="first"
corenet_t0a$width=a$t0a

corenet_t0b=corenet
corenet_t0b$arrows=rep("both",nrow(corenet))
corenet_t0b$arrows[which(bigcorenetsign$t0b==1)]="last"
corenet_t0b$arrows[which(bigcorenetsign$t0b==-1)]="first"
corenet_t0b$width=a$t0b


corenet_t03w=corenet
corenet_t03w$arrows=rep("both",nrow(corenet))
corenet_t03w$arrows[which(bigcorenetsign$t03wt==1)]="last"
corenet_t03w$arrows[which(bigcorenetsign$t03wt==-1)]="first"
corenet_t03w$width=a$t03wt
 
corenet_t03a=corenet
corenet_t03a$arrows=rep("both",nrow(corenet))
corenet_t03a$arrows[which(bigcorenetsign$t03a==1)]="last"
corenet_t03a$arrows[which(bigcorenetsign$t03a==-1)]="first"
corenet_t03a$width=a$t03a

corenet_t03b=corenet
corenet_t03b$arrows=rep("both",nrow(corenet))
corenet_t03b$arrows[which(bigcorenetsign$t03b==1)]="last"
corenet_t03b$arrows[which(bigcorenetsign$t03b==-1)]="first"
corenet_t03b$width=a$t03b

corenet_t36w=corenet
corenet_t36w$arrows=rep("both",nrow(corenet))
corenet_t36w$arrows[which(bigcorenetsign$t36wt==1)]="last"
corenet_t36w$arrows[which(bigcorenetsign$t36wt==-1)]="first"
corenet_t36w$width=a$t36wt

corenet_t36a=corenet
corenet_t36a$arrows=rep("both",nrow(corenet))
corenet_t36a$arrows[which(bigcorenetsign$t36a==1)]="last"
corenet_t36a$arrows[which(bigcorenetsign$t36a==-1)]="first"
corenet_t36a$width=a$t36a


corenet_t36b=corenet
corenet_t36b$arrows=rep("both",nrow(corenet))
corenet_t36b$arrows[which(bigcorenetsign$t36b==1)]="last"
corenet_t36b$arrows[which(bigcorenetsign$t36b==-1)]="first"
corenet_t36b$width=a$t36b

corenet_t0w$color=rep("wt",nrow(corenet))
corenet_t0a$color=rep("-ndh-2a",nrow(corenet))
corenet_t0b$color=rep("-ndh-2b",nrow(corenet))

corenet_t0=rbind(corenet_t0w,corenet_t0a,corenet_t0b)


library(tidygraph)
library(ggraph)
library(grid)

cgraph=tbl_graph(nodes=corenode,edges=corenet_t0)
save(cgraph,file="cgraph_t0.RData")
ggraph(cgraph, layout="stress") +
  geom_edge_parallel(aes(color=color,width=width),arrow=grid::arrow(angle = 15, length = unit(0.005*corenet_t0$width,"npc"),ends=corenet_t0$arrows), start_cap = circle(3, 'mm'), end_cap = circle(3, 'mm'),sep = unit(3, "mm")) +
  geom_node_point() +
  geom_node_text(aes(label = label),nudge_x=0.2,nudge_y=-0.25, repel=T) +
  scale_edge_width(range=c(0.25,3)) +
  theme_graph()



corenet_t03w$color=rep("wt",nrow(corenet))
corenet_t03a$color=rep("-ndh-2a",nrow(corenet))
corenet_t03b$color=rep("-ndh-2b",nrow(corenet))

corenet_t03=rbind(corenet_t03w,corenet_t03a,corenet_t03b)


cgraph=tbl_graph(nodes=corenode,edges=corenet_t03)
save(cgraph,file="cgraph_t03.RData")
ggraph(cgraph, layout="stress") +
  geom_edge_parallel(aes(color=color,width=width),arrow=grid::arrow(angle = 15, length = unit(0.005*corenet_t03$width,"npc"),ends=corenet_t03$arrows), start_cap = circle(3, 'mm'), end_cap = circle(3, 'mm'),sep = unit(3, "mm")) +
  geom_node_point() +
  geom_node_text(aes(label = label),nudge_x=0.2,nudge_y=-0.25, repel=T) +
  scale_edge_width(range=c(0.25,3)) +
  theme_graph()

corenet_t36w$color=rep("wt",nrow(corenet))
corenet_t36a$color=rep("-ndh-2a",nrow(corenet))
corenet_t36b$color=rep("-ndh-2b",nrow(corenet))

corenet_t36=rbind(corenet_t36w,corenet_t36a,corenet_t36b)


cgraph=tbl_graph(nodes=corenode,edges=corenet_t36)
save(cgraph,file="cgraph_t36.RData")
ggraph(cgraph, layout="stress") +
  geom_edge_parallel(aes(color=color,width=width),arrow=grid::arrow(angle = 15, length = unit(0.005*corenet_t36$width,"npc"),ends=corenet_t36$arrows), start_cap = circle(3, 'mm'), end_cap = circle(3, 'mm'),sep = unit(3, "mm")) +
  geom_node_point() +
  geom_node_text(aes(label = label),nudge_x=0.2,nudge_y=-0.25, repel=T) +
  scale_edge_width(range=c(0.25,3)) +
  theme_graph()




