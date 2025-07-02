library(tidyverse)
library(data.table)
library(limma)
library(edgeR)
library(locfdr)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(zoo)
library(scales)
library(ineq)
library(topGO)
library(viridis)
library(ggrepel)
library(wiggleplotr)
library(universalmotif)

options(stringsAsFactors=F)
setwd("~/projects/021219_if3/current/analysis")

samples=c("IP","IP","input","input")
names(samples)=c("ATCACG","CGATGT","GAGCTA","TCTGAC")

binCounts=data.frame(fread("winCounts_w10_s5.gz",sep="\t",header=T))
binCounts=binCounts[rowSums(binCounts[,names(samples[samples=="IP"])])>=10,] #Row sum of at least 10 counts in IP sample to consider the window

rownames(binCounts)=binCounts$id
binCounts$name=paste(binCounts$chr,binCounts$start,binCounts$end,binCounts$strand,sep="_")

binCounts_colnames=colnames(binCounts)[-(1:6)]
binCounts_mat=as.matrix(binCounts[,-(1:6)])

pos=binCounts[,1:6]
rownames(pos)=pos$name
rm(binCounts)

#TMM scaling
expr=DGEList(counts=binCounts_mat)
expr=calcNormFactors(expr,method="TMM")

#Normalization factor for plotting
normFactors=expr$samples$norm.factors
names(normFactors)=rownames(expr$sample)
rm(binCounts_mat)

design=model.matrix(~0+samples)  
colnames(design)=unlist(lapply(strsplit(colnames(design),"samples"),"[[",2))
cont=makeContrasts(IP-input,IP,input,levels=colnames(design))

#Voom-limma to test IP-input significance
v=voom(expr,plot=T,design=design)
fit=lmFit(v,design)
rm(expr)

fit2=contrasts.fit(fit,cont)
rm(fit)
fit2_eb=eBayes(fit2)

tt=topTable(fit2_eb,n=Inf,coef=1,sort.by="none",confint=T) #Master table
tt$locfdr=locfdr(tt$t,bre=120,mlests=c(-2,0.5),df=10)$fdr #Local FDR
tt=cbind(pos[,c(1:4,6)],fit2_eb$coefficients,tt)
rm(fit2)

#Add annotations using TxDb
#txdb=makeTxDbFromEnsembl(organism="Mus musculus",server="useastdb.ensembl.org") #Ensembl txdb
#saveDb(txdb, "./txdb.rds")
txdb = loadDb("./txdb.rds") #Load already saved txdb
seqlevels(txdb)=seqlevels0(txdb)
seqlevels(txdb)=as.character(c(1:19,"X","Y","MT"))
seqlevels(txdb)[seqlevels(txdb)%in%as.character(c(1:19,"X","Y","MT"))]=c(paste("chr",seqlevels(txdb)[seqlevels(txdb)%in%as.character(c(1:19,"X","Y"))],sep=""),"chrM")

#Useful conversions
ensembl2eg=unlist(as.list(org.Mm.egENSEMBL2EG))
eg2ensembl=unlist(as.list(org.Mm.egENSEMBL))
eg2symbol=unlist(as.list(org.Mm.egSYMBOL))
ensembl2symbol=eg2symbol[ensembl2eg]
names(ensembl2symbol)=names(ensembl2eg)

#Overlaps with exons: 5'UTR, ORF, or 3'UTR
pos_gr=GRanges(seqnames=Rle(pos[,1]),IRanges(start=pos[,2],end=pos[,3]),Rle(strand(pos[,6])))
tx=as.data.frame(transcripts(txdb,columns=c("tx_name","gene_id","tx_id","tx_chrom")))
tx$gene_id=unlist(tx$gene_id)
rownames(tx)=tx$tx_name
utr5=fiveUTRsByTranscript(txdb)
utr3=threeUTRsByTranscript(txdb)
cds=cdsBy(txdb)
exons=exonsBy(txdb,by=c("gene"))
longesttx=max(width(transcriptsBy(txdb,by="gene"))) #longest transcript

enst2ensg=unique(tx[,c("tx_name","gene_id")]) #ENST to ENSG mapping

#Map overlaps
pos_gr_utr5=as.data.frame(findOverlaps(pos_gr,utr5))
pos_gr_utr3=as.data.frame(findOverlaps(pos_gr,utr3))
pos_gr_cds=as.data.frame(findOverlaps(pos_gr,cds))

pos_gr_exon=as.data.frame(findOverlaps(pos_gr,exons))
pos_gr_exon$ensg=names(exons)[pos_gr_exon$subjectHits]
pos_gr_exon$symbol=ensembl2symbol[pos_gr_exon$ensg]
pos_gr_exon_dups=subset(pos_gr_exon %>% group_by(queryHits) %>% summarise(freq=length(queryHits)),freq>1)$queryHits
pos_gr_exon_collapsed=data.frame(rbind(subset(pos_gr_exon,queryHits%in%pos_gr_exon_dups) %>% group_by(queryHits) %>% summarise(subjectHits=NA,ensg=paste(ensg,collapse=","),symbol=paste(symbol,collapse=",")),subset(pos_gr_exon,!queryHits%in%pos_gr_exon_dups)))
pos_gr_exon$longest=longesttx[pos_gr_exon$ensg] #longest TX
rownames(pos_gr_exon_collapsed)=pos_gr_exon_collapsed$queryHits

pos_gr_exon_longest=data.frame(pos_gr_exon %>% group_by(queryHits) %>% top_n(longest,n=1))
rownames(pos_gr_exon_longest)=pos_gr_exon_longest$queryHits

utr5overlap=unique(pos_gr_utr5$queryHits)
utr3overlap=unique(pos_gr_utr3$queryHits)
cdsoverlap=unique(pos_gr_cds$queryHits)
exonoverlap=unique(as.data.frame(pos_gr_exon)$queryHits)

#Add the overlap annotations to the master table
tt$cds=F
tt$utr3=F
tt$utr5=F
tt$exon=F
tt$cds[cdsoverlap]=T
tt$utr3[utr3overlap]=T
tt$utr5[utr5overlap]=T
tt$exon[exonoverlap]=T

tt$longest.ensg=NA
tt$longest.symbol=NA
tt$all.ensg=NA
tt$all.symbol=NA

tt_useid=as.character(1:nrow(tt))[as.character(1:nrow(tt))%in%pos_gr_exon$queryHits]
tt[as.numeric(tt_useid),"longest.ensg"]=pos_gr_exon_longest[tt_useid,"ensg"]
tt[as.numeric(tt_useid),"longest.symbol"]=pos_gr_exon_longest[tt_useid,"symbol"]
tt[as.numeric(tt_useid),"all.ensg"]=pos_gr_exon_collapsed[tt_useid,"ensg"]
tt[as.numeric(tt_useid),"all.symbol"]=pos_gr_exon_collapsed[tt_useid,"symbol"]

#Add Ptch1 RACE annotations manually for plotting
tt[tt$name%in%subset(tt,(chr=="chr13")&(end<=63566578)&(start>=63565735))$name,"exon"]=T
tt[tt$name%in%subset(tt,(chr=="chr13")&(end<=63566578)&(start>=63565735))$name,"utr5"]=T
tt[tt$name%in%subset(tt,(chr=="chr13")&(end<=63566578)&(start>=63565735))$name,"longest.symbol"]="Ptch1"
tt[tt$name%in%subset(tt,(chr=="chr13")&(end<=63566578)&(start>=63565735))$name,"longest.ensg"]="ENSMUSG00000021466"
tt[tt$name%in%subset(tt,(chr=="chr13")&(end<=63566578)&(start>=63565735))$name,"all.symbol"]="Ptch1"
tt[tt$name%in%subset(tt,(chr=="chr13")&(end<=63566578)&(start>=63565735))$name,"all.ensg"]="ENSMUSG00000021466"

tt$sig=F #Per window significance column
fdrthres=0.05 #FDR threshold at 0.05
tt[subset(tt,(locfdr<=fdrthres)&(logFC>=2)&(IP>=(-2.5)))$name,"sig"]=T #Combination of minimum IP counts, minimum fold-change, FDR cutoff

tt_exon=subset(tt,exon) #Windows overlapping exons only

#Scatter of significant hits
fig_scatter=ggplot(tt,aes(x=input,y=IP))+theme_classic()+
  geom_point(data=subset(tt,(locfdr>fdrthres)|(logFC<2)|(IP<(-2.5))),shape=19,alpha=0.1,size=0.5,stroke=0)+
  geom_point(data=subset(tt,(locfdr<=fdrthres)&(logFC>=2)&(IP>=(-2.5))),stroke=0,shape=19,size=0.5,alpha=0.3,colour="#d82433")+
  coord_fixed(ratio=1)+scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+theme(legend.title=element_blank())

pdf("./plots/scatter.pdf",width=6,height=6)
print(fig_scatter)
dev.off()

#5UTR/ORF/3UTR proportions
utrprop=rep(NA,nrow(tt_exon))
utrprop[(tt_exon$utr5)&(!tt_exon$utr3)&(!tt_exon$cds)]="5UTR"
utrprop[(tt_exon$utr3)&(!tt_exon$utr5)&(!tt_exon$cds)]="3UTR"
utrprop[(tt_exon$cds)&(!tt_exon$utr5)&(!tt_exon$utr3)]="ORF"
utrprop[is.na(utrprop)]="AMBIGUOUS"
chitest=chisq.test(table(utrprop,(tt_exon$locfdr<=fdrthres)&(tt_exon$logFC>=0))[-3,])
tableDf=data.frame(Type=names(table(utrprop)),All=as.numeric(table(utrprop)),IP=as.numeric(table(utrprop[(tt_exon$locfdr<=fdrthres)&(tt_exon$logFC>=0)])))
tableDf=subset(tableDf,Type%in%c("5UTR","ORF","3UTR"))
tableDf[,-1]=tableDf[,-1]/matrix(rep(colSums(tableDf[,-1]),nrow(tableDf)),nrow(tableDf),2,byrow=T)
tableDf=reshape::melt(tableDf,id.vars=c("Type"))
colnames(tableDf)=c("Type","Sample","Freq")
tableDf$Type=factor(tableDf$Type,levels=(rev(c("5UTR","ORF","3UTR"))))

#5UTR/ORF/3UTR proportions plot
fig_percents=ggplot(tableDf,aes(x=Sample,y=Freq,fill=Type))+theme_classic()+
  geom_bar(position="fill",stat="identity",colour="black")+
  geom_text(data=subset(tableDf,Type=="5UTR"),aes(y=Freq-0.25,label=paste0(round(Freq*100,1),"%")),nudge_y=0.1)+
  scale_y_continuous(expand=c(0,0),labels=scales::percent)+
  scale_fill_manual(values=c("#377eb8","#4daf4a","#e41a1c"))+
  coord_flip()+ylab("Percent")

pdf("./plots/percents.pdf",width=8,height=3)
print(fig_percents)
dev.off()

#Calculate gini index for 10 neighboring windows around each window
tt_exon_siggenes=unique(subset(tt_exon,sig)$longest.ensg)
tt_exon_gini=
  data.frame(subset(tt_exon,(longest.ensg%in%tt_exon_siggenes)&(IP>=(-2.5))) 
             %>% group_by(longest.ensg) 
             %>% mutate(gini=rollapplyr(rescale(2^IP),10,Gini,partial=T,fill=NA,align="center"),
                        peak=rollapplyr(rescale(IP),10,which.max,partial=T,fill=NA,align="center"),
                        wincount=rollapplyr(IP,10,length,partial=T,fill=NA,align="center")))

#Hierarchical selection of top sharp peaks by binning by IP signal
tt_exon_gini_center=
  subset(tt_exon_gini,
         (locfdr<=fdrthres)&(logFC>=2)&(IP>=-0.5)&(!is.na(gini))&(peak>=(0.3*wincount))&(peak<=(0.7*wincount))&(wincount==10)) #Take peaks that are towards the middle
tt_exon_gini_center$ginirank_quantile_all=rank(tt_exon_gini_center$gini)/nrow(tt_exon_gini_center)*100
tt_exon_gini_center$IPbin=ntile(tt_exon_gini_center$IP,10) #Divide into N bins by I
tt_exon_gini_center=data.frame(
  tt_exon_gini_center %>% group_by(IPbin) %>% 
    mutate(ginirank=rank(gini),ginirank_quantile=rank(gini)/length(gini)*100))#Rank of gini index within each bin

rownames(tt_exon_gini_center)=tt_exon_gini_center$name

#Plot IP signal vs Gini index to see this
tt_exon_gini_center$label=""
tt_exon_gini_center.plotlabels=data.frame(subset(tt_exon_gini_center,longest.symbol%in%c("Ptch1","Gli3")) %>% group_by(longest.symbol) %>% dplyr::slice(which.max(ginirank)))$name
tt_exon_gini_center[tt_exon_gini_center.plotlabels,"label"]=tt_exon_gini_center[tt_exon_gini_center.plotlabels,"longest.symbol"]

fig_ipgini=ggplot(tt_exon_gini_center,aes(x=IP,y=(gini),colour=ginirank_quantile))+
  theme_classic()+geom_point(shape=19,size=0.5,alpha=0.75,stroke=0)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_colour_gradient(name="Quantile\nwithin bin",low="black",high="red")+
  ylab("Gini index")+xlab("Log2 IP")+
  geom_text_repel(aes(label=label),min.segment.length=0,colour="black",max.overlaps=Inf)

pdf("./plots/ip_wingini.pdf",width=6,height=4)
print(fig_ipgini)
dev.off()

#Sharp peak cutoff (>=60 percentile in each bin)
tt$sig.sharp=F
tt[subset(tt_exon_gini_center,ginirank_quantile>=60)$name,"sig.sharp"]=T
tt$gini=NA
tt[tt_exon_gini_center$name,"gini"]=tt_exon_gini_center$gini
tt$ginirankpercentile_all=NA
tt[tt_exon_gini_center$name,"ginirankpercentile_all"]=tt_exon_gini_center$ginirank_quantile_all
tt$ginibin=NA
tt[tt_exon_gini_center$name,"ginibin"]=tt_exon_gini_center$IPbin
tt$ginirank=NA
tt[tt_exon_gini_center$name,"ginirank"]=tt_exon_gini_center$ginirank
tt$ginipercentile=NA
tt[tt_exon_gini_center$name,"ginipercentile"]=tt_exon_gini_center$ginirank_quantile

tt_exon=subset(tt,exon)

#Output the tables
write.table(tt,file="./tables/all_windows.tsv",row.names=F,col.names=T,quote=F,sep="\t")
write.table(tt_exon,file="./tables/exononly_windows.tsv",row.names=F,col.names=T,quote=F,sep="\t")
write.table(unique(tt_exon[,c("longest.ensg","longest.symbol")]),file="./tables/detected_genes.tsv",row.names=F,col.names=T,quote=F,sep="\t")

#Coverage plotting
txlist=as.data.frame(transcripts(txdb))$tx_name
exonsByTx=exonsBy(txdb,by="tx",use.names=T)
cdssByTx=cdsBy(txdb,by="tx",use.names=T)
utr5ByTx=unlist(fiveUTRsByTranscript(txdb,use.names=T))

#Patch Gli3 and Ptch1 5' based on RACE for plotting
start(exonsByTx[c("ENSMUST00000141194","ENSMUST00000110510")][[2]][1])=15463705
start(exonsByTx[c("ENSMUST00000141194","ENSMUST00000110510")][[1]][1])=15463235
end(exonsByTx["ENSMUST00000021921"][[1]][1])=63566578

transcript_metadata=as.data.frame(transcripts(txdb,columns=c("tx_name","gene_id","tx_type")))[,c("tx_name","gene_id","strand","tx_type")]
transcript_metadata$gene_id=unlist(transcript_metadata$gene_id)
colnames(transcript_metadata)=c("transcript_id","gene_id","strand","transcript_type")
transcript_metadata$gene_name=ensembl2symbol[transcript_metadata$gene_id]

sample_data=data_frame(sample_id=gsub("\\.","-",binCounts_colnames),condition=factor(samples,unique(samples)),scaling_factor=1/normFactors)
sample_data=mutate(sample_data, track_id = condition, colour_group = condition)

sample_data_plus=sample_data %>% dplyr::mutate(bigWig=paste0(getwd(),"/bw/",sample_id, ".minus.bw"))
sample_data_minus = sample_data %>% dplyr::mutate(bigWig=paste0(getwd(),"/bw/",sample_id, ".plus.bw"))

plotCov=function(plotsymbol, usetx=NULL, useexons=NULL, usecdss=NULL, utr5=F) {
  if(!is.null(usetx)) {
    selected_transcripts=usetx
  } else{
    selected_transcripts=(transcript_metadata %>% filter(gene_name==plotsymbol))$transcript_id
  }
  tx_ids=intersect(selected_transcripts, txlist)
  if(length(tx_ids)==0) {
    return(NA)
  }
  plotstrand=as.character(strand(exonsByTx[tx_ids[1]])[[1]])[1]
  colourpalette=c("#d82433","#000000","#183f6d","#000000")
  
  if(plotstrand=="+") {
    sample_data_use=sample_data_plus
  } else{
    sample_data_use=sample_data_minus
  }
  
  plotheights=c(1,1)
  plotregion=c(as.character(seqnames(range(GenomicRanges::reduce(unlist(exonsByTx[tx_ids]))))),start(range(GenomicRanges::reduce(unlist(exonsByTx[tx_ids]))))-100,end(range(GenomicRanges::reduce(unlist(exonsByTx[tx_ids]))))+100,as.character(strand(range(GenomicRanges::reduce(unlist(exonsByTx[tx_ids])))))) #+-100nt
  plotobj=list()
  
  usecds=intersect(tx_ids,names(cdssByTx))
  if(length(usecds)==0) {
    usecds=NULL
  } else{
    usecds=cdssByTx[usecds]
  }
  
  if(!is.null(useexons)) {
    plotexons=useexons
  } else{
    plotexons=exonsByTx[tx_ids]
  }
  
  if(utr5)
  {
    use_itx=NULL
    plotexons=exonsByTx[tx_ids]
    plotcdss=usecds
    for(itx in names(plotexons)) {
      plotexons[[itx]]=plotexons[[itx]][plotexons[[itx]]$exon_name%in%utr5ByTx$exon_name,]
      plotcdss[[itx]]=usecds[[itx]][1,]
      if(length(plotexons[[itx]])>0)
      {
        use_itx=c(use_itx,itx)
      }
    }
    plotexons=plotexons[intersect(use_itx,intersect(names(plotexons),names(plotcdss)))]
    plotcdss=plotcdss[intersect(use_itx,intersect(names(plotexons),names(plotcdss)))]
  }
  
  if(!is.null(usecdss)) {
    plotcdss=usecdss
  } else{
    if(!utr5){
      plotcdss=usecds
    }
  }
  
  plotobj[["cov"]]=plotCoverage(plotexons,plotcdss,transcript_metadata,sample_data_use,rescale_introns=T,heights=plotheights,transcript_label=F,mean_only=T,alpha=1,fill_palette=colourpalette,coverage_type="area",plot_fraction=1,return_subplots_list=T,region_coords=as.numeric(plotregion[2:3]),flanking_length=c(0,0))
  plotobj[["region"]]=plotregion
  return(plotobj)
}

#ENSMUST00000141194; Gli3
gli3_all=plotCov("Gli3",c("ENSMUST00000141194","ENSMUST00000110510"))
gli3_all_top=gli3_all$cov$coverage_plot+scale_y_continuous(expand=c(0,0),limits=c(0,55))#+scale_x_continuous(expand=c(0,0))
gli3_all_bot=gli3_all$cov$tx_structure
fig_gli3_all=cowplot::plot_grid(gli3_all_top,gli3_all_bot,ncol=1,align="v",rel_heights=c(0.7,0.3))

pdf("./plots/gli3_all.pdf",width=6,height=4)
print(fig_gli3_all)
dev.off()

#Gli3 5'UTR only
gli3_utr5=plotCov("Gli3",c("ENSMUST00000141194","ENSMUST00000110510"),GRangesList(gli3=exonsByTx[["ENSMUST00000110510"]][1:2]),GRangesList(gli3=cdssByTx[["ENSMUST00000110510"]][1:2]))
gli3_utr5_top=gli3_utr5$cov$coverage_plot+scale_y_continuous(expand=c(0,0),limits=c(0,55))
gli3_utr5_bot=gli3_utr5$cov$tx_structure
fig_gli3_utr5=cowplot::plot_grid(gli3_utr5_top,gli3_utr5_bot,ncol=1,align="v",rel_heights=c(0.7,0.3))

pdf("./plots/gli3_utr5.pdf",width=6,height=4)
print(fig_gli3_utr5)
dev.off()

#Ptch1
plotsymbol="Ptch1"
ptch1_all=plotCov("Ptch1","ENSMUST00000021921")
ptch1_all_top=ptch1_all$cov$coverage_plot+scale_y_continuous(expand=c(0,0),limits=c(0,40))
ptch1_all_bot=ptch1_all$cov$tx_structure
fig_ptch1_all=cowplot::plot_grid(ptch1_all_top,ptch1_all_bot,ncol=1,align="v",rel_heights=c(0.7,0.3))

pdf("./plots/ptch1_all.pdf",width=6,height=4)
print(fig_ptch1_all)
dev.off()

#Ptch1 5'UTR only
ptch1_utr5=plotCov("Ptch1","ENSMUST00000021921",GRangesList(Ptch1=exonsByTx[["ENSMUST00000021921"]][1]),GRangesList(Ptch1=cdssByTx[["ENSMUST00000021921"]][1]))
ptch1_utr5_top=ptch1_utr5$cov$coverage_plot+scale_y_continuous(expand=c(0,0),limits=c(0,40))
ptch1_utr5_bot=ptch1_utr5$cov$tx_structure
fig_ptch1_utr5=cowplot::plot_grid(ptch1_utr5_top,ptch1_utr5_bot,ncol=1,align="v",rel_heights=c(0.7,0.3))

pdf("./plots/ptch1_utr5.pdf",width=6,height=4)
print(fig_ptch1_utr5)
dev.off()

#Output sharp peaks for motif analysis
sharpbed=GRanges(subset(tt_exon,sig.sharp)[,c("chr","start","end","name","logFC","strand","longest.ensg")])
start(sharpbed)=start(sharpbed)-25 #+-25nt 
end(sharpbed)=end(sharpbed)+25
sharpbed_reduced=reduce(sharpbed,with.revmap=T) #Merge overlapping windows
sharpbed_reduced_out=as.data.frame(sharpbed_reduced)
sharpbed_reduced_out=data.frame(chr=sharpbed_reduced_out$seqnames,start=sharpbed_reduced_out$start,end=sharpbed_reduced_out$end,name=NA,score=0,strand=sharpbed_reduced_out$strand,longest.ensg=sharpbed[unlist(lapply(sharpbed_reduced_out$revmap,"[[",1))]$longest.ensg,longest.symbol=ensembl2symbol[sharpbed[unlist(lapply(sharpbed_reduced_out$revmap,"[[",1))]$longest.ensg])
sharpbed_reduced_out$name=1:nrow(sharpbed_reduced_out)
write.table(sharpbed_reduced_out[,1:6],file="sharp.bed",row.names=F,col.names=F,quote=F,sep="\t")
write.table(sharpbed_reduced_out,file="sharp.tsv",sep="\t",row.names=F,col.names=T,quote=F)

#Background
bgbed=GRanges(subset(tt,(CI.R<=(-1))&(IP>=(-2.5)&(!sig.sharp)))[,c("chr","start","end","name","logFC","strand")]) #lower CI for fold change less than log2 -1, and IP signal larger than the minimum count, and no sharp peak
start(bgbed)=start(bgbed)-25 #+-25nt
end(bgbed)=end(bgbed)+25
bgbed_reduced=as.data.frame(reduce(bgbed))
bgbed_reduced=data.frame(chr=bgbed_reduced$seqnames,start=bgbed_reduced$start,end=bgbed_reduced$end,name=NA,score=0,strand=bgbed_reduced$strand)
bgbed_reduced$name=1:nrow(bgbed_reduced)
bgbed_reduced=subset(bgbed_reduced,chr!="gi|307829144|gb|GU372691.1|")
write.table(bgbed_reduced,file="bg.bed",row.names=F,col.names=F,quote=F,sep="\t")

#Read in MEME results and annotate
memeres=read_meme("./sharp_meme/meme.txt",readsites=T,readsites.meta=T)
meme_dfs=list()
for(meme in names(memeres$sites))
{
  meme_df=memeres$sites.meta[[meme]]
  meme_df$longest.ensg=sharpbed_reduced_out[memeres$sites.meta[[meme]]$Sequence,"longest.ensg"]
  meme_df$longest.symbol=sharpbed_reduced_out[memeres$sites.meta[[meme]]$Sequence,"longest.symbol"]
  meme_df$chr=sharpbed_reduced_out[memeres$sites.meta[[meme]]$Sequence,"chr"]
  meme_df$start=sharpbed_reduced_out[memeres$sites.meta[[meme]]$Sequence,"start"]
  meme_df$end=sharpbed_reduced_out[memeres$sites.meta[[meme]]$Sequence,"end"]
  meme_df$strand=sharpbed_reduced_out[memeres$sites.meta[[meme]]$Sequence,"strand"]
  meme_df$motif=meme
  meme_dfs[[meme]]=meme_df
}
meme_dfs=do.call(rbind,meme_dfs)
write.table(meme_dfs,file="./tables/meme_motif_sites.tsv",sep="\t",row.names=F,col.names=T,quote=F) #Write annotated MEME result table

#GO term enrichment using topGO
runGO=function(sigwin,allwin,ensg)
{
  sigsel=function(gv) {return(gv==1)}
  
  if(ensg){
    siggenes=rep(0,length(allwin))
    names(siggenes)=allwin
    siggenes[allwin%in%sigwin]=1
  }
  else {
    allgenes=unique(subset(tt,name%in%allwin)$longest.ensg)
    siggenes=rep(0,length(allgenes))
    names(siggenes)=allgenes
    siggenes[subset(tt,name%in%sigwin)$longest.ensg]=1
  }
  
  topGO=new("topGOdata",ontology="BP",allGenes=siggenes,geneSel=sigsel,nodeSize=10,annot=annFUN.org,mapping="org.Mm.eg.db",ID="ensembl")
  resultFC=runTest(topGO,algorithm="classic",statistic="fisher")
  resultFW01=runTest(topGO,algorithm="weight01",statistic="fisher")
  allGO=usedGO(topGO)
  
  allResF=GenTable(topGO,fisherWeight01=resultFW01,fisherClassic=resultFC,topNodes=length(allGO),numChar=75)
  allResF$ser=(allResF$Significant/allResF$Expected) #Enrichment ratio
  allResF$fisherClassic=score(resultFC)[allResF$GO.ID]
  allResF$fisherWeight01=score(resultFW01)[allResF$GO.ID]
  allResF$FCBH=p.adjust(score(resultFC)[allResF$GO.ID],method="BH") #BH FDR estimation
  
  allResF_nz=subset(allResF,Significant>0)
  
  #Get gene names in each term
  gogenes=list()
  gosymbols=list()
  for (gorow in 1:nrow(allResF_nz)){
    goterm=allResF_nz[gorow,1]
    gogenes[[goterm]]=intersect(unlist(genesInTerm(topGO,goterm)),names(siggenes[siggenes==1]))
    gosymbols[[goterm]]=ensembl2symbol[gogenes[[goterm]]]
  }
  
  return(list(allResF_nz, gogenes, gosymbols))
}

#GO terms for sharp peaks vs all expressed genes
go_all=runGO(subset(tt_exon,sig.sharp)$name,tt_exon$name,F)
go_all_out=subset(go_all[[1]],(fisherWeight01<=0.05)&(Significant>=3)&(ser>=2)&(FCBH<=0.05))
go_all_out_melt=melt(go_all_out[,c(2,4,5)])

#Make bar plot
fig_go_all=ggplot(go_all_out_melt,aes(x=Term,y=value,fill=variable))+theme_classic()+geom_bar(stat="identity",position=position_dodge())+scale_y_continuous(expand=c(0,0))+ylab("Count")+scale_fill_manual(values=c("#377eb8","#e41a1c"),name="")+coord_flip()

pdf("./plots/goterm_vs_all_background.pdf",width=7,height=8,bg="transparent") 
print(fig_go_all)
dev.off()

#Output enriched terms as table
go_all_out_table=go_all_out
go_all_out_table$genes=unlist(lapply(go_all[[3]][go_all_out$GO.ID],paste,collapse=","))
colnames(go_all_out_table)[c(9,11)]=c("Observed/Expected","Genes")
go_all_out_table=go_all_out_table[order(go_all_out_table[,9],decreasing=T),]
write.table(go_all_out_table,file="./tables/go_genes_all_background.tsv",row.names=F,col.names=T,quote=F,sep="\t")

#GO terms for sharp peaks with motif
go_sharp_motif2=runGO(unique(subset(meme_dfs,motif=="YUUYYUKYUYYYBSHBYYYB")$longest.ensg),unique(tt_exon$longest.ensg),T)
go_sharp_motif2_out=subset(go_sharp_motif2[[1]],(fisherWeight01<=0.1)&(Significant>=3)&(ser>=2)&(FCBH<=0.05))
go_sharp_motif2_out_melt=melt(go_sharp_motif2_out[,c(2,4,5)])

#Bar plot
fig_go_sharp_motif2=ggplot(go_sharp_motif2_out_melt,aes(x=Term,y=value,fill=variable))+theme_classic()+geom_bar(stat="identity",position=position_dodge())+scale_y_continuous(expand=c(0,0))+ylab("Count")+scale_fill_manual(values=c("#377eb8","#e41a1c"),name="")+coord_flip()

pdf("./plots/goterm_motif2_all_background.pdf",width=7,height=8,bg="transparent") 
print(fig_go_sharp_motif2)
dev.off()

#Output enriched terms as table
go_sharp_motif2_out_table=go_sharp_motif2_out
go_sharp_motif2_out_table$genes=unlist(lapply(go_sharp_motif2[[3]][go_sharp_motif2_out$GO.ID],paste,collapse=","))
colnames(go_sharp_motif2_out_table)[c(9,11)]=c("Observed/Expected","Genes")
go_sharp_motif2_out_table=go_sharp_motif2_out_table[order(go_sharp_motif2_out_table[,9],decreasing=T),]
write.table(go_sharp_motif2_out_table,file="./tables/go_genes_motif2.tsv",row.names=F,col.names=T,quote=F,sep="\t")

write.table(subset(go_sharp_motif2[[1]]),file="./tables/go_genes_motif2_all_terms.tsv",row.names=F,col.names=T,quote=F,sep="\t")


#Ribo-seq analysis
ribo_counts=fread("./ribo_counts.gz")
ribo_counts=subset(ribo_counts,(name=="ORF_COUNT")&(strand=="+")) #Only positive strand and only ORF without start/stop region
ribo_samples=paste(unlist(lapply(strsplit(colnames(ribo_counts)[-(1:6)],"_"),"[[",2)),
                   unlist(lapply(strsplit(colnames(ribo_counts)[-(1:6)],"_"),"[[",3)),unlist(lapply(strsplit(colnames(ribo_counts)[-(1:6)],"_"),"[[",4)),sep="_")
ribo_counts[,-(1:6)]
colnames(ribo_counts)[-(1:6)]=ribo_samples
ribo_counts=ribo_counts[,-c(4,5,6)]

#Combine libraries for each sample
ribo_counts=ribo_counts %>% 
  pivot_longer(!c("chr","start","end"),values_to="count") %>% 
  group_by(chr,start,end,name) %>% summarize(sum=sum(count)) %>% 
  pivot_wider(id_cols=c("chr","start","end"),names_from="name",values_from="sum")

ribo_counts=ribo_counts[rowSums(ribo_counts[,-c(1:3)])>=10,] #Minimum 10 summed row counts
ribo_counts$chr=unlist(lapply(strsplit(ribo_counts %>% pull(chr),"\\."),"[[",1))
ribo_counts$ensg=enst2ensg[ribo_counts%>%pull(chr),2] #ENSG mapping
ribo_counts$rowsum=rowSums(ribo_counts[,4:15])
ribo_counts=(ribo_counts%>%group_by(ensg)%>%slice_max(rowsum,with_ties=F))[,1:15] #Take maximum count ENSG only

ribo_expr=DGEList(counts=ribo_counts[,-c(1:3)])
ribo_expr=calcNormFactors(ribo_expr,method="TMM") #TMM scaling
rownames(ribo_expr)=ribo_counts%>%pull(chr)

ribo_samples=colnames(ribo_counts)[-(1:3)]
ribo_batch=unlist(lapply(strsplit(ribo_samples,"_"),"[[",1)) #Pair each batch of Ribo/RNA
ribo_samples=paste(unlist(lapply(strsplit(ribo_samples,"_"),"[[",2)),unlist(lapply(strsplit(ribo_samples,"_"),"[[",3)),sep="_")
ribo_design=model.matrix(~0+ribo_samples+ribo_batch)
colnames(ribo_design)[1:4]=unlist(lapply(strsplit(colnames(ribo_design)[1:4],"ribo_samples"),"[[",2))
ribo_cont=makeContrasts(mut_Ribo-mut_RNA-WT_Ribo+WT_RNA,mut_RNA-WT_RNA,mut_Ribo-WT_Ribo,levels=colnames(ribo_design)) #Contrasts

#Voom-limma to test IP-input significance
ribo_v=voom(ribo_expr,plot=T,design=ribo_design)

ribo_fit=lmFit(ribo_v,ribo_design)
ribo_fit2=contrasts.fit(ribo_fit,ribo_cont)
ribo_fit2_eb=eBayes(ribo_fit2)

#Riboseq table
ribo_df=data.frame(enst=rownames(ribo_fit2_eb),ensg=enst2ensg[rownames(ribo_v),2],symbol=ensembl2symbol[enst2ensg[rownames(ribo_v),2]])
ribo_df2=cbind(ribo_fit$coefficients[,1:4],ribo_fit2_eb$coefficients,ribo_fit2_eb$t)
colnames(ribo_df2)[5:10]=c("coef_te","coef_rna","coef_ribo","t_te","t_rna","t_ribo")
ribo_df=cbind(ribo_df,ribo_df2)
rm(ribo_df2)

#Locfdr estimates
ribo_df$locfdr_te=locfdr(ribo_df$t_te,bre=150,mlests=c(1,0.75),nulltype=1,pct=0.00001)$fdr #TE
ribo_df$locfdr_rna=locfdr(ribo_df$t_rna,bre=150,mlests=c(0,0.75),nulltype=1,pct=0.0001,df=10)$fdr #RNA only
ribo_df$locfdr_ribo=locfdr(ribo_df$t_ribo,bre=150,mlests=c(1,0.75),nulltype=1,pct=0.001)$fdr #Ribo only

rownames(ribo_df)=ribo_df$ensg

#Intersect with 5'UTR sharp peak genes
ribo_df$sig.sharp=NA
ribo_df[intersect(unique(subset(tt_exon)$longest.ensg),ribo_df$ensg),"sig.sharp"]=F
ribo_df[intersect(unique(subset(tt_exon,sig.sharp&utr5)$longest.ensg),ribo_df$ensg),"sig.sharp"]=T

write.table(ribo_df,file="eif3_ribo.tsv",sep="\t",row.names=F,col.names=T,quote=F)

#Generate "controls" with matching expression level
ribo_df_rna_sel=subset(ribo_df,!(ensg%in%tt_exon_siggenes))
ribo_df_rna_sel_ctrl=NULL
for(ii in ribo_df[which(ribo_df$sig.sharp),"ensg"])
{
  #closest RNA expression to each sharp peak ribo_df record
  iii=which.min(abs(ribo_df_rna_sel$WT_RNA-ribo_df[ii,]$WT_RNA))
  ribo_df_rna_sel_ctrl=rbind(ribo_df_rna_sel_ctrl,ribo_df_rna_sel[iii,])
  ribo_df_rna_sel=ribo_df_rna_sel[-iii,]
  #subset(ribo_df,!ensg%in%tt_exon_siggenes)
}
ribo_df_rna_sel_ctrl$sig.sharp=F

wilcox.test(ribo_df_rna_sel_ctrl$t_te,subset(ribo_df,sig.sharp)$t_te)

#Plot distributions of T values for TE differences
pdf("clip_ribo_enrich.pdf",width=5,height=3)
ggplot(rbind(ribo_df_rna_sel_ctrl,subset(ribo_df,sig.sharp)),aes(x=t_te,colour=sig.sharp))+theme_classic()+geom_density(adjust=1.5)+scale_colour_manual(values=c("#377eb8","#e41a1c"),name="",labels=c("No EIF3 binding","Sharp 5'UTR peaks"))+ylab("Density")+xlab(expression("T-stat, TE"[mut-WT]))+annotate("text",label=expression(p==5.95%*%10^-15),x=4.5,y=0.15)+theme(legend.position=c(0.85,0.85))
dev.off()