# Stage 2 pipeline development
# read in batched results into a list 
# version: 4
# major revisions: remove the unnessary variables, clean D colnames
#library(shape)

options(stringsAsFactors=FALSE)
load("bin/TCGA_result.D"); tcga=res
inFolder=outFolder="workingFolder/"

if (file.access(outFolder)!=0) { dir.create(outFolder) }

# prepare
registry=read.delim("bin/Nanostring_Probe_Registry.txt")
aliasKey=aliasName=c()
for (i in which(registry$Alias!="")) {
	tmp=strsplit(registry$Alias[i],",")[[1]]
	tmp=gsub("[()]","",tmp)
	aliasKey=c(aliasKey,tmp)
	aliasName=c(aliasName,rep(registry$ProbeName[i],length(tmp)))
}

#folders=list.dirs(path=inFolder,full=TRUE)[-1] 
#folders=paste(folders,"/",sep="")
#folders=folders[grep("Huse_4.10.14",folders)]						#############folder selection
ref_directory = read.table("dir.in", sep="\t", stringsAsFactors=FALSE)[[1]]
folders=ref_directory
batch=list()
for (i in 1:length(folders)) {
	folder=folders[i]
    load(file=paste(folder,"/Stage1Result.D",sep=""))
	batch[[i]]=Stage1Result
}

uProbes=registry$ProbeName 

R=D=data.frame(jnk=rep(0,length(uProbes)))
Stems=Samples=c()
for (i in 1:length(batch)) {
	d=batch[[i]]$norm
	r=batch[[i]]$data
	probes=batch[[i]]$normProbes
	samples=make.names(batch[[i]]$nanoSamples)
	samplesR=batch[[i]]$samplesR
	di=match(uProbes,probes)
	for (j in 1:ncol(d)) {
		D[,samples[j]]=d[di,j]
		R[,samples[j]]=r[di,j]		
	}
	Stems=c(Stems,rep(batch[[i]]$stem,ncol(d)))
	Samples=c(Samples,samples)
}
D=D[,-1]+1  # normalized data table
R=R[,-1]+1  # raw data table
DC=D-rowMeans(D,na.rm=T) # normalized, mean-centered data table

# Cancer Genes
cgenes=subset(registry,CancerGene)$ProbeName[1:16]

loci=unique(registry$LocusID); loci=loci[-which(loci=="")]
probes=uProbes
start<-proc.time()
res=list()
for (i in 1:ncol(D)) {
	res[[i]]=list()
### the whole assay failed
	if (all(is.na(D[,i]))) next

### QC Report data
	cm=colMeans(R,na.rm=T)
	res[[i]]$QC=list(cm=cm)
	
### Cancer Genes Report data
	res[[i]]$CG=list()
	res[[i]]$CG$exp=list()	
	for (gg in cgenes) {
		di=which(uProbes==gg)
		vv=as.numeric(D[di,])
		if (!is.finite(vv[i])) next
		res[[i]]$CG$exp[[gg]]=vv
	}
	
### Locus Report data
	res[[i]]$LL=list()
	res[[i]]$LL$exp=list()
	for (ll in loci) {
		gg=registry$ProbeName[which(registry$LocusID==ll)]
		di=which(uProbes %in% gg)
		if (length(di)==1) {
			vv=log10(D[di,i]) 
		} else {
			vv=mean(log10(D[di,i]),na.rm=T) 
		}
		if (!is.finite(vv)) next
		if (length(di)==1) {
			v=na.omit(as.numeric(log10(D[di,-i])))
		} else {
			v=colMeans(log10(D[di,-i]),na.rm=T)
		}
		res[[i]]$LL$exp[[ll]]=list(v=v,vv=vv)
	}
	
### Variant Report data
	res[[i]]$VA=list()
	res[[i]]$VA$vec=list(
		 ex30=as.numeric(D[which(probes=="EGFREX30"),]),
		 evii=as.numeric(D[which(probes=="EGFRVII"),]),
		 eviii=as.numeric(D[which(probes=="EGFRVIII"),]),
		 ekd=as.numeric(D[which(probes=="EGFR"),]),
		 pkd=as.numeric(D[which(probes=="PDGFRA"),]),
		 pd89=as.numeric(D[which(probes=="PDGFRAD89"),])
	)
##List the samples with no EGFR PDGFRA signal picked up
	
### Transcriptome Report data
	res[[i]]$TX=list()
	
	v=log10(D[match(tcga$probes,probes),i])
	gi=which(!is.na(v))
	sset=cbind(v[gi],tcga$rawSet[gi,]) #btc one sample combine with 95 samples from tcga to create a new set
	nset=sset-rowMeans(sset) #normalization by mean subtraction 
	cor(nset,tcga$centroids[gi,])->cc #corelation 
	v1=cc[,1]-cc[,2]
	v2=cc[,3]-cc[,4]
	csample=c(3,4,2,1)[which.max(cc[1,])]
	gi=which(!is.na(nset[,1]))
	#print(i)								
	apply(tcga$centroids[gi,],2, function(x) cor.test(nset[gi,1],x)$conf[1])->cl 	#this line causes issues if 'gi' is empty
	apply(tcga$centroids[gi,],2, function(x) cor.test(nset[gi,1],x)$conf[2])->ch
	cd=(ch-cl)/2
	res[[i]]$TX$TCGA=list(csample=csample,v1=v1,v2=v2,cd=cd)
}
time_consumed<-proc.time()-start
Stage2Results=list(res=res,D=D,R=R,DC=DC,Samples=Samples, Stems=Stems, probes=probes)
save(Stage2Results,file=paste(outFolder,"Stage2Results.D",sep=""))

data.frame(Samples=colnames(D),ex30=res[[1]]$VA$vec$ex30,evii=res[[1]]$VA$vec$evii,eviii=res[[1]]$VA$vec$eviii,ekd=res[[1]]$VA$vec$ekd,pkd=res[[1]]$VA$vec$pkd,pd89=res[[1]]$VA$vec$pd89)->VA
write.csv(VA,file="bin/VA.csv")
