# 0.50:  three independant pipeline stages / major revision
#		To Do: TCGA plot colors are wrong
#
#
#

#ls## parser pipeline for Nanostring raw data
#library(gdata)
#library(plyr)
options(stringsAsFactors=FALSE)
rm(list=ls()) # clean the workspace before stage 1 runs; use "list=" rather than "list<-"

folders=dir("workingFolder",full=TRUE) # char vector


# prepare
registry=read.delim("bin/Nanostring_Probe_Registry.txt")
aliasKey=aliasName=c()
for (i in which(registry$Alias!="")) {
	tmp=strsplit(registry$Alias[i],",")[[1]]
	tmp=gsub("[()]","",tmp)
	aliasKey=c(aliasKey,tmp)
	aliasName=c(aliasName,rep(registry$ProbeName[i],length(tmp)))
}
folders=paste(folders,"/",sep="")
#folders=folders[grep("Huse_4.10.14",folders)]						#############folder selection
ref_directory = read.table("dir.in", sep="\t", stringsAsFactors=FALSE)[[1]]
folders=ref_directory

# first pass:  read and parse files
nanoKey=c()
batch=list(); 
i=1
for(folder in folders){
files=as.list(dir(folder,pattern="\\.txt$"))
	for (file in files) {
		filedir=paste(folder,file,sep="")
		cat("Reading in ", filedir,"…\n")
		tmp=read.delim(filedir)
		if (!all(colnames(tmp)[1:3]==c("Code_Class","Name","Accession") )) {
			print("Error: problem with headers for file.  Please check.")
			stop()
		}
		tmp=subset(tmp,Code_Class!="")
		probes=tmp$Name
		class=tmp$Code_Class
		head=scan(filedir,nlines=1,what="character",sep="\t")
		samples=head[-c(1:3)]
		nanoKey=c(nanoKey,samples)
		samplesR=make.names(samples)
		data=array(0,dim=c(nrow(tmp),ncol(tmp)-3))
		for (si in 1:ncol(data)) data[,si]=as.numeric(tmp[,si+3])
		colnames(data)=samplesR

		bi=which(probes %in% aliasKey) 
		if (length(bi)>0) {
			cat("Aliasing: ")
			print(probes[bi])
			probes[bi]=aliasName[match(probes[bi],aliasKey)]
		}
		
		gi=which(class=="Endogenous")  
		bi=which(!(probes[gi] %in% registry$ProbeName ))
		if (length(bi) > 0) {
			print("Warning:  the following probes are not found in the registry and are ignored")
			print(probes[gi][bi])
			data=data[-gi[bi],]
			probes=probes[-gi[bi]]
			class=class[-gi[bi]]
		}
		stem=unlist(lapply(strsplit(folder,"\\/"),function(x) x[[length(x)]]))

		tmpl=list()			
		tmpl$class=class
		tmpl$probes=probes
		tmpl$data=data
		tmpl$file=file
		tmpl$folder=folder
		tmpl$stem=stem
		tmpl$samples=samples
		tmpl$samplesR=samplesR
    }
    
		batch[[i]]=tmpl
		i=i+1
}

#stop()
# second pass:  get all samples and probes 
allSamples=c()
for (i in 1:length(batch)) {
	samples=batch[[i]]$samples
	allSamples=c(allSamples,samples)
}
allProbes=c()
for (i in 1:length(batch)) {
    gi=which(batch[[i]]$class=="Endogenous")
	probes=batch[[i]]$probes[gi]
	allProbes=c(allProbes,probes)
}

vector=allSamples

## check probes: all should pass
uProbes=unique(allProbes)
i=which(!(uProbes %in% registry$ProbeName ))
if (length(i) > 0) {
	print("Error:  the following probes are not found in the registry:")
	print(uProbes[i])
	stop()
}


if (TRUE) {  ## wait to address duplicate samples in Stage2 pipeline
## third pass: 
##  resolve duplicate sample names (uniquify)
	i=which(duplicated(allSamples))
	dupSamples=unique(allSamples[i])
	allSamples=c()
	for (i in 1:length(batch)) {
		folder=batch[[i]]$folder
		stem=gsub("RawData","",folder)
		stem=gsub("/","",stem)
		file=batch[[i]]$file
		nanoSamples=samples=batch[[i]]$samples
		sis=which(samples %in% dupSamples)
		if (length(sis)>0) {
			for (si in sis) {
				nanoSamples[si]=paste(samples[si],"_",stem,sep="")
			}
		}
		batch[[i]]$nanoSamples=nanoSamples
		allSamples=(c(allSamples,nanoSamples))
	}
}


## fourth pass: normalize the data
for (i in 1:length(batch)) {
	d=batch[[i]]$data
	class=batch[[i]]$class
	probes=batch[[i]]$probes
	
	##Background, colMean+colSD will be substracted from each reading
	ci=which(class=="Negative")	
    nm=colMeans(d[ci,],na.rm=T)
	nsd=apply(d[ci,],2,sd,na.rm=T)##for a matrix '1' indicates rows, '2' columns
	nfactor=round(nm+1*nsd)
	nd=t(t(d)-nfactor); nd[nd<0]=0
	
	ci=which(class=="Positive")
	pc=nd[ci,]/rowMeans(nd[ci,],na.rm=T)
	pcf=apply(pc,2,median,na.rm=T)
	pcf[pcf==0]=NA
	pd=t(t(nd)/pcf)
	
	tmp=subset(registry,Control)$ProbeName
	ci=which(probes %in% tmp)
	ec=pd[ci,]/rowMeans(pd[ci,],na.rm=T)
	ecf=apply(ec,2,median,na.rm=T)
	ed=t(t(pd)/ecf)
	
	eset=which(class=="Endogenous")  
	batch[[i]]$normProbes=probes[eset]
	batch[[i]]$norm=ed[eset,]
	
}


## fifth pass: write out results
for (i in 1:length(batch)) {
	folder=batch[[i]]$folder
    Stage1Result=batch[[i]]
    save(Stage1Result,file=paste(folder,"Stage1Result.D",sep=""))
}

