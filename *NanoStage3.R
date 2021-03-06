# disable all the pdf() in R source code
# generate .Rnw file
# run Sweave(.Rnw) inside of R to generate .tex file
# or run texi2dvi(.tex, pdf=TRUE) to generate .pdf
# need bash to call pdflatex on .tex file to generate .pdf
# stage 3
# v3_4, create .Rnw file rather than .tex file.
# v3_7, ".." in file names and sample names are taken care of.
# v3_8, zscored the locus expression data by creating res[[i]]$LL$exp[[ll]]$zs twisted in v3_9
# v3_9, run NanoStage2_1 (z-scored LL raw data)
# Ambitious version: linking NanoString, Sequenom and IHC using the BTC_ID key
# v3_10, Make Locus section removable, Nanostring_Probe_Registry.txt is revised
# v3_10, fixed the variant results
# v3_10, also make color plots for variants
# v3_11, v3_10 is broken. 1) C-term change color to black from grey 2) pd89 QC value cut off is log(10) >=2; and we only record two categories 3) ++ 2 red; + 8 grey; - 1 black
# v3_12, separate oncogenes and tumor suppressor, didn't work well with .Rnw file
# v3_13, generate color.D, the report stage is transferred to stage4
# v3_14, build new object CGE cancer Gene Expression table, color code all the variant threshold, add Two classifications
# v3_16, fix the BTC cohort plotting

library(gplots)
library(shape)
library(xtable)
library(tools) # for function texi2dvi(.tex,pdf=TRUE)
library(plyr) # to use this function: rbind.fill()
library(R2HTML)
library(knitr)
options(stringsAsFactors=F)


nano_caisis=list()
cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

color1=color2=color3=color4=c()
color=as.data.frame(color1,color2,color3,color4)


subclassVec=c("Neural", "Classical", "Mesenchymal","Proneural")


registry=read.delim("bin/Nanostring_Probe_Registry.txt")
og=registry$ProbeName[registry$CancerGeneType=="OG"]
ts=registry$ProbeName[registry$CancerGeneType=="TS"]
registry$LocusID[registry$LocusType=="AMP"]->AmpLL
AmpLL=unique(AmpLL)
registry$LocusID[registry$LocusType=="DEL"]->DelLL
DelLL=unique(DelLL)
#===================For the Final Pipeline===================================#
outFolder="Reports/"
inFolder="workingFolder/"
load("bin/TCGA_result.D"); tcga=res
load("Stage2Results.D")
file.remove("Stage2Results.D")
res=Stage2Results$res
D=Stage2Results$D
R=Stage2Results$R
DC=Stage2Results$DC
probes=Stage2Results$probes
Stems=Stage2Results$Stems
Samples=Stage2Results$Samples
D_S=D[,grep("S_",colnames(DC))]
D_F=D[,grep("F_",colnames(DC))]
R_S=R[,grep("S_",colnames(R))]
R_F=R[,grep("F-",colnames(R))]
DC_F=D[,grep("F_",colnames(DC))]
DC_S=D[,grep("S_",colnames(DC))]
D_CC=cbind(D_S,D_F)
R_CC=cbind(R_S,R_F)
DC_CC=cbind(DC_S,DC_F)
#D=D_CC
#DC=DC_CC
#cor(DC[match(tcga$probes,probes),],tcga$centroids,use="pairwise")->BTCCentroidCor
#as.numeric(apply(BTCCentroidCor,1, function(x) c(3,4,2,1)[which.max(x)]))-> BTCCentroidCol
c(1,2,3,4)[match(tcga$class,c("N","C","M","P"))]-> TCGACentroidCol

#cor(DC[match(tcga$probes,probes),],tcga$centroids,use="pairwise")->BTCCentroidCor_MC
#as.numeric(apply(BTCCentroidCor_MC,1, function(x) c(3,4,2,1)[which.max(x)]))-> BTCCentroidCol_MC
####################################################
Ccol=list()
V1=list()
V2=list()
Csamples=list()
Cd=list()
for(i in 1:ncol(D)){
	if(is.null(res[[i]]$TX$TCGA$v1)){ 
		Csamples[[i]]=0
		V1[[i]]=0
		V2[[i]]=0
		next
	}
	ccol=TCGACentroidCol[1]
	v1=res[[i]]$TX$TCGA$v1[1]
	v2=res[[i]]$TX$TCGA$v2[1]
	csample=res[[i]]$TX$TCGA$csample
	cd=res[[i]]$TX$TCGA$cd	
	Ccol[[i]]=ccol
	V1[[i]]=v1
	V2[[i]]=v2
	Csamples[[i]]=csample
	Cd[[i]]=cd
}
####################################################

nanoKey=c()
RES=list()
index=0
CGE=rep("NA",16)
ClinicTrial<-c(colnames(D)[grep("F_",colnames(D))],colnames(D)[grep("S_",colnames(D))])
CTpos<-match(ClinicTrial,colnames(D))
subHTMLDir<-"/Reports"
for (i in 1:ncol(D)) {
    #print[i]
    corename=unlist(strsplit(colnames(DC)[i],"_"))[1]
    corename=gsub("X","",corename)
    corename=unlist(strsplit(corename,"[.]"))[1]
    nanoKey=c(nanoKey,corename)
	print(colnames(DC)[i])
	if (all(is.na(D[,i]))) {
		cgexpp=""
		CGE=rbind(CGE,cgexpp)
		next
	}
    cat("Test0")
	
    ## QC Report plot	
	cm=res[[i]]$QC$cm
  	if(is.null(cm)) next
	    
    ## Cancer Genes Report plot. This part is not consistent with .Rnw file.
	gset=names(res[[i]]$CG$exp)
    cgexpp=c()#log values for the whole gset
#cgCF=c(3.7,0.6,2.44,4.1,3.28,4,3.14,2,4.66,2.99,3.99,3.66,3.94,3.19,2.76,2.87)
	if(FALSE){
	    for (j in 1:length(gset)) {
			exp=res[[i]]$CG$exp[[j]]
			v=na.omit(exp[-i])
			
			sv<-sort(v)
			
				if(j==3) {
					tmp<-sv[round(length(sv)*0.05,0)]
					mexp=log10(tmp)# bottom 5%
				}
				else if(j==5) {mexp=log10(mean(exp,na.rm=T)+2*sd(exp,na.rm=T))}
				else if(j==7){
					tmp<-sv[round(length(sv)*0.9,0)]
					mexp=log10(tmp)# top 10%
					}
				else if(j==8) mexp=log10(mean(exp,na.rm=T)-sd(exp,na.rm=T))
				else if(j==9) mexp=log10(mean(exp,na.rm=T)+2*sd(exp,na.rm=T))
				else if(j==10) mexp=log10(mean(exp,na.rm=T)-sd(exp,na.rm=T))
				else if(j==11){
					tmp<-sv[round(length(sv)*0.9,0)]
					mexp=log10(tmp)
				} 
				else if(j==12){
					tmp<-sv[round(length(sv)*0.96,0)]
					mexp=log10(tmp)
				} 
				else if(j==13) mexp=log10(mean(exp,na.rm=T)+sd(exp,na.rm=T))
				else if(j==14){
					tmp<-sv[round(length(sv)*0.97,0)]
					mexp=log10(tmp)
				} # top 4%
				else if(j==15){
					tmp<-sv[round(length(sv)*0.04,0)]
					mexp=log10(tmp)
				} # top 4%
				else if(j==16){
					tmp<-sv[round(length(sv)*0.04,0)]
					mexp=log10(tmp)
				} # top 4%
				else mexp=log10(mean(exp,na.rm=T))
				#print(j)
				#print(mexp)
				
				cgCF=c(cgCF,mexp) 
				
			vv=exp[i]
			xx=log10(vv)
			#pp=round(sum(vv>v,na.rm=T)/length(v)*100,1)
			#pl=round(sum(vv<v,na.rm=T)/length(v)*100,1)
			#tmp=par()$xaxp; xwidth=tmp[2]-tmp[1]; xch=xwidth/4
		    cgexpp=c(cgexpp,xx)
		}
	    cgCF[1]=3.7
	    cgCF[2]=0.6 # also another cut-off value 2.1
	    cgCF[4]=4.1
	    cgCF[6]=4 # around top 10%
	}

cgCF=c(3.6,0.5,2.44,4.0,3.28,4,3.14,2,4.75,2.99,3.99,3.66,3.94,3.19,2.76,2.87)# modified in Oct 2012
	for (j in 1:length(gset)) {
		exp=res[[i]]$CG$exp[[j]]
		v=na.omit(exp[-i])
		vv=exp[i]
		xx=log10(vv)
		cgexpp=c(cgexpp,xx)
	}
    res[[i]]$CGE$self=cgexpp
	CGE=rbind(CGE,cgexpp)
# Cancer Genes Cut Off
    cgGmsg=""
    cgLmsg=""
    
    for(jj in 1:length(gset)){
        if(is.na(cgCF[jj])) next
		if((gset[jj] %in% og) && (cgexpp[jj]>cgCF[jj])){
           	     cgGmsg=paste(cgGmsg,gset[jj],sep=" ")
         }else if((gset[jj] %in% ts) && (cgexpp[jj]<cgCF[jj]))
                cgLmsg=paste(cgLmsg,gset[jj],sep=" ")
         else next
    }
    cgexppog=cgexpp[match(og,gset)]
	cgCFog=cgCF[match(og,gset)]
	cgexppts=cgexpp[match(ts,gset)]
	cgCFts=cgCF[match(ts,gset)]
    ########## Variant############### 
    ########## Variant############### 
	##only show cut-off values##
	eviii_cf=paste("[-]=0~0.01","[+]=0.01~0.1","[++]>0.1")
	evii_cf=paste("[-]=0~0.01","[+]=0.01~0.1","[++]>0.1")
	ecterm_cf=paste("[-]>0.1","[+]=0.01~0.1","[++]=0~0.1")
	pd89_cf=paste("[-]=0~0.01","[+]>0.1","[++]>0.1")
	logicVec=list()
	logicVec=rbind(logicVec,"eviii_cf"=eviii_cf)
	logicVec=rbind(logicVec,"evii_cf"=evii_cf)
	logicVec=rbind(logicVec,"ec-term_cf"=ecterm_cf)
	logicVec=rbind(logicVec,"pd89_cf"=pd89_cf)
	logicVec=as.data.frame(logicVec)
	VcfT<-xtable(logicVec)
	colnames(VcfT)<-c("Cut-off Range")
	################

	ekd=res[[i]]$VA$vec$ekd
	eviii=res[[i]]$VA$vec$eviii
	evii=res[[i]]$VA$vec$evii
	ex30=res[[i]]$VA$vec$ex30
	pd89=res[[i]]$VA$vec$pd89
	pkd=res[[i]]$VA$vec$pkd
    
	eviiirf<-log10(eviii)/log10(ekd)
	eviirf<-log10(evii)
	ex30rf<-log10(ex30)/log10(ekd)
	pd89rf<-log10(pd89)/log10(pkd)
	variantName=c("EGFRVIII","EGFRVII","EGFRC","PDGFRd89")
	variantVec=matrix(c(0,0.01,0.1,0,0.01,0.1,0.1,0.01,0,0,0.01,0.1),nrow=4,ncol=3,byrow=T)
	dataVec=c("viii"=eviii[i]/ekd[i],"vii"=evii[i]/ekd[i],"c-term"=ex30[i]/ekd[i],"pd89"=pd89[i]/pkd[i])
    varRes=c()
    for(va in 1:length(dataVec)){
    	if(is.na(dataVec[va])) next
    	if(va==1){
    		if(log10(ekd[i])<2 || log10(eviii[i])<2){
    			varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
    		}else{
	    		if(dataVec[va]>variantVec[va,3]){
		    		varRes[va]=paste(variantName[va],"[++]",sep="")
		            color[i,va]=2
		        }else if(dataVec[va]>variantVec[va,2]){
		            varRes[va]=paste(variantName[va],"[+]",sep="")
		            color[i,va]=7
		        }else{
		            varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
		        }
    		}
    	}
    	if(va==2){
    		if(log10(ekd[i])<2 || log10(evii[i])<1){
    			varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
    		}else{
	    		if(dataVec[va]>variantVec[va,3]){
		    		varRes[va]=paste(variantName[va],"[++]",sep="")
		            color[i,va]=2
		        }else if(dataVec[va]>variantVec[va,2]){
		            varRes[va]=paste(variantName[va],"[+]",sep="")
		            color[i,va]=7
		        }else{
		            varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
		        }
    		}
    	}
    	
    	if(va==3){
    		if(log10(ekd[i])<2 || log10(ex30[i])< 2){
    			varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
    		}else{
		    	if(dataVec[va]>variantVec[va,1]){
		    		varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
		        }else if(dataVec[va]>variantVec[va,2]){
		            varRes[va]=paste(variantName[va],"[+]",sep="")
		            color[i,va]=7
		        }else{
		            varRes[va]=paste(variantName[va],"[++]",sep="")
		       		color[i,va]=2
		        }
            }
        }
    	
    	if(va==4){ 
    		if(log10(pkd[i])<2 || log10(pd89[i])<1){
    			varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
    		}else{
	    		if(dataVec[va]>variantVec[va,3]){
		    		#varRes[va]=paste(variantName[va],"[++]",sep="")
		    		#color[i,va]=2
		    		varRes[va]=paste(variantName[va],"[+]",sep="")
		            color[i,va]=7
		    	}
		        #else if(dataVec[va]>variantVec[va,2]){
		        #   varRes[va]=paste(variantName[va],"[+]",sep="")
		        #    color[i,va]=8}
		        else{
		            varRes[va]=paste(variantName[va],"[-]",sep="")
		            color[i,va]=1
		        }
	        }
        }
        print("dataVec")
        print(dataVec[va])
	    print(varRes[va])
    }
    
    message=paste(varRes,collapse=" ")  
    ######### Transcriptome Report plot #########
	par(mar=c(6,6,4,2), mfrow=c(1,1))  ###  plot TCGA clouds
	ccol=TCGACentroidCol
	v1=res[[i]]$TX$TCGA$v1
	v2=res[[i]]$TX$TCGA$v2
	csample=res[[i]]$TX$TCGA$csample
	cd=res[[i]]$TX$TCGA$cd	
    #cc=BTCCentroidCor
	#csamples=BTCCentroidCol
	#v1=cc[,1]-cc[,2]
	#v2=cc[,3]-cc[,4]
	gi=which(!is.na(DC[,i]))
	
	apply(tcga$centroids[gi,],2, function(x) cor.test(DC[gi,i],x)$conf[1])->cl
	apply(tcga$centroids[gi,],2, function(x) cor.test(DC[gi,i],x)$conf[2])->ch
	cat("testtest")
	cd=(ch-cl)/2
    cat("Test2")
    ######################################
   	summary=list()
    ReceivedDate=unlist(gsub("BTC_","",Stems[i]))
    summary$meta=c("Summary","Specimen: "=colnames(DC)[i],"Date: "=ReceivedDate)
    summary$cg=list()
    resultMeta=list()
    resultMeta=rbind(resultMeta,"Sample"=colnames(DC)[i])
    resultMeta=resultMeta[-1,]
    fn=gsub("[..]","-",colnames(DC)[i])
    fn=gsub("[.]","-",fn)
    ReceivedDate=gsub("_","/",ReceivedDate)
    resultMeta=rbind(resultMeta,"Processed Date"=ReceivedDate)
    resultMeta=rbind(resultMeta,"QC result "=round(cm[colnames(D)[i]],2))
	resultMeta2=resultMeta
    resultMeta2=rbind(resultMeta2,"Cancer Gene Gain"=cgGmsg)
    resultMeta2=rbind(resultMeta2,"Cancer Gene Loss"=cgLmsg)
    resultMeta2=rbind(resultMeta2,"EGFR, PDGFRa Variants "=message)
	resultMeta2=rbind(resultMeta2,"TCGACentroid Subclass"=subclassVec[Csamples[[i]]])
    #resultMeta2=rbind(resultMeta2,"BTCCentroid Subclass "=subclassVec[csamples[i]])
	#if(subclassVec[Csamples[[i]]]==subclassVec[csamples[i]]) Agree="Agree"
    #else Agree="Not Agree"
    #resultMeta2=rbind(resultMeta2,"Two Classifications"=Agree)
    print(resultMeta2)
	cat("Test3") 
    msg<-data.frame(cgGmsg,cgLmsg)
       
    
	nano_resString=paste("NSG2012:","[",subclassVec[Csamples[[i]]],"] ",message)
    tmp=c(colnames(DC)[i],nano_resString)
    nano_caisis=cbind(nano_caisis,tmp)
	index=index+1
	RES[[i]]=resultMeta2
}
save(color,file="bin/color.D")
#Sweave("cancerGene.Rnw")
#texi2dvi("cancerGene.tex",pdf=T)
#Sweave("variant.Rnw")
#texi2dvi("variant.tex",pdf=T)
#Sweave("variant2.Rnw")
#texi2dvi("variant2t.tex",pdf=T)
#Sweave("TCGA_Plot.Rnw")
#texi2dvi("TCGA_Plot.tex",pdf=T)
#Sweave("BTC_cohortPlot.Rnw")
#texi2dvi("BTC_cohortPlot.tex",pdf=T)
#file.remove(dir(pattern=".tex"))
#file.remove(dir(pattern=".log"))
#file.remove(dir(pattern=".aux"))
#stop()
##########################################
##########################################
##############Data Organization###########
##########################################
##########################################
load("bin/Ted.D")
gsub("-","\\.",TedTable[,2])->Ted_surg
Res=RES[[1]]
for(j in 2:length(RES)){
    Res=cbind(Res,RES[[j]])
}
as.data.frame(Res)->Res
Res<-t(Res)
mid_names=c()
for(i in 1:nrow(Res)){
    mid_names[i]=unlist(strsplit(rownames(Res)[i],"_"))[1]
}
last_names=c()
for(i in 1:nrow(Res)){
    if(mid_names[i] %in% Ted_surg){
        last_names[i]=TedTable[match(mid_names[i],Ted_surg),1]
    }
    else
        last_names[i]=mid_names[i]
}
cbind(Res,mid_names)->Res
cbind(Res,last_names)->Res
#D_BTC<-TedTable[match(Ted_surg,colnames(D)),1]
#rbind(Res,D_BTC)->Res
#colnames(Res)<-D_BTC
write.csv(Res,file="bin/res.csv")
#stop()

load("bin/color.D")
GenerateCombinedReport=F
GeneratePDFReport=F
GenerateHTMLReport=T
LLReport=F
subHTMLDir<-paste(getwd(),"Reports",sep="/")
#InterestList=c(333,98,277,170,266,203,88)
write.csv(t(Res),file="bin/Res.csv")
for (i in 1:ncol(D)) {
	cm=res[[i]]$QC$cm
    CGE_self=res[[i]]$CGE$self
	if(is.null(cm)) 
	next
	BT<-xtable(RES[[i]],Caption="Summary",align=c("l","c"))
	fn=gsub("[..]","-",colnames(DC)[i])
	fn=gsub("[.]","-",fn)
	RnwFilename=paste(fn,".Rhtml",sep="")
	if(GenerateHTMLReport){
		figPath=paste("figure/",fn,sep="")
		file.copy("bin/nanoString_template.Rhtml",RnwFilename,overwrite=T)
		#Sweave(RnwFilename,driver=RweaveHTML)
		knit(RnwFilename)
		#file.rename(dir(pattern=fn),paste(subHTMLDir,dir(pattern=fn),fn,sep="/"))
	}
	nano_resString=paste("NSG2012:","[",subclassVec[Csamples[[i]]],"] ",message)
	tmp=c(colnames(DC)[i],nano_resString)
	nano_caisis=cbind(nano_caisis,tmp)
}
system("mv  --backup=numbered figure out")
system("mv *.html out")
system("rm *.Rhtml")
