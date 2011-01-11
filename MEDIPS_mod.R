source("~/r/statsBMC.R")
library(MEDIPS)

calculateFDR<-function(mp=NULL,method='fdr'){
	for(i in c('p.value.wilcox','p.value.ttest')){
		p<-mp[,i]
		adj<-p.adjust(p,method=method)
		label<-strsplit(i,split=".",fixed=T)
		mp[,paste(method,label[[1]][3],sep=".")]<-adj
	}
	
	mp

	}
	
modifiedSigSelect<-function(frames,fdr=.01,up=2,down=.5,control=T){
	f.1<-frames[frames$rpm_A>0|frames$rpm_B>0,]
	f.2<-f.1[!f.1$fdr.ttest>fdr&!f.1$fdr.wilcox>fdr,]
	f.3<-f.2[!f.2$ratio<up|!f.2$ratio>down,]
	if(control){
		f.4<-f.3[f.3$ratio>1,]
		f.5<-f.4[is.finite(f.4$ratio),]
		}
	else {
		f.4<-f.3[f.3$ratio<1,]
		f.5<-f.4[f.4$ratio>0,]
		}	
	f.5
}
	
makeRatios<-function(frames,ann,...){
	ratio.l<-list()
	ann.f<-read.delim(ann.f,header=F)
	
	###Select significant regions
	for(i in 0:1){
		p<-modifiedSigSelect(frames,...,control=i)
		p.m<-MEDIPS.mergeFrames(p)
		p.a<-MEDIPS.annotate(region=p.m,anno=ann)
		ratio.l[[i+1]]<-table(p.a[,4])
		}
	
	ratio.d<-data.frame(row.names=union(names(ratio.l[[1]]),names(ratio.l[[2]])))
	
	for (i in 1:2){
		m<-matcher(names(ratio.l[[i]]),rownames(ratio.d))
		ratio.d[,i]<-1
		ratio.d[m[,2],i]<-ratio.l[[i]][m[,1]]
		}

	m<-matcher(ann.f[,4],rownames(ratio.d))
	ratio.d<-cbind(ratio.d[m[,2],],ann.f[m[,1],1:3])
	colnames(ratio.d)<-c("down","up","chr","start","end")
	ratio.d$span<-with(ratio.d,(end-start)/1000)
	for(i in c('down','up')){
		ratio.d[,paste(i,"1kb",sep=".")]<-ratio.d[,i]/ratio.d$span		}
	ratio.d<-with(ratio.d,ratio.d[up.1kb>=.1|down.1kb>=.1,])
	ratio.d[,'up.down.log2']<-with(ratio.d,log(up.1kb/down.1kb,2))
	ratio.d<-ratio.d[ratio.d$up.down.log2>=1|ratio.d$up.down.log2<=-1,]
	ratio.d
	}
	
	