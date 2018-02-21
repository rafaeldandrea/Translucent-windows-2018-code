KmeansGap=function(dat,nozeros=FALSE,multiD=0,numnulls=100,verbose=1,mink=1,maxk=NULL,nulltype='shuffle',nstartingpoints=100,YanYe=0,plot=0,line90=1){
	set.seed(0)
	if(class(dat)!='data.frame') stop('Valid input to KmeansGap must be a data frame')
	if(!multiD & !'trait'%in%names(dat)) stop('Input data frame for the 1d case must have a "trait" column')
	if(!'N'%in%names(dat)) stop('Input data frame must have a "N" column')
	
	## prepare data for analysis
	dat=dat[apply(dat,1,function(v) all(!is.na(v))),] 
	if(nozeros) dat=dat[dat$N>0,]
	dat$N=dat$N/min(1,min(dat$N[dat$N>0])) 
	if(!multiD) dat=dat[order(dat$trait),]
	
	traitdata=if(multiD) subset(dat,select=-N) else data.frame(trait=dat$trait)
	
	if(is.null(maxk)) maxk=round(min(20,sum(dat$N>0)/2))
	kvec=mink:maxk
	
	library(plyr)
	rep.row=function(r,a) colwise(function(x) rep(x,a))(r)
	kms=function(mat,n,k){
		if(ncol(mat)>1) com=rep.row(r=mat,a=n) else com=rep(mat$trait,n)
		if(k==1){ 
			z=kmeans(com,centers=1,iter.max=100)
		}else{
			nullcenters=apply(unique(mat),2,function(v) sample(seq(min(v),max(v),l=k)))
			z=try(kmeans(com,centers=nullcenters,iter.max=100),silent=TRUE)
			if(!is.null(attr(z,'condition'))) z=kmeans(com,centers=k,nstart=nstartingpoints,iter.max=100)
		}
		return(z)
	}
	
	Dispersion=function(dtf,weight){
		if(weight==0) return(dtf$tot.withinss)							## Unweighted sum of dist^2
		if(weight==1) return(with(dtf,sum(withinss/(2*size))))			## Tibshirani et al 2001
		if(weight==2) return(with(dtf,sum(withinss/(2*size*(size-1)))))	## Correction from Yan & Ye 2007
	}
	
	Wk=sapply(kvec,function(k){ 	## Wk = within-cluster sum of squares (ie squared distances); D in Tibshirani et al 2001
		mod=kms(mat=traitdata,n=dat$N,k)
		Dispersion(dtf=mod,weight=YanYe)
	})
	
	nullWk=sapply(kvec,function(k){
		if(verbose) cat(".",if(k%%10==0 | k==maxk) paste(k,"\n"))
		if(nulltype=='shuffle'){ 
			ndat=sapply(seq(numnulls),function(run) sample(dat$N))
			return(apply(ndat,2,function(N){ 
				mod=kms(mat=traitdata,n=N,k)
				Dispersion(dtf=mod,weight=YanYe)
			}))
		}else if(nulltype=='spline'){
			ndat=rnorm(length(dat$N),mean=dat$N,sd=sd(dat$N)); ndat[ndat<0]=0
			return(apply(ndat,2,function(N){ 
				mod=kms(mat=traitdata,n=N,k)
				Dispersion(dtf=mod,weight=YanYe)
			}))
		}else{
			return(apply(ndat,2,function(N){ 
				tmp=apply(as.matrix(traitdata),2,function(v) runif(length(v),min=min(v),max=max(v)))
				mod=kms(mat=as.data.frame(tmp),n=N,k)
				Dispersion(dtf=mod,weight=YanYe)
			}))
		}
	})
	
	logWk=log(Wk)										## vector length(kvec) 
	lognullWk=log(nullWk)								## matrix numnulls x length(kvec)
	ElogWk=apply(lognullWk,2,mean,na.rm=TRUE)			## vector length(kvec)
	Gapk=ElogWk-logWk									## vector length(kvec)
	nullGapk=t(ElogWk-t(lognullWk))						## matrix numnulls x length(kvec)
	
	maxgap=max(Gapk,na.rm=TRUE)							## scalar
	maxnullgap=apply(nullGapk,1,max,na.rm=TRUE)			## vector numnulls
	kmaxnullgap=apply(nullGapk,1,which.max)				## vector numnulls
	kmaxnullquant=merge(data.frame(kmax=kvec),ddply(data.frame(kmax=kmaxnullgap,max=maxnullgap),.(kmax),function(v) quantile(v$max,.95)),by='kmax',all.x=TRUE)
	kmaxnullquant90=merge(data.frame(kmax=kvec),ddply(data.frame(kmax=kmaxnullgap,max=maxnullgap),.(kmax),function(v) quantile(v$max,.9)),by='kmax',all.x=TRUE)
	khat=kvec[which.max(Gapk)]							## scalar
	z.score=(maxgap-mean(maxnullgap))/sd(maxnullgap)	## scalar
	p.value=sum(maxnullgap>maxgap)/length(maxnullgap)	## scalar
	
	mngap=apply(nullGapk,2,mean,na.rm=TRUE); 			## mean  of the null gaps --- vector length(kvec)
	sdgap=apply(nullGapk,2,sd,na.rm=TRUE); 				## stdev of the null gaps --- vector length(kvec)
	nullquant=apply(nullGapk,2,function(vec) as.numeric(quantile(vec,.95,na.rm=TRUE)))	## 95% quantile of the null gaps --- vector length(kvec)
	nullquant90=apply(nullGapk,2,function(vec) as.numeric(quantile(vec,.9,na.rm=TRUE)))	## 90% quantile of the null gaps --- vector length(kvec)
	maxnullquant=quantile(maxnullgap,.95)				## 95% quantile of the max null gaps --- scalar
	maxnullquant90=quantile(maxnullgap,.9)				## 90% quantile of the max null gaps --- scalar
	
	if(plot){		
		plot(0,t='n',xlim=c(min(kvec),max(kvec)),ylim=c(0,max(Gapk,maxnullquant)),las=1,xlab='No. clusters',ylab='Gap')
		polygon(x=c(kvec,rev(kvec)),y=c(mngap,rev(nullquant)),col='grey90',border=NA)
		polygon(x=c(kvec,rev(kvec)),y=c(mngap,rev(nullquant90)),col='grey80',border=NA)
		lines(kvec,rep(maxnullquant,length(kvec)),col=2)
		if(line90) lines(kvec,rep(maxnullquant90,length(kvec)),col=2,lty=2)
		points(kvec,Gapk,t='o',pch=20,lwd=2)
		polygon(x=c(kvec,rev(kvec)),y=c(mngap,rep(-maxnullquant,length(kvec))),col='white',border=NA)
		box()
	}
	
	if(YanYe==0) disp='D'; if(YanYe==1) disp='D/2n'; if(YanYe==2) disp='D/(2n(n-1))'
	
	return(list(
		data=data.frame(
			k=kvec,
			gap=Gapk,
			Egap=mngap,
			sdgap=sdgap,
			nullquant=nullquant,
			nullquant90=nullquant90,
			kmaxnullquant=kmaxnullquant[,2],
			kmaxnullquant90=kmaxnullquant90[,2],
			logWk=logWk,
			ElogWk=ElogWk
		),
		khat=khat,
		maxgap=maxgap,
		maxnullquant=maxnullquant,
		maxnullquant90=maxnullquant90,
		mink=mink,
		maxk=maxk,
		z.score=z.score,
		p.value=p.value,
		dispersion=disp
	))
}
