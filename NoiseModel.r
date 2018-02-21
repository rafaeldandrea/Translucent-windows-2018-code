## load package deSolve for ode simulation
library(deSolve)

## types of noise: 
## 1. process noise: complementary or essential resources in the main text model. core, tail, or general in the model in Appendix D
## 2. measurement noise
noisetypes=c('complementary','essential','core','tail','general','measurement')

## scaling factors: scale noise levels in the kernel
scaling.factors=data.frame(noise=c('complementary','essential','core','tail','general','measurement'),scaling=c(30,1,8,250,2,.4))

## number of species 
S=200 								

## species niche/trait values
set.seed(0)
trait=sort(runif(S))
d=as.matrix(dist(trait))
d[d>max(d)/2]=max(d)-d[d>max(d)/2]

## kernel width
w=.1

## noise-free competition kernel
A0=exp(-(d/w)^4); A0=A0/mean(A0)*.2

## Function to generate the competition kernel.
## Inputs: species x-niche/trait values, type of noise, degree of noise, simulation number (run).
## Outputs: A = noisy kernel in the case of process noise, noise-free kernel in the case of measurement noise
## (cont) A0 = noise-free kernel in the case of process noise, noisy kernel in the case of measurement noise; 
## (cont) d = distances on x-niche/trait axis (on a circular axis);
## (cont) prox = proxy trait values in the case of measurement noise; traity = y-niche/trait values;
## (cont) spear = Spearman's rho, correlation between kernel and distances on x-axis;
## (cont) pi = P(Atail > Amean) - P(Acore > Amean), quantifies core-tail structure in the kernel
Kernel=function(noiselevel,noisetype,traitvalues=trait,run=0){
	stopifnot(noisetype%in%noisetypes)
	
	## tie random number generator to run
	set.seed(run)
	
	## number of species
	S=length(traitvalues)
	
	## scaled noise level
	a0=noiselevel*with(scaling.factors,scaling[match(noisetype,noise)])
	
	## distances on x-niche/trait axis (axis is circular)
	d=as.matrix(dist(traitvalues)); d[d>max(d)/2]=max(d)-d[d>max(d)/2]
	
	## A0 = noise-free kernel in the case of process noise, noisy kernel in the case of measurement noise
	A0=exp(-(d/w)^4); A0=A0/mean(A0)*.2
	
	## threshold distance separating core from tail
	dcoretail=d[which.min(abs(A0-mean(A0)))]	
	
	## traity = y-niche/trait values; will contribute to kernel in proportion to noise level
	traity=NULL
	
	## A = noisy kernel in the case of process noise, noise-free kernel in the case of measurement noise
	if(noisetype%in%c('complementary','essential')){ traity=runif(S); dy=as.matrix(dist(traity)); dy[dy>max(dy)/2]=max(dy)-dy[dy>max(dy)/2]}
	if(noisetype=='complementary'){ wyinverse=a0; A=exp(-sqrt((d/w)^2+(dy*wyinverse)^2)^4)}		## from 1 to 10
	if(noisetype=='essential'){ wy=a0; A=pmax(exp(-(d/w)^4),exp(-(dy/wy)^4))}				## from .01 to .15
	
	if(noisetype%in%c('core','tail','general')){ 
		sigma = if(noisetype=='core') a0*exp(-(d/.1)^4) else if(noisetype=='tail') a0*(exp(d^4)-1) else if(noisetype=='general') a0
		A=matrix(pmax(0,rnorm(S*S,mean=exp(-(d/w)^4),sd=sigma)),S,S)
		while(TRUE){ i=which(diag(A)==0); if(length(i)==0) break; diag(A)[i]=pmax(0,rnorm(length(i),mean=1,sd=a0))}
	}
	
	## prox = proxy trait values in the case of measurement noise
	prox=NULL
	if(noisetype=='measurement'){ 
		prox=rnorm(S,mean=traitvalues,sd=a0); d=as.matrix(dist(prox)); d[d>max(d)/2]=max(d)-d[d>max(d)/2]		
		A=A0; A0=exp(-(d/w)^4); A0=A0/mean(A0)*.2
	}
	
	## normalize A
	A=A/mean(A)*.2
	
	## spear = Spearman's rho; pi = P(Atail > Amean) - P(Acore > Amean)
	spear=cor(as.numeric(d),as.numeric(A),method='spearman')											## -1 (perfect) to 0 (null)
	pi=sum(A[d>=dcoretail]>mean(A))/sum(d>=dcoretail)-(sum(A[d<=dcoretail]>mean(A))/sum(d<=dcoretail))	## -1 (perfect) to 0 (null)
	
	return(list(
		A=A,
		A0=A0,
		d=d,
		prox=prox,
		traity=traity,
		spear=spear,
		pi=pi
	))
}

Run=function(traits,compkernel,immig,maxtime=1e4,plot=1,logplot=1){	
	## initial abundances
	x0=rep(1e-3,nrow(compkernel))+runif(nrow(compkernel),min=0,max=1e-3/100)	
	
	## time series
	times=seq(0,maxtime,by=.1)
	
	## Lotka-Volterra dynamic equations
	LV=function(times,N,parms) list(pmax(N,0)*(1-as.numeric(parms$Kernel%*%pmax(N,0)))+parms$m)
	
	## Simulate model
	x=ode(y=x0,times=times,func=LV,parms=list(Kernel=compkernel,m=immig),method='lsoda'); n=x[nrow(x),-1]; n[n<1e-5]=0
	
	## Calculate richness at each time point
	s=apply(x[,-1],1,function(v) sum(v>1e-5))
	
	## Plot results
	if(plot) if(!logplot) plot(trait,n,t='h',las=1) else plot(trait[n>0],log(1+n[n>0]/min(n[n>0])),t='h',las=1)
	
	return(list(data=data.frame(trait=traits,abundance=n),richness=data.frame(time=x[,1],S=s)))	
}
