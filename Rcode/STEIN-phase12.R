library(magrittr)
library(Iso)

get.oc <- function(target, pE.true,pT.true,  psi1, psi2, ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=10)
{
	peestimate<-function(yE,n){
	ndose<-length(yE)
	lik<-rep(0,ndose)
	pe<-(yE+0.05)/(n+0.1)
	p.e<-matrix(NA,ndose,ndose)
	for (i in 1:ndose){
		if (i==1) {x<-seq(ndose,1,by=-1)} else {x<-c(1:(i-1),seq(ndose,i))}
		#x<-x
		p.e[i,]<-ufit(pe,lmode=i,x=x,w=n+0.5)[[2]]
		lik[i]<-prod(dbinom(yE,n,p.e[i,]))		
	}
	lik<-lik/sum(lik)
	pe<-t(p.e)%*%lik+0.01*seq(1,ndose)
	return(pe)}

	## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
	get.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95)
	{
 
# if the user does not provide p.saf and p.tox, use the default values
		if(is.na(p.saf)) p.saf=0.75*target;
		if(is.na(p.tox)) p.tox=1.25*target;

### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
		npts = ncohort*cohortsize;
		ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
		for(n in (1:ncohort)*cohortsize)
		{
			error.min=3;
			for(m1 in 0:(n-1))
			{
				for(m2 in (m1+1):n)
				{
 
						error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
						error2 = 1-pbinom(m1, n, p.saf);
						error3 = pbinom(m2-1, n, p.tox);
					
					error=error1+error2+error3;
					if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
				}
			}
			ntrt = c(ntrt, n);
			b.e = c(b.e, cutoff1);
			b.d = c(b.d, cutoff2);
			
			elimineed=0; # indicating whether elimination is needed
			if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
			else
			{
				for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
				{
					if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
				}
				if(elimineed==1) { elim = c(elim, ntox); }
				else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
			}
		}
		for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
		boundaries = rbind(ntrt, elim, b.d, b.e);
		rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
								 "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
		colnames(boundaries) = rep("", ncohort);
		
		return(boundaries);
	}
	
targetT<-0.3
set.seed(30);
ndose=length(pE.true)	
npts = ncohort*cohortsize;
YT=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
YE=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
dselect = rep(0, ntrial); # store the selected dose level
sel=rep(0,ndose);
pts=rep(0,ndose);
dlt=rep(0,ndose);
eff=rep(0,ndose);
ntox=0;
neff=0;
acr=0;
exc=0;
temp=get.boundary(targetT, ncohort, cohortsize,cutoff.eli=0.95) 	
b.e=temp[4,];   # escalation boundary
b.d=temp[3,];   # deescalation boundary
b.elim=temp[2,];  # elimination boundary
psi<-log((1-psi1)/(1-psi2))/log(psi2*(1-psi1)/psi1/(1-psi2))
#print(cutoff.e)
	################## simulate trials ###################
	for(trial in 1:ntrial)
	{

		yT<-yE<-rep(0, ndose);    ## the number of DLT at each dose level
		n<-rep(0, ndose);    ## the number of patients treated at each dose level
		earlystop=0;         ## indiate whether the trial terminates early
		d=startdose;         ## starting dose level
		elimi = rep(0, ndose);  ## indicate whether doses are eliminated
		elimiE=  rep(0,ndose);
		
		for(i in 1:ncohort)  
		{  			
			### generate toxicity outcome
			wT = sum(runif(cohortsize)<pT.true[d])
			yT[d] = yT[d] + wT;
			wE = sum(runif(cohortsize)<pE.true[d])
			yE[d] = yE[d] + wE;
			n[d] = n[d] + cohortsize;
			nc = n[d]/cohortsize;

			if(!is.na(b.elim[nc]))
			{
				if(yT[d]>=b.elim[nc]) 
				{      
					elimi[d:ndose]=1;
					if(d==1) {earlystop=1; break;} 
				}
				
			}
			if(n[d]>=3 && pbeta(0.2,yE[d]+1,n[d]-yE[d]+1)>0.95) {elimiE[d]=1;}
			
 

			if (yT[d]>=b.d[nc] && d!=1) {d_opt=max(which(elimi[1:(d-1)]==0))
			} else if (yT[d]>=b.d[nc] && d==1) {d_opt=d
			} else{
			       if (yE[d]/n[d]>psi) {d_opt=d
			} else if (yT[d]> b.e[nc]) {
			posH<-c(pbeta(psi,yE[max(1,d-1)]+1,n[max(1,d-1)]-yE[max(1,d-1)]+1),
			pbeta(psi,yE[d]+1,n[d]-yE[d]+1))-c(0,0.001)
            d_opt<-max(d-2+which.min(posH),1);
			} else if (yT[d]<=b.e[nc] && d==1) {
			       if (elimi[d+1]==0) {
			posH<-c(pbeta(psi,yE[d]+1,n[d]-yE[d]+1),
			(n[d+1]>0)*pbeta(psi,yE[d+1]+1,n[d+1]-yE[d+1]+1))-c(0,0.001)	
			d_opt=d-1+which.min(posH)} else {d_opt=d}
			} else if (yT[d]<=b.e[nc] && d==ndose) {
			posH<-c(pbeta(psi,yE[d-1]+1,n[d-1]-yE[d-1]+1),
			pbeta(psi,yE[d]+1,n[d]-yE[d]+1))-c(0,0.001)
			d_opt=d-2+which.min(posH)
		    } else if (yT[d]<=b.e[nc] && d!=1 && d!=ndose) {
		    	       if (elimi[d+1]==0){
		    	posH<-c(pbeta(psi,yE[d-1]+1,n[d-1]-yE[d-1]+1),pbeta(psi,yE[d]+1,n[d]-yE[d]+1),
				(n[d+1]>0)*pbeta(psi,yE[d+1]+1,n[d+1]-yE[d+1]+1))-c(0,0.001,0.002)       	
				d_opt=d-2+which.min(posH)		    	       	
		    	       } else {
			posH<-c(pbeta(psi,yE[d-1]+1,n[d-1]-yE[d-1]+1),
			pbeta(psi,yE[d]+1,n[d]-yE[d]+1))-c(0,0.001)
			d_opt=d-2+which.min(posH)		    	       	
		    	       	
		    	       }

			}
            }	
            if (elimiE[d_opt]==1) {earlystop=1; break} 
			d<-d_opt


		}

		YT[trial,]=yT;
		YE[trial,]=yE;
		N[trial,]=n;
        ntox=ntox+sum(yT)/ntrial
		neff=neff+sum(yE)/ntrial
		if (earlystop==0){
		pT<-(yT+0.05)/(n+0.1)
		pE<-(yE+0.05)/(n+0.1)
		pT<-pava(pT,n+0.1)+0.001*seq(1,ndose)

        pE<-peestimate(yE,n)
		pT[elimi==1]<-20
		pT[n==0]<-20
		pE[n==0]<-0
		pE[elimi==1]<-0
		u<-pE-0.33*pT-1.09*(pT>targetT) 	
		d_opt<-which.max(u)	
		dselect[trial]=d_opt
		sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
		#if (d_opt==2) {print(yT);print(yE);print(n)}
		} else {dselect[trial]<-99}	
		pts<-pts+n/ntrial
		dlt<-dlt+yT/ntrial
		eff<-eff+yE/ntrial	
	}	
	pts<-round(pts,1)
	dlt<-round(dlt,1)
	
	# results:
	# pT.true: true toxicity rate
	# pE.true: true efficacy rate
	# sel: selection percentage of each dose
	# pts: number of patients treated at each dose
	# dlt: number of DLTs observed at each dose
	# eff: number of efficacy outcomes observed at each dose
	# ntox: total number of DLTs
	# neff: total number of Effs
	
	results=list(pT.true=pT.true,pE.true=pE.true,sel=sel,pts=pts,dlt=dlt,eff=eff,ntox=ntox,neff=neff)
	return(results)	
}



target = 0.2 # target toxicity rate
cohortsize = 1 # cohort size
ncohort = 12 # number of cohorts
psi1<-0.35 # lowest acceptable efficacy rate
psi2<-0.65 # very promising efficacy rate that leads to dose retainment

pT.true1 <- c(0.1, 0.2, 0.3)
pT.true2 <- c(0.05, 0.22, 0.38)

pT.true3 <- c(0.2, 0.3, 0.4)
pT.true4 <- c(0.18, 0.3, 0.45)

pT.true5 <- c(0.07, 0.13, 0.21)
pT.true6 <- c(0.04, 0.1, 0.2)
#pT.true<-c(0.01,0.05,0.10,0.15,0.30)
pE.true1<-c(0.4, 0.5, 0.6)
pE.true2<-c(0.4, 0.5, 0.4)
pE.true3<-c(0.6, 0.5, 0.4)

res <- get.oc(target, pE.true3, 
              pT.true5, psi1,psi2, ncohort, cohortsize, startdose=2, cutoff.eli=0.95, ntrial=1000)
res$sel %>% round(1)
res$pts %>% round(1)
sum(res$pts) %>% round(1)
res$ntox %>% round(1)
res$neff %>% round(1)

