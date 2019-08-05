##
##  Simulate population for non-nested analysis
##

library(sampling)
library(MASS)
library(survey)
library(svylme)

overlapsim<-function(N1,N2,n1,n2,overlap_f,
                     truebeta1=1, truebeta2=3,truesigma2=1,
                     truetau2_11=0.8,truetau_12=0.6,truetau2_22=1){
#setting: notation
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)
population$PSU<-population$long

overlap=ceiling(N2*overlap_f)

model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}

population<-model_cluster(population, overlap)
T=length(unique(population$cluster))

#Model: parameter from random slope model  model

PairCov<-matrix(c(truetau2_11, truetau_12, truetau_12, truetau2_22), nrow=2, byrow=T)

truevalue<-c(truebeta1,truebeta2, truesigma2, truetau2_11, truetau_12, truetau2_22)
names(truevalue)<-c("beta1", "beta2", "sigma2", "tau2_11", "tau_12", "tau2_22")

##Population data
re=mvrnorm(n=T, mu = c(0,0), Sigma = PairCov) #generate vector of random effect (a, b)
population$a<-re[,1][population$cluster]
population$b<-re[,2][population$cluster]

population$x<-rnorm(N1*N2)+rnorm(T)[population$cluster]
population$y<-with(population, truebeta1+a+truebeta2*x+b*x+rnorm(N1*N2,s=sqrt(truesigma2)))
population$r=with(population, x*(y-truebeta1-truebeta2*x))
population$ID_unit=with(population, 1:(N1*N2))

#uninformative two-stage sampling design (first-stage: SRSWOR, Second-stage:SRSWOR)
#n1=ceiling(N1/10) ##number of sampling cluster in the first stage (sample level)
#n2=ceiling(N2/10) ##umber of elements in each sampling cluster ( sample level)

# Using sampling package for two-stage sampling (First-stage: SRSWOR, Second-stage: SRSWOR ) 
##uninformative two-stage  sampling design (First-stage: SRSWOR, Second-stage: SRSWOR) and extracts the observed data
##first-stage
FirststageSRSWOR=srswor(n1, N1)
FirststageSRSWORSample=subset(population, population$PSU%in% which(FirststageSRSWOR==1))

#second-stage
SecondstageSRSWOR=unlist(lapply(rep(n2,n1), function(v) return(srswor(v, N2))))
TwostageSRSWORSample<-FirststageSRSWORSample[c(which(SecondstageSRSWOR==1)),] 
    TwostageSRSWORSample$p1<-n1/N1
        TwostageSRSWORSample$p2<-n2/N2

#informative two-stage sampling design (first-stage: SRSWOR, Second-stage:SRSWOR)
##number of elements in each sampling cluster
param=c(1.5, 0.45)

n2informative= function(r, sc, param, N2){
   a=rep(NA, length=length(unique(population$sc)))
   b=rep(NA, length=length(unique(population$sc)))
   for (i in unique(sc)){
      a[i]=mean(r[sc==i])
      b[i]=2*ceiling((param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))*N2/2)
   }
   b
}

##informative two-stage  sampling design (SRSWOR)[second-stage is informative] and extracts the observed data
###second-stage
    n2pop=n2informative(population$r,population$PSU, param ,N2)
    n2is=n2pop*FirststageSRSWOR
    SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
    TwostageSRSWORSampleis=FirststageSRSWORSample[c(which(SecondstageSRSWORis==1)), ]
    TwostageSRSWORSampleis$p1<-n1/N1
    TwostageSRSWORSampleis$p2<-rep(n2is[n2is!=0],n2is[n2is!=0],)/N2
    list(ignorable=TwostageSRSWORSample, informative=TwostageSRSWORSampleis)
}



rr<-replicate(1000,{
    cat(".")

a<-overlapsim(400,500,20,25,overlap_f=0.8)
design_ign<-svydesign(id=~PSU+ID_unit,prob=~p1+p2,data=a$ignorable)
m_ig<-svy2lme(y~(1+x|cluster)+x,design=design_ign)


    
design_inf<-svydesign(id=~PSU+ID_unit,prob=~p1+p2,data=a$informative)
m_inf<-svy2lme(y~(1+x|cluster)+x,design=design_inf)

c(ig=c(coef(m_ig), m_ig$opt$par, m_ig$s2, SE(m_ig)), inf=c(coef(m_inf), m_inf$opt$par, m_inf$s2, SE(m_inf)))
})

