rm(list=ls())
library(MASS)
library(mvtnorm)
library(parallel)
library(foreach)
library(doParallel)
library(abind)
library(ggplot2)
library(mnonr)
library(BinOrdNonNor)

####1.distance function####
#calculate entropy-based distance of two cat points x and y
distance_f=function(x,y,entropy,maxentr,d_con,d_ord,d_nom){
    #x,y are vectors(cat~ord~nom)
    x_con=x[1:d_con]
    x_cat=x[(d_con+1):(d_con+d_ord+d_nom)]
    y_con=y[1:d_con]
    y_cat=y[(d_con+1):(d_con+d_ord+d_nom)]
    #continous part difference
    con_diff=(x_con-y_con)^2
    #cat part difference(ord~nom)
    cat_diff=rep(NA,(d_ord+d_nom))
    #ordinal part 
    for (i in 1:d_ord) {
        cat_diff[i]=(sum(entropy[x_cat[i]:y_cat[i],i])/maxentr[i])^2
    }
    #nominal part 
    for (i in (d_ord+1):(d_ord+d_nom)) {
        cat_diff[i]=((entropy[x_cat[i],i]+entropy[y_cat[i],i])/maxentr[i])^2
    }
    #mixed
    ED=sqrt(sum(con_diff)+sum(cat_diff))
    return(ED)
}

####2.find the UCL of different (m,n) setting ####
#2.1 given UCL£¬calculate RL
#2.1.1plotting statistic L_wa=WRS^2+AB^2
#mean&cov of WRS/AB
mm=m-1#size of sub-refdata(exclude one basepoint)
N=mm+n
#WRS
miu1=n*(N+1)/2;sig1=mm*n*(N+1)/12
#AB
if(N %% 2==0){
    miu2=n*N/4
}else{
    miu2=n*(N^2-1)/(4*N)
}
if(N %% 2==0){
    sig2=mm*n*(N^2-4)/(48*(N-1))
}else{
    sig2=mm*n*(N+1)*(N^2+3)/(48*N^2)
}
rl_f=function(m,n,miu1,miu2,sig1,sig2,plist,me,var,skew,kurt,cor,d_con,d_ord,d_nom,UCL){
    #m IC reference sample
    d_cat=d_ord+d_nom
    output=capture.output({
        refdata=genBinOrdNN(m,plist,me,var,skew,kurt,0,d_cat,d_con,cor)
    })
    refdata_cat=refdata[,1:d_cat]
    refdata_con=refdata[,(d_cat+1):(d_cat+d_con)]
    #con part:transform 
    samplecov=cov(refdata_con)
    trans=diag(1/sqrt(diag(samplecov)))
    refdata_con_trans=refdata_con %*% trans
    #cat part:entropy (ord~nom)
    prob_cat=apply(refdata_cat,2,function(x){prop.table(table(x))})
    entropy=-prob_cat*log(prob_cat)
    k_vec=sapply(plist, length)+1
    maxentr=log(k_vec)
    #basepoint
    basepointid=sample(1:m,1)
    basepoint_cat=refdata_cat[basepointid,]
    basepoint_con=refdata_con_trans[basepointid,]
    basepoint=c(basepoint_con,basepoint_cat)
    refsubdata_cat=refdata_cat[-basepointid,]
    refsubdata_con_trans=refdata_con_trans[-basepointid,]
    refsubdata_trans=cbind(refsubdata_con_trans,refsubdata_cat)
    #distance of sub-reference data
    refdist=apply(refsubdata_trans,1,function(x) {distance_f(x,basepoint,entropy,maxentr,d_con,d_ord,d_nom)})
    
    
    #RL for this reference sample
    plotsta=0
    RL=0
    while (plotsta<UCL && RL<=700) {
        RL=RL+1
        #testdata
        output=capture.output({
            testdata=genBinOrdNN(n,plist,me,var,skew,kurt,0,d_cat,d_con,cor)
        })
        testdata_cat=testdata[,1:d_cat]
        testdata_con=testdata[,(d_cat+1):(d_cat+d_con)]
        testdata_con_trans=testdata_con %*% trans
        
        testdata_trans=cbind(testdata_con_trans,testdata_cat)
        #distance of test data
        testdist=apply(testdata_trans,1,function(x){distance_f(x,basepoint,entropy,maxentr,d_con,d_ord,d_nom)})
        #rank of test data
        pooldist=c(testdist,refdist)
        testrank=rank(pooldist)[1:n]
        
        #plotting statistic
        #WRS
        T1=sum(testrank)
        #AB
        T2=sum(abs(testrank-((m-1+n)+1)/2))
        #lepage-style statistics
        plotsta=(T1-miu1)^2/sig1+(T2-miu2)^2/sig2
    }
    return(RL)
}

#2.1.2plotting statistic L_wm=WRS^2+M^
#mean&cov of WRS/Mood
mm=m-1#size of sub-refdata(exclude one basepoint)
N=mm+n
#WRS
miu1=n*(N+1)/2;sig1=mm*n*(N+1)/12
#mood
miu2=n*(N^2-1)/12
sig2=m*n*(N+1)*(N^2-4)/180
rl_f=function(m,n,miu1,miu2,sig1,sig2,plist,me,var,skew,kurt,cor,d_con,d_ord,d_nom,UCL){
    #m IC reference sample
    d_cat=d_ord+d_nom
    output=capture.output({
        refdata=genBinOrdNN(m,plist,me,var,skew,kurt,0,d_cat,d_con,cor)
    })
    refdata_cat=refdata[,1:d_cat]
    refdata_con=refdata[,(d_cat+1):(d_cat+d_con)]
    #con part:transform 
    samplecov=cov(refdata_con)
    trans=diag(1/sqrt(diag(samplecov)))
    refdata_con_trans=refdata_con %*% trans
    #cat part:entropy (ord~nom)
    prob_cat=apply(refdata_cat,2,function(x){prop.table(table(x))})
    entropy=-prob_cat*log(prob_cat)
    k_vec=sapply(plist, length)+1
    maxentr=log(k_vec)
    #basepoint
    basepointid=sample(1:m,1)
    basepoint_cat=refdata_cat[basepointid,]
    basepoint_con=refdata_con_trans[basepointid,]
    basepoint=c(basepoint_con,basepoint_cat)
    refsubdata_cat=refdata_cat[-basepointid,]
    refsubdata_con_trans=refdata_con_trans[-basepointid,]
    refsubdata_trans=cbind(refsubdata_con_trans,refsubdata_cat)
    #distance of sub-reference data
    refdist=apply(refsubdata_trans,1,function(x) {distance_f(x,basepoint,entropy,maxentr,d_con,d_ord,d_nom)})
    
    #RL for this reference sample
    plotsta=0
    RL=0
    while (plotsta<UCL && RL<=700) {
        RL=RL+1
        #testdata
        output=capture.output({
            testdata=genBinOrdNN(n,plist,me,var,skew,kurt,0,d_cat,d_con,cor)
        })
        testdata_cat=testdata[,1:d_cat]
        testdata_con=testdata[,(d_cat+1):(d_cat+d_con)]
        testdata_con_trans=testdata_con %*% trans
        
        testdata_trans=cbind(testdata_con_trans,testdata_cat)
        #distance of test data
        testdist=apply(testdata_trans,1,function(x){distance_f(x,basepoint,entropy,maxentr,d_con,d_ord,d_nom)})
        #rank of test data
        pooldist=c(testdist,refdist)
        testrank=rank(pooldist)[1:n]
        
        #plotting statistic
        mm=m-1#size of sub-refdata(exclude one basepoint)
        N=mm+n
        #WRS
        T1=sum(testrank)
        #mood
        T2=sum((testrank-(N+1)/2)^2)
        
        #lepage-style statistics
        plotsta=(T1-miu1)^2/sig1+(T2-miu2)^2/sig2
    }
    return(RL)
}

#2.1.3plotting statistic cucconi
mm=m-1#size of sub-refdata(exclude one basepoint)
N=mm+n
miu1=n*(N+1)*(2*N+1)
miu2=N+1
sig1=sqrt(m*n*(N+1)*(2*N+1)*(8*N+11)/5)
sig2=2*(N^2-4)/((2*N+1)*(8*N+11))-1
rl_f=function(m,n,miu1,miu2,sig1,sig2,plist,me,var,skew,kurt,cor,d_con,d_ord,d_nom,UCL){
    #m IC reference sample
    d_cat=d_ord+d_nom
    output=capture.output({
        refdata=genBinOrdNN(m,plist,me,var,skew,kurt,0,d_cat,d_con,cor)
    })
    refdata_cat=refdata[,1:d_cat]
    refdata_con=refdata[,(d_cat+1):(d_cat+d_con)]
    #con part:transform 
    samplecov=cov(refdata_con)
    trans=diag(1/sqrt(diag(samplecov)))
    refdata_con_trans=refdata_con %*% trans
    #cat part:entropy (ord~nom)
    prob_cat=apply(refdata_cat,2,function(x){prop.table(table(x))})
    entropy=-prob_cat*log(prob_cat)
    k_vec=sapply(plist, length)+1
    maxentr=log(k_vec)
    #basepoint
    basepointid=sample(1:m,1)
    basepoint_cat=refdata_cat[basepointid,]
    basepoint_con=refdata_con_trans[basepointid,]
    basepoint=c(basepoint_con,basepoint_cat)
    refsubdata_cat=refdata_cat[-basepointid,]
    refsubdata_con_trans=refdata_con_trans[-basepointid,]
    refsubdata_trans=cbind(refsubdata_con_trans,refsubdata_cat)
    #distance of sub-reference data
    refdist=apply(refsubdata_trans,1,function(x) {distance_f(x,basepoint,entropy,maxentr,d_con,d_ord,d_nom)})
    
    #RL for this reference sample
    plotsta=0
    RL=0
    while (plotsta<UCL && RL<=700) {
        RL=RL+1
        #testdata
        output=capture.output({
            testdata=genBinOrdNN(n,plist,me,var,skew,kurt,0,d_cat,d_con,cor)
        })
        testdata_cat=testdata[,1:d_cat]
        testdata_con=testdata[,(d_cat+1):(d_cat+d_con)]
        testdata_con_trans=testdata_con %*% trans
        
        testdata_trans=cbind(testdata_con_trans,testdata_cat)
        #distance of test data
        testdist=apply(testdata_trans,1,function(x){distance_f(x,basepoint,entropy,maxentr,d_con,d_ord,d_nom)})
        #rank of test data
        pooldist=c(testdist,refdist)
        testrank=rank(pooldist)[1:n]
        
        #plotting statistic
        #cucconi
        plotsta=(((6*sum((testrank)^2)-miu1)/sig1)^2+((6*sum((miu2-testrank)^2)-miu1)/sig1)^2-2*sig2*((6*sum((testrank)^2)-miu1)/sig1)*((6*sum((miu2-testrank)^2)-miu1)/sig1))/(2*(1-sig2^2))
    }
    return(RL)
}

#test 
t0=Sys.time();t0
rl0=foreach(1:sim,.combine = "c")%dopar%rl_f(m,n,miu1,miu2,sig1,sig2,plist,me,var,skew,kurt,cor,d_con,d_ord,d_nom,11.02)
t1=Sys.time();t=t1-t0;t
quantile(rl0, probs = c(0.25,0.5,0.75))
#

#2.4 search the UCL
t0=Sys.time();t0
uper=15;lower=7
mrl=mrl0+10
while (abs(mrl-mrl0)>1) {
    UCL=(uper+lower)/2
    print(UCL)
    rl=foreach(1:sim,.combine = "c")%dopar%rl_f(m,n,miu1,miu2,sig1,sig2,plist,me,var,skew,kurt,cor,d_con,d_ord,d_nom,UCL)
    mrl=median(rl)
    print(mrl)
    if(abs(mrl-mrl0)<1 | uper-lower<0.01){
        break
    }else{
        if(mrl-mrl0>0){
            uper=UCL
        }else{
            lower=UCL
        }
    }
    
}
t1=Sys.time();t=t1-t0;t
UCL

####3.MRL####
t0=Sys.time();t0
rl=foreach(1:sim,.combine = "c")%dopar%rl1_f(m,n,miu1,miu2,sig1,sig2,plist,me,var,plist,me,var,skew,kurt,cor,d_con,d_ord,d_nom,UCL)
t1=Sys.time();t=t1-t0;t
quantile(rl, probs = c(0.25,0.5,0.75))