#'ICARBHSampler
#'  This function fits a piecewise hazard using a hierarchical model with a ICAR dependence applied to the hazard heights
#' @param Y - This is a n-vector containing patient survival times
#' @param I - This is a n-vector containing patient censoring indicators (0 for censored patient)
#' @param B - Number of iterations to run the sample
#' @param hyper - Vector of hyperparameters. In order, this contains a1, b1 which are the inverse gamma hyperparameters on sigma^2. Phi which is the hyperparameter on the mean number of split points. Jmax which is the maximum allowed number of split points.  cl1 which is a tuning parameter greater than 0. J1 is the starting number of split points for the MCMC. Finally, clam1 which is between 0 and 1 and characterizes the spatial dependency of the baseline hazard heights.
#' @return Returns a list containing the posterior samples of the split points, split point locations, log hazard rates and hierarchical samples
#'@import mvtnorm
#'@import graphics
#'@import stats
#' @examples
#' ####This generates random survival data
#'  Y=rexp(100,1/20)
#'  I=rbinom(100,1,.5)
#'  ###Sets hyperparameters
#'  a1=.7
#'  b1=.7
#'  phi=3
#'  Jmax=20
#'  cl1=.25
#'  clam1=.5
#'  J1=3
#'  ###Combines the hyperparameters in to a vector
#'  hyper=c(a1,b1,phi,Jmax,cl1,J1,clam1)
#'  ###Set Number of iterations
#'  B=100
#'  ###Run the Sampler
#'  X=ICARBHSampler(Y,I,B,hyper)
#'  X
#'
#' @references
#' Lee, K. H., Haneuse, S., Schrag, D. and Dominici, F. (2015), Bayesian semiparametric analysis of semicompeting risks data: investigating hospital readmission after a pancreatic cancer diagnosis. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 253-273.
#'
#' @export
ICARBHSampler=function(Y,I,B,hyper){

  a1=hyper[1]
  b1=hyper[2]
  phi=hyper[3]
  Jmax=hyper[4]
  cl1=hyper[5]
  J1=hyper[6]
  clam1=hyper[7]

  G1=J1+1
  if(J1>Jmax){
    cat("Number of starting split points must be less than the maximum allowed (Jmax)")
  }
  else{



###InPut Hyperparameters from vector


Mulam1=rep(0,B)
Siglam1=rep(1,B)



m1=max(Y)
#########################S Matrices!!!
#Reset up lam and S1 matrices
s1=matrix(rep(NA,B*(Jmax+2)),nrow=B)
s1[1,1:(J1+2)]=sort(seq(0,m1,length.out = J1+2))


lam1=matrix(rep(NA,B*(Jmax+1)),nrow=B)

lam1[1,1:(J1+1)]=rep(0,J1+1)


split=rep(1,B)





##Likelihood
LK = function(Y,I,J1,s,lam){

  LOGBH=0
  G1=J1+1
  for(k in 1:G1){
    Del=pmax(0,pmin(Y,s[k+1])-s[k])



    LOGBH=LOGBH-sum(Del*exp(lam[k]))

    zu=Y<=s[k+1]
    zl=Y>s[k]
    LOGBH=LOGBH+sum(zu*zl)*lam[k]
  }
  return(LOGBH)

}
##baselinehazardsamplingfunction


base_samp = function(Y,I,J1,lam1,s1,Sigma){

  #lambdaB
  #######
  J1=J1
  m1=max(Y[I==1])
  G1=J1+1

  ##Matrix Setup
  if(J1>0){SigLam1=Sigma}else{SigLam1=as.matrix(1)}




  ##




  for(m in 1:(J1+1)){

    lam=lam1
    lam=lam[is.na(lam)==FALSE]
    lambda=lam

    lam[m]=lambda[m]+runif(1,-cl1,cl1)


    if(J1==0){
      do=log(dnorm(lambda[m],Mu,sqrt(Sig)))
      dn=log(dnorm(lam[m],Mu,sqrt(Sig)))
    }else{


      do=dmvnorm(lambda,rep(Mu,J1+1),Sig*SigLam1)
      dn=dmvnorm(lam,rep(Mu,J1+1),Sig*SigLam1)
    }

    Likeo=LK(Y,I,J1,s1,lambda)

    Liken=LK(Y,I,J1,s1,lam)



    U=log(runif(1,0,1))
    alphalam=Liken-Likeo+dn-do

    if(is.na(alphalam)==FALSE){

      if(U<alphalam){
        lam1[m]=lam[m]
      }
    }


  }




  ###BD



  ###Random Perturbation###
  U1=runif(1,0,1)
  #####

  s=s1
  s=s[!is.na(s)]

  if(length(s)<(Jmax+1)){
    Birth=runif(1,0,m1)

    s2=sort(c(s,Birth))

    for(k in 2:(J1+2)){
      if(Birth>s1[k-1] && Birth<s1[k]){
        Ind=k-1
      }
    }

    lam=rep(0,J1+2)

    if(Ind==1 | Ind==J1+1){
      if(Ind==1){
        lam[Ind]=lam1[Ind] - ((s1[Ind+1]-Birth)/(s1[Ind+1]-s1[Ind]))*log((1-U1)/U1)
        lam[Ind+1]=lam1[Ind] + ((Birth-s1[Ind])/(s1[Ind+1]-s1[Ind]))*log((1-U1)/U1)
        lam[(Ind+2):length(lam)]=lam1[(Ind+1):(J1+1)]
      }else{
        lam[Ind]=lam1[Ind] - ((s1[Ind+1]-Birth)/(s1[Ind+1]-s1[Ind]))*log((1-U1)/U1)
        lam[Ind+1]=lam1[Ind] + ((Birth-s1[Ind])/(s1[Ind+1]-s1[Ind]))*log((1-U1)/U1)
        lam[1:(Ind-1)]=lam1[1:(Ind-1)]
      }
    }else{
      lam[Ind]=lam1[Ind] - ((s1[Ind+1]-Birth)/(s1[Ind+1]-s1[Ind]))*log((1-U1)/U1)
      lam[Ind+1]=lam1[Ind] + ((Birth-s1[Ind])/(s1[Ind+1]-s1[Ind]))*log((1-U1)/U1)
      lam[1:(Ind-1)]=lam1[1:(Ind-1)]
      lam[(Ind+2):length(lam)]=lam1[(Ind+1):(J1+1)]
    }

    lam=lam[!is.na(lam)]

    lambda=lam1
    lambda=lambda[!is.na(lambda)]


    ##


    Lo=LK(Y,I,J1,s1,lam1)



    if(J1>0){
      do=log(dpois(J1,phi))+log(dmvnorm(lambda,rep(Mu,length(lambda)),SigLam1*Sig))
    }else{
      do=log(dpois(J1,phi))+log(dnorm(lambda,Mu,sqrt(Sig)))
    }

    prior=((2*J1+3)*(2*J1+2)*(Birth-s1[Ind])*(s1[Ind+1]-Birth))/((m1^2)*(s1[Ind+1]-s1[Ind]))

    G1=G1+1
    J1=J1+1


    Ln=LK(Y,I,J1,s2,lam)



    ##Make SigLamB



    W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
    Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


    length1=diff(s2[!(is.na(s2))])


    if(J1<2){
      if(J1==1){
        W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
        W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
        Q1[1,1]=2/(2*length1[1]+length1[2])
        Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
        SigLam2=solve(diag(J1+1)-W1)%*%Q1





      }else{

        SigLam2=as.matrix(1)
      }
    }else{


      for(j in 2:J1){
        W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
        W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
        Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
      }


      Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
      Q1[1,1]=2/(2*length1[1]+length1[2])
      W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
      W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


      SigLam2=solve(diag(J1+1)-W1)%*%Q1

    }



    dn=log(dpois(J1,phi))+log(dmvnorm(lam,rep(Mu,length(lam)),Sig*SigLam2))




    alpha=Ln-Lo+dn-do-log(U1*(1-U1)) + log(prior)

    if(is.na(alpha)){
      J1=J1-1
      G1=G1-1
    }else{

      U=log(runif(1,0,1))

      if(U<alpha){
        lam1[1:(J1+1)]=lam
        s1[1:(J1+2)]=s2
        SigLam=SigLam2
      }else{
        J1=J1-1
        G1=G1-1
        SigLam=SigLam1
      }

    }


  }


  #########################################################
  ###################Death sBmpler#########################
  ##########################################################

  U1=runif(1,0,1)

  if(J1>0){


    if(J1==1){
      Ind=2
    }else{

      Ind=sample(2:(J1+1),1)
    }


    s=s1
    s2=s[-Ind]

    lam=lam1
    lambda=lam[!is.na(lam)]

    lam=lam[!is.na(lam)]
    lam=lam[-Ind]

    lam[Ind-1]=((s1[Ind]-s1[Ind-1])*lam1[Ind-1]+(s1[Ind+1]-s1[Ind])*lam1[Ind])/(s1[Ind+1]-s1[Ind-1])



    #############################################
    ####Sets up SigLamB matrix for old density###
    #############################################


    W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
    Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)


    length1=diff(s1[!(is.na(s1))])



    if(J1<2){
      if(J1==1){
        W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
        W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
        Q1[1,1]=2/(2*length1[1]+length1[2])
        Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
        SigLam1=solve(diag(J1+1)-W1)%*%Q1

        do=log(dpois(J1,phi))+log(dmvnorm(lambda,rep(Mu,length(lambda)),SigLam1*Sig))




      }else{


        do=log(dpois(J1,phi))+log(dnorm(lambda,Mu,sqrt(Sig)))

      }
    }else{

      for(j in 2:J1){
        W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
        W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
        Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
      }


      Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
      Q1[1,1]=2/(2*length1[1]+length1[2])
      W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
      W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


      SigLam1=solve(diag(J1+1)-W1)%*%Q1

      do=log(dpois(J1,phi))+log(dmvnorm(lambda,rep(Mu,length(lambda)),SigLam1*Sig))


    }
    #############################################
    #############################################


    Lo=LK(Y,I,J1,s1,lam1)



    prior=((m1^2)*(s1[Ind+1]-s1[Ind-1]))/((2*J1+1)*(2*J1)*(s1[Ind]-s1[Ind-1])*(s1[Ind+1]-s1[Ind]))


    G1=G1-1
    J1=J1-1




    Ln=LK(Y,I,J1,s2,lam)


    ###Make siglam matrix

    if(J1>0){

      W1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)
      Q1=matrix(rep(0,(J1+1)*(J1+1)),nrow=J1+1)



      length1=diff(s2[!(is.na(s2))])






      if(J1<2){
        if(J1==1){
          W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
          W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])
          Q1[1,1]=2/(2*length1[1]+length1[2])
          Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
          SigLam2=solve(diag(J1+1)-W1)%*%Q1


          dn=log(dpois(J1,phi))+log(dmvnorm(lam,rep(Mu,length(lam)),SigLam2*Sig))



        }
      }else{


        for(j in 2:J1){
          W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
          W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
          Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
        }


        Q1[J1+1,J1+1]=2/(length1[J1]+2*length1[J1+1])
        Q1[1,1]=2/(2*length1[1]+length1[2])
        W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
        W1[J1+1,J1]=(clam1*(length1[J1+1]+length1[J1]))/(length1[J1]+2*length1[J1+1])


        SigLam2=solve(diag(J1+1)-W1)%*%Q1

        dn=log(dpois(J1,phi))+log(dmvnorm(lam,rep(Mu,length(lam)),SigLam2*Sig))


      }
      ####

    }else{

      dn=log(dpois(J1,phi))+log(dnorm(lam,Mu,sqrt(Sig)))
      SigLam2=as.matrix(1)

    }


    alpha=Ln-Lo+dn-do+log(prior)+log(U1*(1-U1))

    if(is.na(alpha)){
      J1=J1+1
      G1=G1+1
    }else{

      U=log(runif(1,0,1))


      if(U<alpha){
        s1=c(s2,NA)
        lam1[1:(J1+1)]=lam
        lam1[(J1+2):Jmax]=rep(NA,Jmax-J1-1)
        SigLam=SigLam2
      }else{
        J1=J1+1
        G1=G1+1
        SigLam=SigLam1
      }
    }

    ####End else
  }
  ##




  return(list(J1,s1,lam1,SigLam))


}



##This function returns a posterior MCMC sample, the whole thing, so size B and treatment vector Trt


JA=J1
sA=s1
lamA=lam1

###First Make Siglam Matrix ###
if(JA>0){
  W1=matrix(rep(0,(JA+1)*(JA+1)),nrow=JA+1)
  Q1=matrix(rep(0,(JA+1)*(JA+1)),nrow=JA+1)
  length1=rep(0,JA+1)
  for(j in 1:length(length1)){
    length1[j]=sA[1,j+1]-sA[1,j]
  }
  if(JA<2){
    if(JA==1){
      W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
      W1[JA+1,JA]=(clam1*(length1[JA+1]+length1[JA]))/(length1[JA]+2*length1[JA+1])
      Q1[1,1]=2/(2*length1[1]+length1[2])
      Q1[JA+1,JA+1]=2/(length1[JA]+2*length1[JA+1])
      SigLam1=solve(diag(JA+1)-W1)%*%Q1
    }
  }else{
    for(j in 2:JA){
      W1[j,j-1]=(clam1*(length1[j]+length1[j-1]))/(length1[j-1]+2*length1[j]+length1[j+1])
      W1[j,j+1]=(clam1*(length1[j]+length1[j+1]))/(length1[j-1]+2*length1[j]+length1[j+1])
      Q1[j,j]=2/(length1[j-1]+2*length1[j]+length1[j+1])
    }
    Q1[JA+1,JA+1]=2/(length1[JA]+2*length1[JA+1])
    Q1[1,1]=2/(2*length1[1]+length1[2])
    W1[1,2]=(clam1*(length1[1]+length1[2]))/(2*length1[1]+length1[2])
    W1[JA+1,JA]=(clam1*(length1[JA+1]+length1[JA]))/(length1[JA]+2*length1[JA+1])
    SigLam1=solve(diag(JA+1)-W1)%*%Q1
  }
}else{SigLam1=as.matrix(1)}
###Set SiglamA=Siglam1
SigLam1=SigLam1







for(b in 2:B){
  ##Outline:
  ##Hierarchical Lambda
  ## clambda
  ##Lambda
  ## Add
  ## Delete




  iter="Hier"

  ##Lambda1 Hierarchical Sampler
  ##Mulam
  L1=lam1[b-1,1:G1]



  ##Use 2 cores to do baseline samplers

  iter[2]="Mu"
  ##Lambda1 Hierarchical Sampler
  ##Mulam

  if(J1>0){

    Mulam1[b]=rnorm(1,(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%L1)/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1))),sqrt(Siglam1[b-1]/(t(as.matrix(rep(1,J1+1)))%*%solve(SigLam1)%*%as.matrix(rep(1,J1+1)))))


    Siglam1[b]=1/rgamma(1,a1+(J1+1)/2,b1+.5*(t(as.matrix(rep(Mulam1[b],J1+1))-L1)%*%solve(SigLam1)%*%(as.matrix(rep(Mulam1[b],J1+1))-L1)))


    ##Siglam

    iter[2]="Sigma"
  }else{



    Mulam1[b]=rnorm(1,lam1[b-1,1],sqrt(Siglam1[b-1]))


    Siglam1[b]=1/rgamma(1,a1+1/2,b1+.5*(Mulam1[b]-lam1[b-1,1])^2)



  }







  Mu=Mulam1[b]
  Sig=Siglam1[b]


  ##Try Sequentially First ##

  iter="baselineA"

  X=base_samp(Y,I,J1,lam1[b-1,],s1[b-1,],SigLam1)
  J1=X[[1]]
  G1=J1+1
  lam1[b,]=X[[3]]
  s1[b,]=X[[2]]
  SigLam1=X[[4]]


  split[b]=J1





  ###End For Loop
}


par(mfrow=c(3,1))

hist(split, main="Histogram of Number of SplitPoints")


plot(1:B,Mulam1,type="l",ylab="Mulam", main="Trace of Mulam")
plot(1:B,Siglam1,type="l",ylab="Siglam", main="Trace of Siglam")



z=list(split,s1,lam1,Mulam1,Siglam1)

return(z)




}


}








