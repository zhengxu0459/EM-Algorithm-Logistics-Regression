library(PoissonBinomial)
library(MASS)
library(mvtnorm)

Aggregate_Logit_Reg_Minus_Log_Likelihood = 
  function(Para, DataX, AggregatedDataY, Add_Intercept = TRUE, Mode = 2 ){
  Theta = Para
  if (Add_Intercept == TRUE){
    DataX = cbind(1,DataX)
  }
  DataX = as.matrix(DataX)
  nx = dim(as.matrix(DataX))[1]
  ny = length(AggregatedDataY)
  Aggregate = nx/ny
  Eta = DataX %*% Theta
  Prob = 1/(1+exp(-Eta))
  Prob_Matrix = matrix(Prob, ncol= Aggregate,  byrow = TRUE)
  if (Mode == 1){
    cat("Mode=",Mode,"\n") # For debug purpose.
    Loglike = rep(NA, ny)
    for (i in 1:ny){
      Loglike[i] = dpbinom( AggregatedDataY[i] , probs = Prob_Matrix[i,], log = TRUE)
    }
  }
  if (Mode==2){
      helper_dpbinom = function(vector, log=TRUE){
        return( dpbinom(vector[Aggregate+1],vector[1:Aggregate],log=TRUE ))
      }
      Temp = cbind(Prob_Matrix,AggregatedDataY)
      Loglike = apply(Temp,1,helper_dpbinom,log=TRUE)
  }
  return(-sum(Loglike))
}

Aggregate_Logit_Reg = function(DataX, AggregatedDataY, 
                               Add_Intercept = TRUE, Start=NULL, Mode = 2,
                               Optimization_Method = "Nelder-Mead"){
  if (is.null(Start)==TRUE){
    p = dim(as.matrix(DataX))[2] + (Add_Intercept + 0.0)
    Start = rep(0,p) + 0.01
  }
  value = Aggregate_Logit_Reg_Minus_Log_Likelihood(Start,DataX,AggregatedDataY,
                                                   Add_Intercept = Add_Intercept, Mode = Mode)
  if (is.infinite(value)){
    p = dim(as.matrix(DataX))[2] + (Add_Intercept + 0.0)
    Start = rep(0,p) + 0.01
  }
  w = optim(par = Start, Aggregate_Logit_Reg_Minus_Log_Likelihood,
            DataX=DataX, AggregatedDataY=AggregatedDataY, Add_Intercept = Add_Intercept, Mode = Mode,
            method = Optimization_Method)
  return(w)
}

Log_Conditional_Probability_For_EM = 
  function(Para, DataX, AggregatedDataY, Add_Intercept = TRUE, Mode = 2, Regulate = TRUE ){
  Theta = Para
  if (Add_Intercept == TRUE){
    DataX = cbind(1,DataX)
  }
  DataX = as.matrix(DataX)
  nx = dim(as.matrix(DataX))[1]
  ny = length(AggregatedDataY)
  Aggregate = nx/ny
  Eta = DataX %*% Theta
  Prob = 1/(1+exp(-Eta))
  Prob_Matrix = matrix(Prob, ncol= Aggregate,  byrow = TRUE)
  helper_dpbinom = function(vector, log=TRUE){
	return( dpbinom(vector[length(vector)],vector[1:(length(vector)-1)],log=TRUE ))
  }
  Temp = cbind(Prob_Matrix,AggregatedDataY)
  Loglike = apply(Temp,1,helper_dpbinom,log=TRUE)
  LeaveOneOut.Loglike = matrix(NA,nrow=ny,ncol=Aggregate)
  for (j in 1:Aggregate){
	TempUse = Temp[,-j]
	TempUse[,dim(TempUse)[2]] = TempUse[,dim(TempUse)[2]] - 1
	LeaveOneOut.Loglike[,j] = apply(TempUse,1,helper_dpbinom,log=TRUE)
  }
  temp = log(Prob_Matrix)+LeaveOneOut.Loglike-
	matrix(Loglike,nrow=length(Loglike),ncol=Aggregate)
  temp.vector = as.vector(t(temp))
  
  return(temp.vector) 
}

Aggregate_Logit_EM_Estimator = function(DataX, AggregatedDataY, 
                               Add_Intercept = TRUE, Start=NULL, Mode = 2,MaxIteration = 100, 
                               RelativeTolerateXNorm = 1e-10, RelativeTolerateY = 1e-10,
                               TolerateXNorm = 1e-10, TolerateY = 1e-10,
                               Verbose = 3, 
                               KeepHistory = TRUE
                               ){
  nx = dim(as.matrix(DataX))[1]
  ny = length(AggregatedDataY)
  Group.Size = nx/ny
  
  if (is.null(Start)==TRUE){
    data.aggregate.x = matrix( apply( matrix(DataX,Group.Size,length(AggregatedDataY)*dim(DataX)[2]), 
                                      2, sum),
                             length(AggregatedDataY), dim(DataX)[2] )
    data.aggregate.mean.x = data.aggregate.x / Group.Size
    glm.fit.1=glm( cbind(AggregatedDataY,Group.Size-AggregatedDataY)~data.aggregate.mean.x, family="binomial")
    Start = glm.fit.1$coefficients
    value = Aggregate_Logit_Reg_Minus_Log_Likelihood(Start,DataX,AggregatedDataY,
                                                   Add_Intercept = Add_Intercept, Mode = Mode)
  }
  if (Add_Intercept == FALSE){
    stop("Add_Intercept should be set as TRUE in the current version.")
  }
  theta.old = Start
  for (Iteration in 1:MaxIteration){
    if (Verbose>=3){
      print(Iteration)
    }
    log.weight = Log_Conditional_Probability_For_EM(Para = theta.old, DataX = DataX, 
                                                    AggregatedDataY = AggregatedDataY, 
                                                    Add_Intercept = TRUE, Mode = 2 )
    PseudoY = c(rep(1,nx),rep(0,nx))
    PseudoWeight = c(exp(log.weight),1-exp(log.weight))
    PseudoX = rbind(DataX,DataX)
    glm.fit.2=glm( PseudoY~PseudoX, family="binomial", weights = PseudoWeight)
    theta.new = glm.fit.2$coefficients
    absolute.distance.in.x = sqrt(sum((theta.new - theta.old)^2))
    relative.distance.in.x = 1 
    value.new = Aggregate_Logit_Reg_Minus_Log_Likelihood(theta.new, DataX, 
                                                         AggregatedDataY, Add_Intercept = TRUE, Mode = 2 ) 
    value.old = Aggregate_Logit_Reg_Minus_Log_Likelihood(theta.old, DataX, 
                                                         AggregatedDataY, Add_Intercept = TRUE, Mode = 2 ) 
    absolute.distance.in.y = abs(value.new - value.old) 
    relative.distance.in.y = 1 
    if ((absolute.distance.in.x<TolerateXNorm )|(relative.distance.in.x<RelativeTolerateXNorm )) {
     break
    }
    if ((absolute.distance.in.y<TolerateY )|(relative.distance.in.y<RelativeTolerateY )) {
     break
    }
    theta.old = theta.new 
  }
  Convergence = TRUE
  if (Iteration == MaxIteration){
    Convergence = FALSE
  }
  if (Verbose >=2){
    print (list(theta.old=theta.old, theta.new=theta.new, Iteration = Iteration, 
                Minus_Log_Likelihood_Value = value.new,
              MaxIteration = MaxIteration, Convergence = Convergence, history=NULL))
  }
  return(list(theta.old=theta.old, theta.new=theta.new, Iteration = Iteration, 
              Minus_Log_Likelihood_Value = value.new, 
              MaxIteration = MaxIteration, Convergence = Convergence, history=NULL))
}


