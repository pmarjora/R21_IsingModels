library(tidyverse)
library(binaryLogic)  # for using binary numbers
library(Rmpfr)
library(mcmc)  # R's MCMC library
library(coda)

set.seed(456)
NumberOfSamples <- 100
NumberOfSites <- 20   # Can go up to about 12 at the moment
LengthOfMCMCRuns <- 50000
MCMC_Flip_Prob <- 0.25


MethData<-as_tibble(mat.or.vec(nr=NumberOfSamples,nc=NumberOfSites))
Locations=as_tibble(seq(0,1,1/(NumberOfSites-1)))


for (i in 1:NumberOfSamples){
  for (j in 1:NumberOfSites){
    MethData[i,j]<-sample(c(0,1),size=1)
  }
}


EnergyCalculator <- function(MethSequence,a,c){
  Exponent <- 0
  for (i in 1:length(MethSequence)){
    Exponent <- Exponent + a*(2*as.numeric(MethSequence[i])-1)
    if (i>1){
      Exponent <- Exponent + c*(2*as.numeric(MethSequence[i])-1)*(2*as.numeric(MethSequence[i-1])-1)
    }
  }
  Energy <- exp(Exponent)
  return (Energy)
}



CalcNormalizingConstant <- function(NSites,a,c){
  z <- 0
  for (i in 1:(2^NSites)){
    # work our way through each possible methylation pattern by converting their binary representations into a vector
    BinaryRep<-as.binary(i)
    z <- z + EnergyCalculator(BinaryRep,a,c)
    #cat("\n ",z,"  ",EnergyCalculator(BinaryRep,a,c))
  }
  return (z)
}

# function to produce samples from the IsingModel for given parameter values. It uses MCMC
MCMCsampler <- function(NSites,a,c,LengthOfRun,SamplingFreq,BurnIn,PFlip)
{
  Seq <- sample(c(0,1),NSites,replace=TRUE)
  EntireRun <- mat.or.vec(nr=LengthOfRun,nc=NSites)
  EntireRun[1,] <- Seq
  for (i in 2:LengthOfRun)
  {
    CurrentSeq <- Seq
    # Propose new state. To do this we flip each site with prob PFlip
    P<-runif(NSites)
    for (j in 1:NSites){
      if (P[j]<PFlip){Seq[j] <- 1-Seq[j]}
    }
    
    # Calculate the Hastings Ratio
    HR <- EnergyCalculator(Seq,a,c)/EnergyCalculator(CurrentSeq,a,c)
    #cat("\nHR=",HR,"    Prop: ",Seq)
    Q<-runif(1)
    # Reject the new state with prob 1-HR
    if (Q>HR){Seq<-CurrentSeq}#else{cat("  Accepted")}
    #cat("\n Seq=",Seq)
    EntireRun[i,] <- Seq
    
    
    # output the current state if appropriate
    if (i>BurnIn){
      if ((i-BurnIn-1)%%SamplingFreq==0){
        # cat("\n It ",i,"   Seq=",Seq)
      }
    }
    
  }
  EntireRun<-mcmc(EntireRun)
  return(EntireRun)
}

DoAllMCMCRuns<-function(A,C,NSites,NLength){
  
  data1 <- MCMCsampler(NumberOfSites,IsingParam_A,IsingParam_C,LengthOfMCMCRuns,100,20,MCMC_Flip_Prob)
  data2 <- MCMCsampler(NumberOfSites,IsingParam_A,IsingParam_C,LengthOfMCMCRuns,100,20,MCMC_Flip_Prob)
  data3 <- MCMCsampler(NumberOfSites,IsingParam_A,IsingParam_C,LengthOfMCMCRuns,100,20,MCMC_Flip_Prob)
  data4 <- MCMCsampler(NumberOfSites,IsingParam_A,IsingParam_C,LengthOfMCMCRuns,100,20,MCMC_Flip_Prob)
  
  MCMC1 <- mcmc(data1)
  MCMC2 <- mcmc(data2)
  MCMC3 <- mcmc(data3)
  MCMC4 <- mcmc(data4)
  
  Combined <- mcmc.list(list(data1,data2,data3,data4))
  return(Combined)
  
  
}

#NormConst<-CalcNormalizingConstant(NumberOfSites,a=2,c=1)

# Simulate some Ising model data

# Scenario 1  a=0, c=0
# Scenario 2 a=0.5  c=0
# scenario 3  a=0  c=0.5
# scenario 4  a=0.5   c=0.5
# Scenario 5 a=0.5  c=-0.5
# Scenario 6  a=1     c=0

IsingParam_A <- 0
IsingParam_C <- 1
Scen<- "Scenario7"
AllRuns <- DoAllMCMCRuns(IsingParam_A,IsingParam_C,NumberOfSites,LengthOfMCMCRuns)
SamplingFrequency <- 400
BurnIn <- 500

# gelman tests
#gelman.plot(AllRuns) # for plots
gelman.diag(AllRuns) # for diagnostic values

ops<- par()
par(mfrow=c(3,3))
for (j in 1:9){acf(AllRuns[[1]][,j],lag.max=400)}
ops<- par()

CreateDataSet <- function(MCMCoutput,NSamples,SampFreq,BurnInTime){
  Sample<-mat.or.vec(nr=NSamples,nc <- length(MCMCoutput[1,]))
  LengthOfRun <- length(MCMCoutput[,1])
  i <- BurnInTime
  nsamp <- 0
  while ((i<LengthOfRun)&&(nsamp<NSamples)){
    i <- i+1
    if (i %% SampFreq == 0){
      nsamp <- nsamp +1
      Sample[nsamp,] <- MCMCoutput[i,]
    }
  }
  if (nsamp<NSamples){
    cat("\nNote: There were not enough lines of MCMC output to complete the sampling.")
  }
  return (Sample)
}

CreatePooledData <- function(Samples){
   return (colMeans(Samples))
}

# create the data
DataSet1 <- CreateDataSet(AllRuns[[1]],NumberOfSamples,SamplingFrequency,BurnIn)
DataSet2 <- CreateDataSet(AllRuns[[2]],NumberOfSamples,SamplingFrequency,BurnIn)
DataSet3 <- CreateDataSet(AllRuns[[3]],NumberOfSamples,SamplingFrequency,BurnIn)
DataSet4 <- CreateDataSet(AllRuns[[4]],NumberOfSamples,SamplingFrequency,BurnIn)

AllData <- list("S1"=DataSet1,"S2"=DataSet2,"S2"=DataSet3)

Pool1<-CreatePooledData(DataSet1)
Pool2<-CreatePooledData(DataSet2)
Pool3<-CreatePooledData(DataSet3)
Pool4<-CreatePooledData(DataSet4)


write.table(Pool1,file=paste(Scen,"Pool1.txt"),quote=FALSE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ")
write.table(Pool2,file=paste(Scen,"Pool2.txt"),quote=FALSE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ")
write.table(Pool3,file=paste(Scen,"Pool3.txt"),quote=FALSE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ")
write.table(Pool4,file=paste(Scen,"Pool4.txt"),quote=FALSE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ")
write.table(DataSet1,file=paste(Scen,"IndividualLevelData1.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=" ")  
write.table(DataSet2,file=paste(Scen,"IndividualLevelData2.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=" ")  
write.table(DataSet3,file=paste(Scen,"IndividualLevelData3.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=" ")                                    
write.table(DataSet4,file=paste(Scen,"IndividualLevelData4.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=" ")                                    


