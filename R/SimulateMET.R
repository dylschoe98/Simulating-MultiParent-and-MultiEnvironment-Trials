# ==== Generate simulated phenotypic values in a MET ==== #
#' @param dt a simulated data.table with columns for 'Genotype','Env',and 'Rep' from the function `SimEnvironments`
#' @param gt a simulated data.table with columns for 'Genotype','BreedVal', and 'GenoVal' from `SimPhenotypes` 
#' @param Verror a numeric vector of which each value is betwen between 0 and 1 corresponding to the within plot error variance. If a vector is provided, a unique variance is provided to each environment by randomly sampling out each of the variances per environment.   
#' @param u a numeric value for the overall mean 
#' @param Ve an optional numeric value for the variance of the environment effect
#' @param Vblock an optional numeric value for a variance component due to the block term. Can be used as a term for replication in the case of an RCBD where the block is replication.
#' @param Vfield an optional numeric value for a variance component due to the field effect (i.e., meant to represent row and column effects). 
#' @return Data.table containing
#' \describe{
#' \item{SimulatedMET}{A data.table with the original'Genotype','Env',and 'Rep' assignments plus a column called 'PhenoVal' for the simulated phenotypes.
#' }

library(data.table)
library(dplyr)

SimMET<-function(dt,gt,Verror,u,Ve=0,Vblock=0,Vfield=0){
  
#An internal function to create a correlated random vector
CorVariable <- function(y, rho){ #https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
    x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
        y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}
  
#Correct column names
colnames(dt)[1:3]<-c('Genotype','Env','Rep')
colnames(gt)<-c('Genotype','Env','BreedVal','GenoVal')

#Correcting for values that could be NULL
ENV<-copy(Ve)
Block<-copy(Vblock)
Field<-copy(Vfield)

# Number of environments
nEnv<-length(unique(dt$Env))
#..Create names for Rep nested within the environment 
  dt$RE<-paste0(dt$Rep,'_',dt$Env)
    re<-unique(dt[Rep >1]$RE) #can be more than 1 rep
    
Verror<-copy(Verror)
  REPS<-unique(dt$Rep)
#...adding a replicate effect such that it is correlated      
    NoiseList<-list()
      env<-unique(dt$Env)
for (j in 1:length(env)){
  
#Randomly sample out the Mean for the environment if provided 
if(length(u) >1){  
  EnvMean<-sample(x = u,1)    
    dt[Env==env[j],EnvAverage:=EnvMean]
}else{
  dt[Env==env[j],EnvAverage:=u]
}  
#..Count the number of genotypes in each rep
    countsRep<-dt[Env==env[j],.N, by=.(Genotype)]
#...sample out a unique error term for each environment 
      if(length(Verror) >1){ErrorSelection<-sample(Verror,1)}else{ErrorSelection<-copy(Verror)} #samples out one error variance 
    # print(ErrorSelection)  
  for (k in 1:length(REPS)){
#..All genotypes in the first rep
    if(REPS[k] ==1){ 
      
#Error  
Vnoise<-(((ErrorSelection)/(1-ErrorSelection))*var(gt[Genotype %in% dt[Env==env[j] & Rep==1]$Genotype]$BreedVal)) #PEV divided by the genetic variance from the additive values
#..noise for the first rep      
        noiseR1<- data.table(Genotype=c(dt[Env==env[j] & Rep==1]$Genotype),Env=c(unique(dt[Env==env[j]]$Env)),Rep=c(1),
                         noise=c(rnorm(nrow(dt[Env==env[j] & Rep==1]), mean = 0, sd = sqrt(Vnoise))))  
NoiseList[[paste0(env[j],'_',REPS[k])]]<-noiseR1
        }else{
          if(REPS[k]==max(REPS[k])){
  Vnoise<-(((ErrorSelection)/(1-ErrorSelection))*var(gt[Genotype %in% dt[Env==env[j] & Rep==REPS[k] & Genotype %in% countsRep[N==REPS[k]]$Genotype]$Genotype]$BreedVal)) #PEV divided by the genetic variance from the additive values
#..All genotypes that are not completely replicated
      NoiseList[[paste0(env[j],'_',REPS[k])]]<- data.table(Genotype=c(dt[Env==env[j] & Rep==REPS[k] & Genotype %in% countsRep[N==REPS[k]]$Genotype]$Genotype),Env=c(unique(dt[Env==env[j]]$Env)),
                                                                  Rep=c(REPS[k]),noise=c(rnorm(nrow(dt[Env==env[j] & Genotype %in% countsRep[N==REPS[k]]$Genotype & Rep==REPS[k]]), mean = 0, sd = sqrt(Vnoise))))
        }else{
Vnoise<-(((ErrorSelection)/(1-ErrorSelection))*var(gt[Genotype %in% dt[Env==env[j] & Rep==REPS[k]]$Genotype]$BreedVal)) #PEV divided by the genetic variance from the additive values
#..All genotypes that are completely replicated
  NoiseList[[paste0(env[j],'_',REPS[k])]]<- data.table(Genotype=c(dt[Env==env[j]]$Genotype),Env=c(unique(dt[Env==env[j]]$Env)),
                                                                 Rep=c(REPS[k]),noise=c(rnorm(nrow(dt[Env==env[j] & Rep==REPS[k]]), mean = 0, sd = sqrt(Vnoise))))
      }
    }
  }
}
noiseFinal<-rbindlist(NoiseList)
#..Combine with original data.table
    FinalDT<-cbind(dt,noiseFinal[,c('noise')])
    
#Creating Random Effects
#..environment    
  ENV<- rnorm(nEnv, mean = 0, sd = sqrt(ENV)) #Environment effect
    names(ENV)<-unique(dt$Env)
#..block variance  
  Block<- rnorm(length(re), mean = 0, sd = sqrt(Block)) #Rep effect    
    names(Block)<-re
#..field  
  Field<- rnorm(nEnv, mean = 0, sd = sqrt(Field)) #Field effect     
    names(Field)<-unique(dt$Env)
        
#Assigning Random Effects to the data.table
  print(paste0('Starting Environment and Field Effect Simulation'))  
  
for (k in 1:nrow(FinalDT)){
#..Simulated environmental effects      
  FinalDT[k,EnvEffect:=ENV[names(ENV)== FinalDT[k,]$Env]]
#..Simualted field effects    
rVal<-Block[names(Block)== FinalDT[k,]$RE]
  if(length(rVal)==0){FinalDT[k,repEffect:=0]
      }else{FinalDT[k,repEffect:=rVal]}
    FinalDT[k,fieldEffect:=Field[names(Field)== FinalDT[k,]$Env]]
#..print progress    
      print(paste0(round(k/nrow(FinalDT)*100, digits = 2),'% Complete'))  
}

# Combine the genetic values with the non-genetic values
  finalvals<- merge(FinalDT,gt, by.x=c('Genotype','Env'),by.y = c('Genotype','Env'))
    
#Generate Final phenotypes 
#..For the genotypic value, subtract off the additive value as the genotypic value is the BV + Vg
  finalvals$PhenoVal<- finalvals$EnvAverage + (finalvals$BreedVal+(finalvals$GenoVal-finalvals$BreedVal)+finalvals$EnvEffect +finalvals$repEffect + finalvals$fieldEffect+finalvals$noise) 
    print(paste0('Simulation Complete!'))
    
#Save the results in a table
return(SimulatedMET=finalvals[,c('Genotype','Env','Rep','PhenoVal')])
}  
