# ==== Assign Simulated individuals to environments in METs given a set of experimental units and speicifc design  ==== #
#' @param NumEU a numeric value corresponding to the overall number of EU to assign to each environment 
#' @param PopSize an optional numeric value corresponding to the overall number of individuals to sample from when assigning genotypes to environments. Default is NULL.
#' @param GenoName an optional vector with length equal to PopSize that contains the genotype names to be assigned to the lines in the MET 
#' @param nEnv a numeric value corresponding to the number of environments in the MET
#' @param nRep a numeric value corresponding to the number of replicates in each environment.  
#' @param percentRep a numeric value between 0 and 1 corresponding to the percent replication of the final replicate. 
#' @param RepAssignment either 'Random' or 'Balanced' (default) if the number of observation per genotype should be even across the total number of environments.
#' @return Data.table containing
#' \describe{
#' \item{SimulatedDT}{A data.table with 3 columns called 'Genotype','Env', and 'Rep' for the simulated genotype name, environments, and replicates.
#' }

library(data.table)
library(dplyr)

SimEnvironments<-function(NumEU,nEnv=2,nRep=2,percentRep=1,PopSize=NULL,GenoName=NULL,RepAssignment='Balanced'){

#Identify the number of individual possible based on EU provided
  NumIndiv<-round(NumEU / (max(nRep) -1 + percentRep))
  
#For scenairos where the population size is not specified   
  if(is.null(PopSize)){PopSize<-NumIndiv}

#The rep assignment argument can be only be Random or Balanced
  if(!RepAssignment %in% c('Random','Balanced')){
    stop('RepAssignment is not equal to Random or Balanced.')
  }
  
#The length of the genotype names has to be equal to population size if provided
  if(!is.null(GenoName) & PopSize !=length(GenoName)){  
    stop("If GenoName is provided, the length must be equal to the total population size.")}  
gName<-GenoName  
#Generate genotype names
  if(!is.null(gName)){  
      SimulatedDT<-data.table(Genotype=c(gName))
    }else{
      Genos<-c(Geno=c(rep(paste0('Geno',seq(1,PopSize,1)),each=1)))
        SimulatedDT<-data.table(Genotype=c(Genos))
}

#Assign Genotypes to fields 
nRep2<-percentRep*(NumIndiv)
EnvList<-list()
EnvNames<-c(Env=c(rep(paste0('Env',seq(1,nEnv,1)),each=1)))
for (i in 1:length(EnvNames)){
      REPS<-seq(1,nRep,1)

#For scenairos where someone accidentily assigned too small of a population for the number of individuals possible given their EU            
if(PopSize < NumIndiv){ #need a population size equal to or greater than that of the number of unique lines per environment.
  stop('Population is smaller than number of individuals given alloted EUs per environemnt. Increase population size!')}
      
#..Subset out the desired number of individuals per environment  
SampleEnv<-sample_n(SimulatedDT,NumIndiv)
SampleEnv$Env<-EnvNames[i]
  for (j in 1:length(REPS)){
    if(REPS[j] !=max(REPS)){
      if(REPS[j]==1){
        SampleEnv$Rep<-REPS[j] #first replicate
          EnvList[[paste0(REPS[j],'_',EnvNames[i])]]<-SampleEnv
            }else{
#..Repeat the desired number of individuals per rep for scenairos where all individuals are repeated
  SampleRep<-copy(SampleEnv)
    SampleRep$Rep<-REPS[j]
EnvList[[paste0(REPS[j],'_',EnvNames[i])]]<-SampleRep
        }
      }
#For the maximum replicate    
    else{
        if(RepAssignment=='Random'){
#..Repeat the desired number of individuals per rep for scenairos where no individuals are repeated for the final rep randomly
  SampleRep<-sample_n(SampleEnv,percentRep*(NumIndiv)) #number of individuals in an environment 
    SampleRep$Rep<-REPS[j]
EnvList[[paste0(REPS[j],'_',EnvNames[i])]]<-SampleRep
#...For scenairos where the replication assignments should be 'Balanced' across environments
        }else{
          if(i==1){
    SampleRep<-sample_n(SampleEnv,percentRep*(NumIndiv)) #number of individuals in an environment 
      SampleRep$Rep<-REPS[j]
#..end the first environment    
  }else{
#Sample for rep 2 after all individuals have been assigned to the second rep    
#..Combine the genotypes and select by minimum plots
    EnvGeno<-rbindlist(EnvList)
     countRem<-EnvGeno[,.N,by=.(Genotype)]
#...order and select the smallest values     
      countRemMin<-countRem[order(countRem$N,decreasing=F)[1:nRep2]]    
        SampleRep<-SampleEnv[Genotype %in% countRemMin$Genotype]
          SampleRep$Rep<-REPS[j]
            }
        SampleRep$Rep<-REPS[j]
          EnvList[[paste0(REPS[j],'_',EnvNames[i])]]<-SampleRep
        }
      } #end the maximum replicate list
    } #end the replicate assigning 
  } #end the environment list
#Combine the data        
  EnvGenoVals<-rbindlist(EnvList)
return(EnvGenoVals)
}  

  