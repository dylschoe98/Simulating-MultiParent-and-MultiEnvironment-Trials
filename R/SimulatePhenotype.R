# The purpose of the current script is to simulate breeding values, genotypic values, and phenotypes in the context of a single environment or breeding values and genotypic values in the context of METs.
# ==== Simulate a genetic model from a normal Distribution ==== #
#' @param NumIndiv an optional numeric value for the the number of individuals to simulate. If not provided, 100 individuals are assumed. However, if environmental assignments are provided via the argument 'envDT', the number of genotypes per environment is used and not the 'NumIndiv'.
#' @param envDT an optional data.table as generated from `SimEnvironments` for specifying unique genotypic/additive values per environment in the case of METs with columns provided in order of 'Genotype','Env', and 'Rep' if a uniue genotypic value should be created per environment.
#' @param X a optional numerical marker matrix assigning the major and minor allele. the allele assignments must be numeric values and are recommended to be -1 and 1 based on the theory of realzied additive relationship matrix as described by Endelman and Jannink, 2011. The rownames of the individuals must be spefied and NAs are allowed but imputation is encouraged.
#' @param nQTL a numeric value specifying the number of additive QTL associated with the trait of interest. 
#' @param a numeric value corresponding to the variance for simulating addtive effects. By default, a is 1 as effects are simulated from a mean of 0 with variance of 1.
#' @param Effect a numeric value for simulating the mean of the marker effects. By default, the mean is 0 as effects are simulated from a mean of 0 with variance of 1.
#' @param Poly a numeric value between 0 and 1, representing a percentage, the multiply the smallest QTL effect by to use as the polygenic background effect.
#' @param Vg an optional numeric value for the the percent of non-additive variance to simulate. If not provided, all of the variance is assumed to be due to additive genetic variation.
#' @param u an optional numeric value for the overall mean. If not provided, the mean is assumed to be 0. 
#' @param Verror an optional numeric value between 0 and 1 for the percent of resiudal variance to simualte. If not provided, the genotypic/breeding values will be equal to the phenotypic value in the context of a single environment.   
#' @param QTLnames an optional character vector giving the marker names the QTLs are assigned to. The argument can only be applied if a genetic marker matrix is provided. The argument can only be applied if a genetic marker matrix is provided.
#' @param QTLallele an optional character argument that must be 'Major' or 'Minor' corresponding to if the QTLs are associated with the major or minor allele respectively. if either nQTL, VarQTL, or QTLnames is specified.
#' @param GxE either 'Effect' (default) or 'Marker' for if GxE should be modeled by allowing unique QTL effect estimates per environment or unique QTL (i.e., different markers) and effect estimates 
#' @return Data.table containing
#' \describe{
#' \item{SimulatedGenoMod}{A data.table with a column called 'PhenoVal','BreedVal', and 'GenoVal' for the simulated phenotypes, additive genetics values (i.e., breeding values), and genotypic values respectively if the argument envDT is set to null.. If no genotypic value is provided then the genotypic value and breeding value will be equal. If environment assigments are provided, then the phenotypes must be generated from the `SimMET` function.  
#' }

library(data.table)
library(dplyr)


SimPhenotypes<-function(NumIndiv=100,nQTL,Effect=0,a=1,envDT=NULL,GxE='Effect',u=0,Verror=0,X=NULL,Poly=0,Vg=0,
                          QTLnames=NULL,QTLallele=NULL){

#An internal function to randomly subset out columns of a matrix
SelectMatrix_n <- function(Mat,n){ #https://stackoverflow.com/questions/39441531/get-a-random-column-from-r-matrix
    Cols <- sample(seq_len(ncol(Mat)), size=n)
return(Mat[, Cols, drop=FALSE])}  

#If a marker matrix is provided, the rownames must be provided
if(!is.null(X) & is.null(rownames(X))){
  stop('If a marker matrix is provided, the rownames must be provided to the matrix.')
}

#To specify QTLs, will need to have a marker matrix provided
if( !is.null(QTLnames) & is.null(X) | !is.null(QTLallele) & is.null(X)){
  stop("A genetic marker matrix must be provided to give QTL names.")}
  
# Only allow Effect or Marker as for GxE Modeling
if(!GxE %in% c('Effect','Marker')){
  stop('GxE argument is not recognizable: must be Effect or Marker')
}
  
#If QTLs are provided, specifying if they are associated with the major or minor allele. 
  stopifnot(QTLallele %in% c(NULL,"Major","Minor"))

#Identify the number of environments and create an 
if(is.null(envDT)){
    if(is.null(X)){
  N<-NumIndiv  
    Genos<-c(Geno=c(rep(paste0('Line',seq(1,N,1)),each=1)))
      envDT<-data.table(Genotype=c(Genos),Env=c('Env1'))
        }else{envDT<-data.table(Genotype=c(rownames(X)),Env=c('Env1'))}  
          nEnv<-unique(envDT$Env)
      }else{nEnv<-unique(envDT$Env)
  colnames(envDT)<-c('Genotype','Env','Rep')  
#..Error must be 0 in the context of multiple-environments as phenotypes are generated using the `SimMET` function  
  if(Verror !=0){
    stop('If multiple environments are present, residual error is assigned using the function SimMET')
  }
}  

#To enure the polygenic background effect is a percentage 
  if(Poly >1 | Poly < 0){
    stop('The percentage for the polygenic background effect cannot be less than 0 or greater than 1.')
}

#To specify a polygenic background effect will need to have marker matrix
  if(is.null(X) & Poly !=0){
    stop('To specify a polygenic background effect, a marker matrix must be provided.')
  }
  
#Number of genotypes  
nQTL<-nQTL

#..Genetic value  
  u<-u #overall mean 
  
#Obtain marker effects
  uQTLVal<-rnorm(n = nQTL,mean = Effect,sd = a)
  
#Iterate the process for the number of desired environments 
  SimulatedList<-list()
for (i in 1:length(nEnv)){
    N<-length(unique(envDT[Env ==nEnv[i]]$Genotype))
#Genetic value from Markers  
  if(!is.null(X)){
#..matching genotypes within an environment   
    if(!is.null(envDT)){
      envX<-X[rownames(X) %in% envDT[Env ==nEnv[i]]$Genotype,drop=F,]
#...if only a single environment is provided        
  }else{envX<-copy(X)}
    
#..Number of genetic markers
  NumMarkers<-ncol(envX) - nQTL
  
#..Identifying any markers assigned to QTL   
  if(is.null(QTLnames)){#for scenaiors where QTL are randomly assigned to markers
#...Select if the same markers should be used in each environment     
    if(GxE=='Effect'){
        if(i==1){
#...Only need to subset out the random QTL the first environment if just want to change the effects     
      QTL_cols<-SelectMatrix_n(envX,n = nQTL)
        X_polygenic<-envX[,!colnames(envX) %in% colnames(QTL_cols), drop=FALSE]
          X_QTL<-envX[,colnames(envX) %in% colnames(QTL_cols), drop=FALSE]
              }  
#...Subset out a unique set of QTLs if on the marker and effect level       
          }else{
      QTL_cols<-SelectMatrix_n(envX,n = nQTL)
        X_polygenic<-envX[,!colnames(envX) %in% colnames(QTL_cols), drop=FALSE]
          X_QTL<-envX[,colnames(envX) %in% colnames(QTL_cols), drop=FALSE]
      }
    }else{X_polygenic<-envX[,!colnames(envX) %in% QTL_cols, drop=FALSE] #for scenairos where the QTL are assigned to a specific marker
           X_QTL<-envX[,colnames(envX) %in% QTLnames, drop=FALSE]
}
  if(QTLallele=='Minor'){
#marker coding can be an arbitrary NUMERIC assignment. Therefore, the program needs to identify which allele per locus is the major/minor.
  MinorDesignation<-as.numeric(min(names(table(X_QTL[,1]))))
    MajorDesignation<-as.numeric(max(names(table(X_QTL[,1]))))
#..Change the allele designation    
  X_QTL[X_QTL==MinorDesignation]<-MajorDesignation+1  
    X_QTL[X_QTL==MajorDesignation]<-MinorDesignation-1  
#...Change back to initial starting value designation but opposite 
    X_QTL[X_QTL==MajorDesignation+1]<- (MajorDesignation+1)-1
      X_QTL[X_QTL==MinorDesignation-1 ]<- (MinorDesignation-1)+1
        }
#Calculate the additive values based on the number of QTL
  gebvQTL <- as.vector(crossprod(t(X_QTL),uQTLVal))
    gebvQTL<-mean(gebvQTL) - gebvQTL
#..Marker effects for the polygenic background effects (all small effects of the remaining markers assumed to be equal)
    Vpolygenic<-a*Poly
#....Use the number of remaining markers with a variance equal to the minimum a value divided by the number of markers  
  uPoly <- rnorm(ncol(X_polygenic),mean = Effect,sd = Vpolygenic)
        
#Assign the effect if the value is a 0 or 1
#..Calculate the additive values effects from the polygenic background component
  Polyval <- as.vector(crossprod(t(X_polygenic),uPoly))  
#..Names of the additive values
      names(Polyval)<-rownames(envX)
 
Addval<- Polyval+gebvQTL 
  }else{#If no Marker matrix is provided 
    
# Build a random marker matrix   
M <- matrix(rep(0,N*nQTL),N,nQTL)
  for (m in 1:N){
    M[m,] <- ifelse(runif(nQTL)<0.5,-1,1)}
      Addval <- as.vector(crossprod(t(M),uQTLVal))  
}  
#..error  
  Vnoise<-(((Verror)/1-Verror)*var(Addval)) #PEV divided by the genetic variance from the additive values
    noise<-rnorm(N,mean=0,sd=sqrt(Vnoise))
    
#Correcting for values that could be NULL
genoVal<-Vg

#Genetic effect (dominance, epistasis, SCA)    
if(Vg !=0){genoVal<- rnorm(N, mean = 0, sd = sqrt(genoVal))} 

#Generate Final phenotypes 
  phenoVal<-u + (Addval+genoVal+noise) 
  
#Save the results in a table
if(!is.null(X)){  
  SimulatedList[[i]]<-data.table(Genotype=c(rownames(envX)),Env=c(nEnv[i]),PhenoVal=c(phenoVal),BreedVal=c(Addval),GenoVal=c(Addval+genoVal))
        }else{
    Genos<-c(Geno=c(rep(paste0('Line',seq(1,N,1)),each=1)))
      SimulatedList[[i]]<-data.table(Genotype=c(Genos),Env=c(nEnv[i]),PhenoVal=c(phenoVal),BreedVal=c(Addval),GenoVal=c(Addval+genoVal))
    }  
  } #ends the environment loop
  SimulatedDT<-rbindlist(SimulatedList)
if(Verror==0){SimulatedDT<-SimulatedDT[,-c('PhenoVal')]}  
  if(length(nEnv)==1){SimulatedDT<-SimulatedDT[,-c('Env')]}
  
#Save the results in a table
  if(is.null(X)){return(SimulatedGenoMod=SimulatedDT)}else{return(SimulatedGenoMod=SimulatedDT)}
}




