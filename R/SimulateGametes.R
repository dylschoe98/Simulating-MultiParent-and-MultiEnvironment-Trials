#The purpose of the current script is to simualte DH gametes
# ==== Simulate Gametes Derived from Biparental and Multiparent Populations  ==== #
#' @param ParentalNum number of founders to simulate. 
#' @param cM The recombination frequency given as a number between 0.001 and 100 to represent the percent of recombination
#' @param Chrom_Num Number of chromosome for the species. The default is 10 for maize.
#' @param Num_Geno The number of DH to simulate.The default is 500.
#' @param Num_Markers Number of genetic markers to simulate per chromosome. The default is 1,000. The same number of markers are assigned to each chromosome.
#' @return Data.table containing:
#' \describe{
#' \item{BP_Genotype}{ a data.table with the first column corresponding to 'Genotype' for names of the genotypes provided as Geno# (i.e., Geno1) ranging from 1 to the number of genotypes provided in the argument `Num_Geno`}
#' \item{SummaryList}{ a data.table with columns corresponding to 'Genotype', 'Chrom', 'cM_Frequency','cm_Count' to give the genotype name and recombination frequency and count of the number of haplotype switches on a given chromosome. Recombination frequency based on the number of haplotype switches divided by the number of markers provided per chromosome.}
#' }
library(data.table)
library(ggplot2)
library(rlist)


myfunc <- function(x){
  x <- unlist(x)
  values <- table(x)
  max_allele <- which.max(values)
  min_allele <- which.min(values)
  ref_allele <- names(values)[max_allele]
  nonref_allele <- names(values)[min_allele]
  x[x == ref_allele] <- 2
  x[x == nonref_allele] <- 0
  x <- as.numeric(x)
  return(x)
}

SimGametes<-function(ParentalNum,cM,Chrom_Num=10,Num_Geno=500,Num_Markers=1000){

#Number of Potential Parents  
  ParentalNum
gametes<-LETTERS[1:ParentalNum]

#Number of Chromosome
Chrom_Num
chrom<-seq(1,Chrom_Num,1)

#Make a table with genotypes
Num_Geno
Geno_dt<-data.table(Geno=c(rep(paste0('Geno',seq(1,Num_Geno,1)),each=1)))
Genos<-rep(paste0('Geno',seq(1,Num_Geno,1)),each=1)

#Make a table with the markers per chromosome
Num_Markers
Markers<-rep(paste0('Geno',seq(1,Num_Markers,1)),each=1)

#Recombination frequency
cM
Probs<-seq(0.001,100,0.001)

MatrixList<-vector('list',length=max(chrom))
SummaryList<-list()
for (i in 1:length(chrom)){
  genoMat<-matrix(nrow = Num_Geno,ncol = Num_Markers,dimnames = list(c(Genos)))
    colnames(genoMat)<-rep(paste0('M',chrom[i],'_',seq(1,Num_Markers,1)))
  for (j in 1:nrow(genoMat)){ #for the genotypes
     RecombCount<-Randvals<-vector(length = Num_Markers) 
    for (k in 1:ncol(genoMat)){ #for the markers
      
#initial gamete      
  if(k==1){InitialGamete<-sample(gametes,1) 
    genoMat[j,k]<-InitialGamete
  }else{
    
#cM% chance of observing XO at each region
  Recom<-sample(Probs,1) 
    Randvals[k]<-Recom
  
#Identifying a X0 taking place at a given frequency
  if(Recom <=cM){ #for XOs scenairos
RecombCount[k]<-k
  NumRecom<-RecombCount[RecombCount !=0]
  
#Save the average % X0 per genotype
  SummaryList[[paste0(chrom[i],'_',j)]]<-data.table(Genotype=c(Genos[j]),Chrom=c(chrom[i]),
                                                    cM_Frequency=c(length(NumRecom)/Num_Markers*100),cm_Count=c(length(NumRecom)))
#..For scenaiors representing the 1st X0
    if(k==min(NumRecom)){
      XOGamete<-sample(gametes[gametes !=InitialGamete],1)
        genoMat[j,k]<-XOGamete
    }else{
#..For scenaiors representing the 2nd, 3rd, 4th, etc. X0       
    previousGamete<-genoMat[j,k-1]
      genoMat[j,k]<-sample(gametes[gametes !=previousGamete],1)
      }      
#For scenairos when a X0 didn't take place per marker  
  }else{ 
    genoMat[j,k]<- genoMat[j,k-1]
        }
      }#end the else statement for 1st marker
    } #ending the column (i.e., marker) loop
  } #ending the row (i.e., gamete) loop
  print(paste0('Chromosome ',chrom[i], " made!"))  
if(i==1){      
  MatrixList[[i]]<-cbind(data.table(Genotype=c(rownames(genoMat))),as.data.table(genoMat))
    }else{MatrixList[[i]]<-as.data.table(genoMat)
  } #else statement for storing the chromosome
}    
  BP_Genotype<-list.cbind(MatrixList)
return(list(BP_Genotype,rbindlist(SummaryList)))}
 
