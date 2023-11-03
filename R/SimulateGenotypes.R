# ==== Simulate genotypes based on parental genotypes or randomly smapling a given parental state per locus  ==== #
#' @param Gametes a data.table of simulated gametes by markers with cells corresponding to the parental value at each locus per simulated progeny as obtained from `SimGametes`. The individuals must be spefied and NAs are not allowed.
#' @param ParentMat a numerical marker matrix of parental genotypes that has 1's and 0 where 1 is the major allele and 0 is the minor allele. The rownames of the individuals must be spefied and NAs are not allowed.
#' @param nParents a numerical value for the number of parents associated with the population. Must be set to NULL if the parent matrix is provided.
#' @param MAF a numerical vector of potential MAF at each locus.


SimGenotypes<-function(Gametes,ParentMat=NULL,nParents=2,MAF=NULL){        
  gameteVals<-as.matrix(Gametes[,-c(1)])
    rownames(gameteVals)<-Gametes$Genotype
    
#If the parent matrix is provided
if(!is.null(ParentMat)){   
#If the nParents is provided with the parent matrix fail
  if(!is.null(nParents)){
    stop('Cannot provide both number of parents and a parental matrix.')
  }
#If the gamete and the parental matrix is present, they must have the same number of columns 
  if(!is.null(ParentMat) & ncol(gameteVals) != ncol(ParentMat)){
    stop('Number of markers in the parents and gametes must be equal.')
  }
#If the rownames of parents is not provided
  if(is.null(rownames(ParentMat))){
    stop('Rownames of the parent matrix must be provided.')
  }
}
    
#Matrix to store the results
  gameteMat<-matrix(nrow = nrow(gameteVals),ncol=ncol(gameteVals),
                    dimnames = list(c(rownames(gameteVals)),c(colnames(gameteVals))))
#..Generating numerical values
  numericMat<-matrix(nrow = nrow(gameteVals),ncol=ncol(gameteVals),
                   dimnames = list(c(rownames(gameteVals)),c(colnames(gameteVals))))

#For scenairos where a parental matrix is not provided
if(is.null(ParentMat)){  
  InfoList<-list()
  for (i in 1:ncol(gameteMat)){
    gameteSub<-gameteVals[,i]
    
  if(nParents==2){AllelePair<-c('A/T','C/G')} else{ 
    AllelePair<-c('A/C','A/G','A/T','C/G','T/C','T/G')
}
  
#Sample a pair of alleles
  pairSamp<-sample(AllelePair,1)  
#..split the allele pair into two 
    Alleles<-unlist(strsplit(pairSamp,"/"))
      MINOR<-sample(Alleles,1)
      
#Sample a haplotype at random
  hapSamp<-sample(gameteSub,1)  
  
#Sample for haplotype sharing
# sharePool<-c(seq(1/nParents),1)
  shareSamp<-sample(sharePool,1)
  
if(shareSamp=='Yes'){
  shareHapSamp<-sample(gameteSub[gameteSub !=hapSamp],1)
    hapSampShare<-c(hapSamp,shareHapSamp)
}else{hapSampShare<-hapSamp
  shareHapSamp<-'No'
} 
#Sample the major allele and save the results
  MarkerInfo<-data.table(Marker=c(colnames(gameteMat)[i]), 
                      AllelePair=c(pairSamp),haplotype=c(hapSamp),MinorAllele=c(MINOR),Share=c(shareHapSamp))
InfoList[[i]]<-MarkerInfo
  
#Assign the haplotypes to the allele designation
gameteSub[gameteSub %in% hapSampShare]<-MINOR
  gameteSub[gameteSub !=MINOR]<-Alleles[Alleles !=MINOR]
#..add to the matrix
  gameteMat[,i]<-gameteSub
    
#Convert the allele designation to numeric matrix
gameteSub[gameteSub==MINOR]<-0 #minor allele
  gameteSub[gameteSub !="0"]<-1 #major allele
#..add to the matrix
    numericMat[,i]<-gameteSub
  } #end initial forloop
# For the presence of the parental matrix
}else{
  for (i in 1:ncol(gameteVals)){
    for (j in 1:nrow(gameteVals)){
      Parent<-gameteVals[j,i]
        gameteVals[j,i]<- as.numeric(parMat[,i][names(parMat[,i])== Parent])
      }
    print(i/ncol(gameteVals)*100)  
    } 
#Save the Numerical and the genetic objects with the genotype names
    NumDT<-cbind(Gametes[,c(1)],as.data.table(gameteVals)) 
  }
#If the parent matrix is provided
if(is.null(ParentMat)){   
#Save the Numerical and the genetic objects with the genotype names
  NumDT<-cbind(Gametes[,c(1)],as.data.table(numericMat)) 
    GeneticDT<-cbind(Gametes[,c(1)],as.data.table(gameteMat)) 
      InfoDT<-rbindlist(InfoList)
#Returns all information if a parental matrix is not provided    
return(GenotypeList=list(NumDT,GeneticDT,InfoDT))
  }else{return(NumDT)}  
}

