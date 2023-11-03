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


# #MAGIC Populations
# MAGIC_GenotypeList<-SimGametes(ParentalNum = 6,cM = 0.3,Num_Markers = 1000) # The Chrom_Num = 10,Num_Geno = 500,Num_Markers = 1000 are default
# MAGIC_Genotype<-MAGIC_GenotypeList[[1]]
# MAGIC_Genotype<-fread("Results/MAGICparentalGenotypes.csv")
#   fwrite(MAGIC_Genotype,"Results/MAGICparentalGenotypes.csv")
# #..Save a parental matrix with designation from WI-SS-MAGIC
#   MAGIC_Genotype<-fread("Results/MAGICparentalGenotypes.csv")
# #..Convert names in gametes to parent names 
#   MAGIC_Genotype[MAGIC_Genotype=='A']<-'B73'
#     MAGIC_Genotype[MAGIC_Genotype=='B']<-'B84'
#       MAGIC_Genotype[MAGIC_Genotype=='C']<-'LH145'
# MAGIC_Genotype[MAGIC_Genotype=='D']<-'NKH8431'
#   MAGIC_Genotype[MAGIC_Genotype=='E']<-'PHB47'
#     MAGIC_Genotype[MAGIC_Genotype=='F']<-'PHJ40'
# #....Save copies for next step
#     fwrite(MAGIC_Genotype,"Results/MAGICparentalGenotypesWI-SS-MAGIC.csv")
#       fwrite(MAGIC_Genotype,"~/Simulations/G2Fsimulation/2-GenotypeGeneration/Data/MAGICparentalGenotypesWI-SS-MAGIC.csv")
#         
# #..Average number of X0 in MAGIC
#   MAGIC_summary<-MAGIC_GenotypeList[[2]]
#     fwrite(MAGIC_summary,"Results/MAGICGenotypesSummary.csv")
# sumXO_MAGIC<-MAGIC_summary[,lapply(.SD,sum), by=.(Genotype),.SDcols=c(4)]
# 
# #..Format the data
#   melt_magic<-melt(MAGIC_Genotype,id.vars = 'Genotype',measure.vars = 2:ncol(MAGIC_Genotype))
#     melt_magic$chr<-as.numeric(gsub("M(.+)_.*", "\\1", melt_magic$variable))
# melt_magic$pos<-as.numeric(gsub(".*_", "", melt_magic$variable))
#   setnames(melt_magic,c('value'),c('Haplotype'))
#   
# #MAGIC Populations V2
# MAGIC_GenotypeList_v2<-SimGametes(ParentalNum = 6,cM = 0.5,Num_Markers = 1000) # The Chrom_Num = 10,Num_Geno = 500,Num_Markers = 1000 are default
#   MAGIC_GenotypeV2<-MAGIC_GenotypeList_v2[[1]]
#   fwrite(MAGIC_GenotypeV2,"Results/MAGICv2parentalGenotypes.csv")
# #..Average number of X0 in MAGIC
#   MAGICv2_summary<-MAGIC_GenotypeList_v2[[2]]
#     fwrite(MAGICv2_summary,"Results/MAGICGv2enotypesSummary.csv")
# sumXO_MAGICv2<-MAGICv2_summary[,lapply(.SD,sum), by=.(Genotype),.SDcols=c(4)]
#   
# #..Format the data
#   melt_magic2<-melt(MAGIC_GenotypeV2,id.vars = 'Genotype',measure.vars = 2:ncol(MAGIC_Genotype))
#     melt_magic2$chr<-as.numeric(gsub("M(.+)_.*", "\\1", melt_magic2$variable))
#     melt_magic2$pos<-as.numeric(gsub(".*_", "", melt_magic2$variable))
#       setnames(melt_magic2,c('value'),c('Haplotype'))
#       
# #Plot Asthetics
# cbPallete<-c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#CC6677", "#DDCC77", 
#                      "#117733", "#332288","#44AA99", "#999933", "#882255","#999999","#D55E00", "#CC79A7","#661100", "#6699CC", 
#                           "#888888","#999999")
# themeGamete<-function(x) theme(plot.title=element_text(hjust=0.5,size=18,face="bold"),strip.text.y = element_text(size=17),strip.text.x = element_text(size=16), 
#                                   axis.text.x = element_text(angle=45, hjust=1,size=15),axis.text.y = element_text(angle = 0, hjust = 1,size = 3),
#                                         plot.subtitle = element_text(face = "italic",size = 16),axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18))
# #Plot Biparental Population
# ggplot(melt_bp[chr==1],aes(x=pos,y=Genotype))+geom_tile(aes(x=pos,y=Genotype,fill=Haplotype)) +theme_bw()+themeGamete()+
#   ggtitle('Simulated Biparental Double Haploids',subtitle = 'Chromosome 1 Shown') + xlab('Position') + ylab("Genotype") +
#     scale_fill_manual(values = cbPallete)+scale_x_continuous(expand = c(0, 0),breaks=seq(0,1000,100))
# ggsave("Results/SimulatedDH_Biparental_XO.png",width = 10,height = 15)  
# 
# #Plot MAGIC Population
# ggplot(melt_magic[chr==1],aes(x=pos,y=Genotype))+geom_tile(aes(x=pos,y=Genotype,fill=Haplotype)) +theme_bw()+themeGamete()+
#   ggtitle('Simulated MAGIC Double Haploids',subtitle = 'Chromosome 1 Shown') + xlab('Position') + ylab("Genotype") +
#     scale_fill_manual(values = cbPallete)+scale_x_continuous(expand = c(0, 0),breaks=seq(0,1000,100)) 
# ggsave("Results/SimulatedDH_MAGIC_XO.png",width = 10,height = 20)  
# 
# #..Plot the Genotypes
# genos<-unique(melt_magic$Genotype)
# themeGameteSolo<-function(x) theme(plot.title=element_text(hjust=0.5,size=18,face="bold"),strip.text.y = element_text(size=17),strip.text.x = element_text(size=16), 
#                                axis.text.x = element_text(angle=45, hjust=1,size=15),axis.text.y = element_text(angle = 0, hjust = 1,size = 15),
#                                     plot.subtitle = element_text(face = "italic",size = 16),axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18))
# 
# for (i in 1:length(genos)){
#   geno_sub<-melt_magic[Genotype==genos[i]]
#     geno_sub$chr<-as.factor(geno_sub$chr)
# #Export the plot  
#   ggplot(geno_sub,aes(x=chr,y=pos))+geom_tile(aes(x=chr,y=pos,fill=Haplotype)) +theme_classic()+themeGameteSolo()+
#     ggtitle(paste0(geno_sub[i],': ',sumXO_MAGIC[Genotype==genos[i]]$cm_Count,' Total Crossovers'),subtitle = 'Simulated MAGIC Double Haploids') + xlab('Chromosome') + ylab("Position") +
#     scale_fill_manual(values = cbPallete)+scale_y_continuous(expand = c(0, 0),breaks=seq(0,1000,100)) 
# ggsave(filename = paste0("Results/GameteHaplotypes/MAGIC_",genos[i],'.png'),width = 8,height = 8)
# }
# 
# save.image("SavePoints/1-SimulatedGametes.RData")
# 
# 
# 
# 
# 
# 
#   
# 
#   
#   
# 
# 
