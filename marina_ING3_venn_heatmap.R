# TODO: Add comment
# 
# Author: SwapnaJoshi
###############################################################################

library(limma)
overExp<-read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/otherPIprojects/Marina/ING3 Overexpr Array.fold 1.csv", sep="," )
row.names(overExp)<-overExp[,1]

save(overExp, file="C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/overExp.rda")


silencing<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/ING3 Silencing Array.Fold 2.csv", sep="," )

row.names(silencing)<-silencing[,1]
save(silencing, file="C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/silencing.rda")


############################ Venn diagram
load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/overExp.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/silencing.rda")

library(VennDiagram)
set1<-row.names(overExp)
set2<-row.names(silencing)

universe <- sort( union(set1, set2) )
Counts <- matrix(0, nrow=length(universe), ncol=2)
colnames(Counts) <- c("ING3 Overexpression", "ING3 Silencing")
for (i in 1:length(universe))

{
	
	Counts[i,1] <- universe[i] %in% set1
	
	Counts[i,2] <- universe[i] %in% set2
	
	
}

vennCounts(Counts)
#   ING3 Overexpression ING3 Silencing Counts
# 1                   0              0      0
# 2                   0              1    567
# 3                   1              0    841
# 4                   1              1     51
# attr(,"class")
# [1] "VennCounts"

plot.new()
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/vennDiagram_counts.png")
draw.pairwise.venn( 618,892,51, c("ING3 Silencing","ING3 Overexpression"), scaled=TRUE, fill=c("blue","red"),
				main = "Venn Diagram")

dev.off()


#51 genes in common

#their labels

common<-sort(intersect (set1, set2) )


silencing.com<-subset(silencing, row.names(silencing)%in%common)

overexpr.com<-subset(overExp, row.names(overExp)%in%common)

write.table(silencing.com, file="C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/silencing.com.csv", sep="\t")

write.table(overexpr.com, file="C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/overexpr.com.csv", sep="\t")



#Heatmap of over expressed


#library(affy)
library("simpleaffy")

setwd("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/otherPIprojects/Marina/raw_data/OverExp/")
OverExpPheno<-read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/otherPIprojects/Marina/raw_data/OverExp/OverExpPheno.csv", sep=",")
SilPheno<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/SilPheno.csv", sep=",")
OverExpPheno<-as.matrix(OverExpPheno)


OverExpDat1 <- read.affy(covdesc="OverExpPheno.txt")   
OverExpDat1


setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/SILENCING/")

SilDat1 <- read.affy(covdesc="SilPheno.txt")   
SilDat1


setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/OverExp/")


rmaOver <- rma(OverExpDat1)

#Background correcting
#Normalizing
#Calculating Expression

rmaOver
# load colour libraries
library(RColorBrewer)
# set colour palette
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/boxplot_qc.png")		
boxplot(rmaOver, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.rma
library(affyPLM)
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/boxplot_normalized_qc.png")		 
boxplot(rmaOver, col=cols)
dev.off()
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/den.vs.log.intensity.histogram_qc.png")		
hist(rmaOver, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/den.vs.log.intensity.histogram_norm_qc.png")		
hist(rmaOver, col=cols)
dev.off()

#From these plots we can conclude that there are no major deviations amongst the 12 chips, and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. To take a closer look at the situation on a per-chip level we can use affyPLM. affyPLM allows us to visualise statistical characteristics of the CEL files.

# Perform probe-level metric calculations on the CEL files:
celfiles.qc <- fitPLM(OverExpDat1) 


# Create an image of GSM24662.CEL:
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/image.cel.png")	
image(celfiles.qc, which=1, add.legend=TRUE)
dev.off()

# # Create an image of GSM524665.CEL
#		 # There is a spatial artifact present
#		 image(celfiles.qc, which=4, add.legend=TRUE)
# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero. GSM524665.CEL is an outlier
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/rle.image.png")	

RLE(celfiles.qc, main="RLE")

dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
# The median standard error should be 1 for most genes.
# GSM524665.CEL appears to be an outlier on this plot too
png("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/image.nuse.png")	

NUSE(celfiles.qc, main="NUSE")
dev.off()


#We can also look at the relationships between the samples using heirarchical clustering:


eset <- exprs(rmaOver)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/clusters.norm_overExp.png")	
plot(clusters)
dev.off()

#Filtering data
load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/celfiles.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/celfiles.rma.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/eset.rda")

#Now we have looked at the data, we can go on to analyse it. The first stage of analysis is to filter out uninformative data such as control probesets and other internal controls as well as removing genes with low variance, that would be unlikely to pass statistical tests for differential expression, or are expressed uniformly close to background detection levels. The modifiers to nsFilter below tell nsFilter not to remove probesets without Entrez Gene identifiers, or have duplicated Entrez Gene identifiers.


celfiles.filtered <- nsFilter(rmaOver, require.entrez=FALSE, remove.dupEntrez=FALSE)


# What got removed and why?
celfiles.filtered$filter.log$numLowVar
# [1] 27307


celfiles.filtered$filter.log$feature.exclude
# [1] 62




#Heatmap


load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/overExp.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/silencing.rda")
save(exp.df.over, file = "marina.OverExpr.cell12.rda")

exp.df.over<-eset[row.names(eset)%in%row.names(overExp),]

dim(exp.df.over)
# [1] 892   4
exp.df.over<-exp.df.over[row.names(overExp),]
match(row.names(exp.df.over), row.names(overExp))
col.j<-jet.colors(75)

overExpZ<- apply(exp.df.over, 1, scale) 
row.names(overExpZ)<-colnames(exp.df.over)
overExpZ.zscore<-t(overExpZ)

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/overExp.hm.zscore1.hm3.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(overExpZ.zscore),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
#		col=col.j,
#		zlim=c(0,1),
		Colv=T,Rowv=T,
		cexRow=1,cexCol=1,
		# keep.dendro = NA,
		main = paste("IBS", dim(overExpZ.zscore)[1],"Probes; ",dim(overExpZ.zscore)[2],"Samples"),
		labRow=row.names(overExpZ.zscore)
)
dev.off()

#which(overExp[,2]%in%"ING3")
# [1] 456

row.order <- hv2$rowInd
overExpZ.zscore1<-overExpZ.zscore[row.order,]
which(row.names(overExpZ.zscore1)%in%"205070_at")
# [1] 132


ing3flank<-overExpZ.zscore1[130:135,]



png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/ing3flank-over1-hm3.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(ing3flank),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
#		col=col.j,
#		zlim=c(0,1),
		Colv=NA,Rowv=NA,
		cexRow=1,cexCol=1,
		# keep.dendro = NA,
		main = paste("IBS", dim(ing3flank)[1],"Probes; ",dim(ing3flank)[2],"Samples"),
		labRow=row.names(ing3flank)
)
dev.off()


#######################################################################################################
#Heatmap of silenced

#library(affy)
library("simpleaffy")

setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/Silencing/")
OverExpPheno<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/OverExp/OverExpPheno.csv", sep=",")
SilPheno<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/Silencing/SilPheno.csv", sep=",")
OverExpPheno<-as.matrix(OverExpPheno)
SilPheno<-as.matrix(SilPheno)

SilDat1 <- read.affy(covdesc="SilPheno.txt")   
SilDat1


setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/SILENCING/")

SilDat1 <- read.affy(covdesc="SilPheno.txt")   
SilDat1


setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/silencing/")


rmaSil <- rma(SilDat1)

#Background correcting
#Normalizing
#Calculating Expression

rmaSil
# load colour libraries
library(RColorBrewer)
# set colour palette
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/boxplot_qc_Sil.png")		
boxplot(rmaSil, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.rma
library(affyPLM)
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/boxplot_normalized_qc_sil.png")		 
boxplot(rmaSil, col=cols)
dev.off()
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/den.vs.log.intensity.histogram_qc_sil.png")		
hist(rmaSil, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/den.vs.log.intensity.histogram_norm_qc_sil.png")		
hist(rmaSil, col=cols)
dev.off()

#From these plots we can conclude that there are no major deviations amongst the 12 chips, and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. To take a closer look at the situation on a per-chip level we can use affyPLM. affyPLM allows us to visualise statistical characteristics of the CEL files.

# Perform probe-level metric calculations on the CEL files:
celfiles.qc <- fitPLM(SilDat1) 


# Create an image of GSM24662.CEL:
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/image.cel.sil.png")	
image(celfiles.qc, which=1, add.legend=TRUE)
dev.off()

# # Create an image of GSM524665.CEL
#		 # There is a spatial artifact present
#		 image(celfiles.qc, which=4, add.legend=TRUE)
# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero. GSM524665.CEL is an outlier
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/rle.image.sil.png")	

RLE(celfiles.qc, main="RLE")

dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
# The median standard error should be 1 for most genes.
# GSM524665.CEL appears to be an outlier on this plot too
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/image.nuse.sil.png")	

NUSE(celfiles.qc, main="NUSE")
dev.off()


#We can also look at the relationships between the samples using heirarchical clustering:


eset <- exprs(rmaSil)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/clusters.norm_Sil.png")	
plot(clusters)
dev.off()

#Filtering data
#load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/celfiles.rda")
#load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/celfiles.rma.rda")
#load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/eset.rda")

#Now we have looked at the data, we can go on to analyse it. The first stage of analysis is to filter out uninformative data such as control probesets and other internal controls as well as removing genes with low variance, that would be unlikely to pass statistical tests for differential expression, or are expressed uniformly close to background detection levels. The modifiers to nsFilter below tell nsFilter not to remove probesets without Entrez Gene identifiers, or have duplicated Entrez Gene identifiers.


celfiles.filtered <- nsFilter(rmaSil, require.entrez=FALSE, remove.dupEntrez=FALSE)


# What got removed and why?
celfiles.filtered$filter.log$numLowVar
# [1] 27307


celfiles.filtered$filter.log$feature.exclude
# [1] 62



save(eset, file="Ing.kmt2d.silencing.eset.rda")

#Heatmap


load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/OverExp.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/silencing.rda")
row.names(silencing)<-silencing[,1]


exp.df.Sil<-eset[row.names(eset)%in%row.names(silencing),]

dim(exp.df.Sil)
# [1] 618   6

exp.df.Sil<-exp.df.Sil[row.names(silencing),]
match(row.names(exp.df.Sil), row.names(silencing))
col.j<-jet.colors(75)

exp.df.Sil.ing3<-exp.df.Sil[,-c(3,4)]

SilZ<- apply(exp.df.Sil.ing3, 1, scale) 
row.names(SilZ)<-colnames(exp.df.Sil.ing3)
SilZ.zscore<-t(SilZ)

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/Sil.ing3.hm.zscore.hm3.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(SilZ.zscore),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
#		col=col.j,
#		zlim=c(0,1),
		Colv=T,Rowv=T,
		cexRow=1,cexCol=1,
		# keep.dendro = NA,
		main = paste("IBS", dim(SilZ.zscore)[1],"Probes; ",dim(SilZ.zscore)[2],"Samples"),
		labRow=row.names(SilZ.zscore)
)
dev.off()

#which(Sil[,2]%in%"ING3")
# [1] 456

row.order <- hv2$rowInd
SilZ.zscore1<-SilZ.zscore[row.order,]
silencing[silencing[,2]%in%"ING3",]

which(row.names(SilZ.zscore1)%in%"205070_at")
# [1] 129


ing3flank<-SilZ.zscore1[125:135,]




png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/ing3flank-hm3-silencing-ing3.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(ing3flank),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
#		col=col.j,
#		zlim=c(0,1),
		Colv=NA,Rowv=NA,
		cexRow=1,cexCol=1,
		# keep.dendro = NA,
		main = paste("IBS", dim(ing3flank)[1],"Probes; ",dim(ing3flank)[2],"Samples"),
		labRow=row.names(ing3flank)
)
dev.off()


#############################

kmt2d.sil<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/sample2_vs_sample1_p0.05_12785genes_for kmt2d.csv", sep=",")
load("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/raw_data/Silencing/Ing.kmt2d.silencing.eset.rda")
row.names(kmt2d.sil)<-kmt2d.sil[,1]

kmt2d.sil1<-subset(kmt2d.sil, kmt2d.sil$p.value.Sample2.vs..Sample1.<0.05 & kmt2d.sil$Fold.Change.Sample2.vs..Sample1.<=-2 |  kmt2d.sil$p.value.Sample2.vs..Sample1.<0.05 & kmt2d.sil$Fold.Change.Sample2.vs..Sample1.>=2)

exp.df.Sil<-eset[row.names(eset)%in%row.names(kmt2d.sil1),]

dim(exp.df.Sil)
# [1] 8845    6


exp.df.Sil<-exp.df.Sil[row.names(kmt2d.sil1),]
match(row.names(exp.df.Sil), row.names(kmt2d.sil1))
col.j<-jet.colors(75)

exp.df.Sil.kmt2d<-exp.df.Sil[,-c(5,6)]

SilZ<- apply(exp.df.Sil.kmt2d, 1, scale) 
row.names(SilZ)<-colnames(exp.df.Sil.kmt2d)
SilZ.zscore<-t(SilZ)


png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/Sil.ktm2d.hm.z.1.25cutoff_rep.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(SilZ.zscore),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
#		col=col.j,
#		zlim=c(0,1),
		Colv=T,Rowv=T,
		cexRow=1,cexCol=1,
		# keep.dendro = NA,
		main = paste("IBS", dim(SilZ.zscore)[1],"Probes; ",dim(SilZ.zscore)[2],"Samples"),
		labRow=row.names(SilZ.zscore)
)
dev.off()

#which(Sil[,2]%in%"ktm2d")
# [1] 456

row.order <- hv2$rowInd
SilZ.zscore1<-SilZ.zscore[row.order,]
row.names(kmt2d.sil1[kmt2d.sil1[,2]%in%"KMT2D",])
# [1] "227528_s_at" "227527_at"  


which(row.names(SilZ.zscore1)=="227528_s_at")
# [1] 1160
which(row.names(SilZ.zscore1)=="227527_at")
# [1] 4775




ktm2dflank<-SilZ.zscore1[4774:4776,]




png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/ktm2dflank-hm-silencing_227527_at.cut0ff1.25.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(ktm2dflank),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
#		col=col.j,
#		zlim=c(0,1),
		Colv=NA,Rowv=NA,
		cexRow=1,cexCol=1,
		# keep.dendro = NA,
		main = paste("IBS", dim(ktm2dflank)[1],"Probes; ",dim(ktm2dflank)[2],"Samples"),
		labRow=row.names(ktm2dflank)
)
dev.off()


save(kmt2d.sil1, file="kmt2d.sil.1.25.fc.rda")
write.table(kmt2d.sil1,file="C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/kmt2d.sil.1.25.fc.csv", sep=",", col.names=NA) 
write.table(kmt2d.sil1,file="C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/kmt2d.sil.2.fc.csv", sep=",", col.names=NA) 

		
##################################################### venn diagram
over_ing3<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/new_venn_diagram_1.5FC_ing_over_sil/ING3 overexpression fold change 1.5.csv", sep=",")
sil_ing3<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/new_venn_diagram_1.5FC_ing_over_sil/ING3 Silencing.fold change 1.csv", sep=",")
sil_ing3_1p5<-subset(sil_ing3, sil_ing3$Fold.Change.Sample3.vs..Sample1.<=-1.5|sil_ing3$Fold.Change.Sample3.vs..Sample1.>=1.5)

#row.names(over_ing3)<-over_ing3[,1]

library(VennDiagram)
set0<-over_ing3[,1]
set0<-gsub("---",NA, set0)
set1<-na.omit(set0)
set01<-sil_ing3_1p5[,2]
set01<-gsub("---",NA, set01)
set2<-na.omit(set01)
universe1 <- sort( union(set1, set2) )
Counts <- matrix(0, nrow=length(universe1), ncol=2)
colnames(Counts) <- c("ING3 Overexpression", "ING3 Silencing")
for (i in 1:length(universe1))

{
	
	Counts[i,1] <- universe1[i] %in% set1
	
	Counts[i,2] <- universe1[i] %in% set2
	
	
}
Counts<-as.data.frame(Counts)
Counts$synbol<-universe1
overlap<-subset(Counts, Counts$"ING3 Overexpression"==1 & Counts$"ING3 Silencing"==1)
overlap$synbol

library(limma)

vennCounts(Counts)

library(colorfulVennPlot)
#proportions: (1841/1946)*100; (71/1946)*100; (34/1946)*100; 
# [1] 94.60432
# [1] 3.64851
# [1] 1.747174




png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/results/venn_ing_over_sil_new.png",height = 4, width = 4, units = 'in',res = 400)
plotVenn2d(rep("",3), radius=resizeCircles(94.6, 3.6, 1.7), Title=NULL, Colors=c("#99CCFF","pink","purple"),resizePlot=0.7
		, labels=c("",""))
dev.off()


# Marina heatmap

#methylation data
load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/nor1.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/can.rda")

meth.dat<-cbind(nor1,can)

methylation.probes<-c("cg05743713","cg15234492","cg01878308","cg00522588","cg13007988","cg16694837")

sel.pr.meth.can<-can[row.names(can)%in%methylation.probes,]
sel.pr.meth.nor<-nor1[row.names(nor1)%in%methylation.probes,]


png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/can_heatmap_methylation_histone_modifiers.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(sel.pr.meth.can),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
		col=col.j,
		zlim=c(0,1),
		Colv=T,Rowv=T,
		cexRow=1,cexCol=1,
		keep.dendro = T,
		main = paste("PanCancerMethylation", dim(sel.pr.meth.can)[1],"Probes; ",dim(sel.pr.meth.can)[2],"Samples"),
		labRow=row.names(sel.pr.meth.can)
)
dev.off()

row.order <- hv2$rowInd
nor<-sel.pr.meth.nor[row.order,]


png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/nor_heatmap_methylation_histone_modifiers.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(nor),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
		col=col.j,
		zlim=c(0,1),
		Colv=T,Rowv=NA,
		cexRow=1,cexCol=1,
		keep.dendro = T,
		main = paste("PanCancerMethylation", dim(nor)[1],"Probes; ",dim(nor)[2],"Samples"),
		labRow=row.names(nor)
)
dev.off()




#expression data 
exp.pca<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/expr.mat.ann.csv", sep=",")
row.names(exp.pca)<-exp.pca[,1]

mRNA.probes<-c("208988_at","226215_s_at","212689_s_at","227527_at","227527_at","229940_at")

sel.pr<-exp.pca[row.names(exp.pca)%in%mRNA.probes,]

#based on the order of methylation probes rearrage the order of expression probes
#after reorganizing, the order is
#cg16694837,cg00522588,cg13007988,cg01878308,cg15234492,cg05743713
#229940_at,227527_at,227527_at,212689_s_at,226215_s_at,208988_at
row.names(sel.pr)
# [1] "208988_at"   "212689_s_at" "226215_s_at" "227527_at"   "229940_at"  

sel.pr<-sel.pr[c(1,3,2,4,4,5),c(2:23)]

sel.pr1<-t(apply(sel.pr,1,scale))
colnames(sel.pr1)<-colnames(sel.pr)

sel.pr.can<-sel.pr1[,grep("Ca",colnames(sel.pr1))]
sel.pr.nor<-sel.pr1[,!colnames(sel.pr1)%in%colnames(sel.pr.can)]
col.j<-jet.colors(75)

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/can_heatmap_expression_histone_modifiers.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(sel.pr.can),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
		col=col.j,
#		zlim=c(0,1),
		Colv=T,Rowv=NA,
		cexRow=1,cexCol=1,
		keep.dendro = T,
		main = paste("PanCancerExpression", dim(sel.pr.can)[1],"Probes; ",dim(sel.pr.can)[2],"Samples"),
		labRow=row.names(sel.pr.can)
)
dev.off()

row.order <- hv2$rowInd
nor<-sel.pr.nor[row.order,]


png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/nor_heatmap_expression_histone_modifiers.png")#, bg="white", res=300, width=3000, height=3000)
hv3<-heatmap.plus(as.matrix(sel.pr.nor),na.rm=TRUE,scale="none",  
#		RowSideColor=col.mat,
#		ColSideColors=col.mat,
		col=col.j,
#		zlim=c(0,1),
		Colv=T,Rowv=NA,
		cexRow=1,cexCol=1,
		keep.dendro = T,
		main = paste("PanCancerExpression", dim(sel.pr.nor)[1],"Probes; ",dim(sel.pr.nor)[2],"Samples"),
		labRow=row.names(sel.pr.nor)
)
dev.off()



#KMT2d
#Pseudo code:
#		1.import cell 2a, cell 2b data (kmt2d) import controls 1a and 2b
#2.import fold change data 
#3.import raw data for expression
#4.select data for 1.25 cut off
#5.make a heatmap

#import fold change data for kmt2d
kmt2d.sil.de<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/sample2_vs_sample1_p0.05_12785genes_for kmt2d.csv", sep=",", check.names=FALSE)

dim(kmt2d.sil.de)
head(kmt2d.sil.de)

kmt2d.sil.de.1.25<-subset(kmt2d.sil.de, kmt2d.sil.de$'Fold-Change(Sample2 vs. Sample1)')



# chromatin modifier heatmap 02012016


panExp<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/Expression_data/exp.analysis/PanCancerGeneExpression.csv", sep=",")
row.names(panExp)<-panExp[,1]

probes1<-c("227527_at","208988_at","1556493_a_at","201548_s_at","220070_at","235339_at","1554555_a_at","219751_at","222566_at")

sel.pr<-panExp[row.names(panExp)%in%probes1,c(2:23)]


sel.pr1<-t(apply(sel.pr,1,scale))
colnames(sel.pr1)<-colnames(sel.pr)

sel.pr.can<-sel.pr1[,grep("Ca",colnames(sel.pr1))]
sel.pr.nor<-sel.pr1[,!colnames(sel.pr1)%in%colnames(sel.pr.can)]

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/can_heatmap_expression_histone_modifiers_new.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(sel.pr.can),na.rm=TRUE,scale="none",  
                  #		RowSideColor=col.mat,
                  #		ColSideColors=col.mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=T,Rowv=T,
                  cexRow=1,cexCol=1,
                  keep.dendro = T,
                  main = paste("PanCancerExpression", dim(sel.pr.can)[1],"Probes; ",dim(sel.pr.can)[2],"Samples"),
                  labRow=row.names(sel.pr.can)
)
dev.off()

row.order <- hv2$rowInd
nor<-sel.pr.nor[row.order,]


png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/nor_heatmap_expression_histone_modifiers_new.png")#, bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(nor),na.rm=TRUE,scale="none",  
                  #		RowSideColor=col.mat,
                  #		ColSideColors=col.mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=T,Rowv=NA,
                  cexRow=1,cexCol=1,
                  keep.dendro = T,
                  main = paste("PanCancerExpression", dim(nor)[1],"Probes; ",dim(nor)[2],"Samples"),
                  labRow=row.names(nor)
)
dev.off()



png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/colorKey.png",bg="white")
my.colors = jet.colors(75)
z=matrix(1:75,nrow=1)
x=1
y=seq(-3.43,3.43,len=75) # supposing 3 and 2345 are the range of your data
image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")

axis(2)
dev.off()

