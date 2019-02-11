#` ---
#` Author Swapna Joshi
#` Title: Marina venn and barplots
#` Date:5262016
#` ---

```{r, message=FALSE} 

library(ggplot2)
library(VennDiagram)

##################### bar plot with color and space
Enrichment<-c(0.87,0.87,0.54, 0, 0.69,0.67,0.33)


KEGG_Pathway <- c("a","b","c","d","e","f","h")
counts <- as.data.frame(cbind(Enrichment,KEGG_Pathway))
counts$KEGG_Pathway <- factor(KEGG_Pathway, 
                       levels=names(rev(sort(table(KEGG_Pathway)))))

# 
counts$val1 <- c("high","high","high", "NA","low","low","low")
rhg_cols <- c(high="#8B0000", low="#CD5C5C", 'NA' = "#e9e9e9")

# dodge <- position_dodge(width = 0.5)
bar <- ggplot(data = counts, aes(x = KEGG_Pathway, y = Enrichment, fill = factor(val1)))  
bar1<-bar +  geom_bar(width = 0.7, stat="identity") + scale_fill_manual(values = rhg_cols) + coord_flip()
# + scale_x_discrete(breaks=counts$KEGG_Pathway[nchar(as.character(counts$KEGG_Pathway))!=1])
ggsave(bar1, file = "C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/phosphorus metabolic process_GO.png")
#+ fig.width=5, fig.height=5
bar1

##################### bar plot with color and space
Enrichment<-c(0.88,  0.86, 0.46, 0.18,0, 0.72, 0.52,0.39)


KEGG_Pathway <- c("a","b","c","d","e","f", "g", "h")
counts <- as.data.frame(cbind(Enrichment,KEGG_Pathway))
counts$KEGG_Pathway <- factor(KEGG_Pathway, 
                              levels=names(rev(sort(table(KEGG_Pathway)))))

# 
counts$val1 <- c("high","high","high","high","NA","low","low", "low")
rhg_cols <- c(high="#191970", low="#6A5ACD", 'NA' = "#e9e9e9")

# dodge <- position_dodge(width = 0.5)
bar <- ggplot(data = counts, aes(x = KEGG_Pathway, y = Enrichment, fill = factor(val1)))  
bar1 <- bar +  geom_bar(width = 0.7, stat="identity") + scale_fill_manual(values = rhg_cols) + coord_flip()
# + scale_x_discrete(breaks=counts$KEGG_Pathway[nchar(as.character(counts$KEGG_Pathway))!=1])
ggsave(bar1, file = "C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/integral to plasma membrane_GO.png" )
#+ fig.width=5, fig.height=5
bar1

# venn diagrams
# grid.newpage()
plot.new()
png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/vennDiagram_chiseq1.png")
draw.pairwise.venn(36897, 28543, 13160, lty = rep("blank",2), fill = c("blue", "red"), cex = c(2,2,2),col = rep("black", 2), lwd = rep(2,2),alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), ext.text =  TRUE)        
dev.off()
#+ fig.width=5, fig.height=5
draw.pairwise.venn(36897, 28543, 13160, lty = rep("blank",2), fill = c("blue", "red"), cex = c(2,2,2),col = rep("black", 2), lwd = rep(2,2),alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), ext.text =  TRUE)        

png("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/vennDiagram_chiseq2.png")
draw.pairwise.venn(19189, 16912, 14837, lty = rep("blank",2), fill = c("blue", "red"), cex = c(2,2,2),col = rep("black", 2), lwd = rep(2,2),alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), ext.text =  TRUE)        
dev.off()
#+ fig.width=5, fig.height=5
draw.pairwise.venn(19189, 16912, 14837, lty = rep("blank",2), fill = c("blue", "red"), cex = c(2,2,2),col = rep("black", 2), lwd = rep(2,2),alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), ext.text =  TRUE)        

# 09-06-2016

Enrichment<-c(0.93,
              0.78,
              0.77,
              0, 1.07,
              0.86,
              0.65,
              0.65
)


KEGG_Pathway <- c("a","b","c","d","e","f", "g", "h")
counts <- as.data.frame(cbind(Enrichment,KEGG_Pathway))
counts$KEGG_Pathway <- factor(KEGG_Pathway, 
                              levels=names(rev(sort(table(KEGG_Pathway)))))

# 
counts$val1 <- c("high","high","high","NA","low","low","low", "low")
rhg_cols <- c(high="#8B0000", low="#CD5C5C", 'NA' = "#e9e9e9")

# dodge <- position_dodge(width = 0.5)
bar <- ggplot(data = counts, aes(x = KEGG_Pathway, y = Enrichment, fill = factor(val1)))  
bar1 <- bar +  geom_bar(width = 0.7, stat="identity") + scale_fill_manual(values = rhg_cols) + coord_flip()
# + scale_x_discrete(breaks=counts$KEGG_Pathway[nchar(as.character(counts$KEGG_Pathway))!=1])
ggsave(bar1, file = "C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/transcription from RNA plot2 promoter.png" )
#+ fig.width=5, fig.height=5
bar1

#################
# Marina analysis 01302017

# heatmap group 1 vs Group 2
# heatmap group 1 vs Group 3


setwd("C:/Users/swapnajoshi/Documents/UCLA_research/Marina_microarray_01082015/marina_analysis_01262017/Swapna_analysis_01302017/")
ge1 <- read.delim("RMA-normalized expression values.txt", sep = "\t", row.names = 1)

de1 <- read.delim("Group 2 vs Group 1 differentially expressed genes.txt", sep = "\t", row.names = 2)
de.gr12 <- de1[de1$p.value.Group.2.vs..Group.1.<0.05, c(2,3,6,8),]; dim(de.gr12)
# [1] 498   5
de.gr13 <- de1[de1$p.value.Group.3.vs..Group.1.<0.05, c(2,3,12,14),]; dim(de.gr13)
# [1] 352   5

ge.de12 <- ge1[row.names(ge1)%in%row.names(de.gr12),c(1:4)]; dim(ge.de12)
# [1] 498   4
ge.de13 <- ge1[row.names(ge1)%in%row.names(de.gr13),c(1:2,5:6)]; dim(ge.de13)
# [1] 352   4
ge12.z<- apply(ge.de12, 1, scale) 
row.names(ge12.z)<-colnames(ge.de12)
ge12.zscore<-t(ge12.z)
# colnames(ge12.zscore)<-substr(colnames(ge12.zscore),3,12)
png("12hm.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.3(as.matrix(ge12.zscore),na.rm=TRUE,scale="none",  
               #		RowSideColor=col.mat,
               #		ColSideColors=col.mat,
               #		col=col.j,
               #		zlim=c(-1.72,1.73),
               Colv=T,Rowv=T,
               cexRow=1,cexCol=1,
               # keep.dendro = NA,
               main = paste("Heatmap", dim(ge12.zscore)[1],"Probes; ",dim(ge12.zscore)[2],"Samples"),
               labRow=row.names(ge12.zscore)
)
dev.off()


