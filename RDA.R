library(vegan)
library(grid)
library(ggplot2)

#Input the sample table
abund_table<-read.csv("30um_CLR_t6sum_CCA.csv", row.names=1, check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)

#Input the metadata 
meta_table<-read.csv("ENV_30_t6ss.csv",row.names=1,check.names=FALSE)

#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table<-meta_table[rownames(abund_table),]

#Filter out any samples taxas that have zero entries 
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table<-subset(meta_table, rowSums(meta_table)!=0)

#Run RDA using all abundances and all environmental variables.
rda_CRL<-rda(abund_table ~ ., data=meta_table)

#create a pdf of the results
#Start by taking the data from the rda
scrs<-scores(rda_CRL,display=c("sp","wa","lc","bp","cn"))
#set up arrows df
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*(multiplier)
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
#set up species df
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")
#set up species colors by genera
df_colors <- read.csv("species_color_codes.csv", row.names = 1)
df_colmatch <- merge(x=df_species, y=df_colors, by="row.names", all.x = TRUE)
#names(df_colmatch)[1] <- "Genus"

#create plot for data from above
p<-ggplot()
#Draw arrows
p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="Black",alpha=1)
#Draw arrow names
df_arrow_names <- df_arrows
p<-p+geom_text(data=as.data.frame(df_arrow_names*1.1),aes(x, y, label = rownames(df_arrows)),color="Black",alpha=0.75)
#Draw Species
p<-p+geom_point(data=df_species,aes(x,y,colour = I(df_colmatch$color)))+scale_shape_manual("",values=2)+labs(colour = "Genera")
#Clean up the plot, make it pretty
p<-p+theme_bw()+theme(panel.grid = element_blank(),
                      text = element_text(size=12, family="sans"))
p<-p+coord_fixed(ratio = 1)
p<-p+xlim(-2,2)
p<-p+ylim(-2,2)
p<-p+geom_hline(yintercept = 0, color = "black")
p<-p+geom_vline(xintercept = 0, color = "black")

#print graph to .eps vector file
postscript("30_sum_RDA_CLR.eps", fonts = "sans")
plot(p)
dev.off()

#Output a text file with the results of ANOVA to check overall significance, axis and term significances, and variance inflation.
sink(file="30_sum_CCA_Sig_BHCorrect.txt", type="output")

rda_CRL
#Percent variance explained by each axis:
100*eigenvals(rda_CRL)/sum(eigenvals(rda_CRL))
anova.cca(rda_CRL, p.adjust.methods="BH")
anova.cca(rda_CRL, by = "axis", p.adjust.methods="BH")
anova.cca(rda_CRL, by = "term", p.adjust.methods="BH")
vif.cca(rda_CRL)

sink()