#
#Code to perform CCA from metadata file and sample file for significant ENV variables.
#

library(vegan)
library(grid)
library(ggplot2)

#Input the sample table
abund_table<-read.csv("Seq_reads.csv",row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)

#Input the metadata 
meta_table<-read.csv("Envs3.csv",row.names=1,check.names=FALSE)

#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table<-meta_table[rownames(abund_table),]

#Filter out any samples taxas that have zero entries 
abund_table<-subset(abund_table,rowSums(abund_table)!=0)

#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)

#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data=meta_table)

#Extract significant and non-significant variables
bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.05]
otherEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)">0.05]

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
otherEnvVariables<-otherEnvVariables[!is.na(otherEnvVariables)]

#We are now going to separately use those environmental variables in cca that were found significant and non-significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=meta_table)",sep="")))
eval(parse(text=paste("sol2 <- cca(abund_table ~ ",do.call(paste,c(as.list(otherEnvVariables),sep=" + ")),",data=meta_table)",sep="")))

#You can use the following to use all the environmental variables
sol_f<-cca(abund_table ~ ., data=meta_table)

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
scrs2<-scores(sol2,display=c("sp","wa","lc","bp","cn"))

#Extract site data first
df_dates<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
colnames(df_dates)<-c("x","y","Date")

p<-ggplot()
#Draw sites
# p<-p+geom_point(data=df_dates,aes(x,y)) + scale_fill_manual(values = c("black"))

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
multiplier2 <- vegan:::ordiArrowMul(scrs2$biplot)

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

#create dfs for drawing environmental factor arrows
df_arrows<- scrs$biplot*(multiplier)*2
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

df_arrows2<- scrs2$biplot*(multiplier)*2
colnames(df_arrows2)<-c("x","y")
df_arrows2=as.data.frame(df_arrows2)

#Draw arrows for significant variables
p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="Black", alpha=1)

#Label arrows
df_arrow_names <- df_arrows
 # p<-p+geom_text(data=as.data.frame(df_arrow_names),aes(x*1.2, y, label = rownames(df_arrows)),color="Black",alpha=0.75)

#Draw arrows for non-significant variables
p<-p+geom_segment(data=df_arrows2, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")), color="Black", alpha=1)
#label arrows
df_arrow_names2 <- df_arrows2
 # p<-p+geom_text(data=as.data.frame(df_arrow_names2),aes(x*1.3, y, label = rownames(df_arrows2)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

# Either choose text or points for species
#p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species)))
p<-p+geom_point(data=df_species,aes(x,y, colour = rownames(df_species)))+scale_shape_manual("",values=2)

#Make the output pretty
p<-p+theme_bw()+theme(panel.grid = element_blank(),
                      text = element_text(size=12, family="sans"))
p<-p+coord_fixed(ratio = 1)
p<-p+xlim(-2,2.75)
p<-p+ylim(-2,2.75)
p<-p+geom_hline(yintercept = 0, color = "black")
p<-p+geom_vline(xintercept = 0, color = "black")

# Save to .eps vector file
postscript("Tara_CCA2.eps", fonts = "sans")
plot(p)
dev.off()
# 
# # Save to pdf
# pdf("TARA_CCA2.pdf")
# print(p)
# dev.off()
# 
# # Save to png raster image
# png("TARA_CCA2.png", width=6, height = 6, units="in", res=600)
# print(p)
# dev.off()

#Print eigenvalues and variance explained to csv
evals.m <- as.matrix(eigenvals(sol_f))
evals<-as.data.frame(evals.m)
evals$Var_exp_percent <- (evals[,1]/sum(evals[,1]))*100
write.csv(evals, file = "E-vals.csv", row.names = TRUE)

#Print P values, etc to csv
pvals<-as.data.frame(abund_table.adonis$aov.tab)
write.csv(pvals, file = "P-vals.csv", row.names = TRUE)
