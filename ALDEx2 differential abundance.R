library(ALDEx2)

#import datasets
abund<-read.csv("abundance.csv", row.names = 1)
treatment<-read.csv("treatments.csv", row.names = 1)

# ALDEx2 expects:
# 'reads': integer counts with columns as samples
# 'conditions': the experimental outcome
# 'denom': the log-ratio transform reference - here "all" for CLR

conditions <- factor(treatment$Treatment,levels = c('0.2 um', '3.0 um'))

tt <- aldex(reads = abund,
 conditions = conditions,
 denom = "all",
 mc.samples=128)

 

# ALDEx2 outputs a data.frame:
# 'we.eBH': the FDR-adjusted p-value
# 'effect': the effect size
# Below, we get the names of genes
# with relatively more abundance
# in the 3.0um group

tt.bh05 <- tt[tt$we.eBH < .05,]

#store genera names and significantly different genera as vectors
three_um_up <- rownames(tt.bh05[tt.bh05$effect > 0,])
pointtwo_um_up <- rownames(tt.bh05[tt.bh05$effect < 0,])
genera<-row.names(tt)

#Create "lifestyle" vector
growthmode<-data.frame("Genus"=row.names(tt),
                       "Attached"=genera %in% three_um_up,
                       "free"=genera %in% pointtwo_um_up,
                       "lifestyle"=1
                       )
#Iterate over loop to calculate lifestyle
r<-1
for (value in growthmode$lifestyle) {
  if(growthmode[r,2]==TRUE){
    growthmode[r,4]="Particle-attached"
  }  else if (growthmode[r,3]==TRUE){
    growthmode[r,4]="Free-living"
  }else {
    growthmode[r,4]="Mixed"
  }
    
  r<-r+1
}


#create a dataframe to store relevant information
outtable<-data.frame("Genus"=row.names(tt),
                     "Effect size"=round(tt$effect, digits = 4), 
                     "Difference (between)"=round(tt$diff.btw, digits = 4), 
                     "Difference (within)"=round(tt$diff.win, digits = 4), 
                     "Expected Benjamini-Hochberg P-value"=round(tt$we.eBH, digits = 4),
                     "Lifestyle"=growthmode$lifestyle
                     )

#Output a table file
write.csv(outtable, "Lifestyle Table.csv")
