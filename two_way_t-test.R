#import data from csv files
attach_abund <-read.csv("3um_top2k_abund.csv", row.names = "OTU_ID", check.names=FALSE)
free_abund <-read.csv("02um_top2k_abund.csv", row.names = "OTU_ID", check.names=FALSE)
taxonomy_attach <- read.csv("3um_taxonomy.csv", row.names = "OTU_ID", check.names = FALSE)
taxonomy_free <- read.csv("02um_taxonomy.csv", row.names = "OTU_ID", check.names = FALSE)

#make a data frame to store information
Ttest_result <- data.frame(OTU_ID = character(), Lifestyle = character(), pvalue = character(), Avg_attach_abund = character(), Avg_free_abund = character(), Taxonomy = character())


# for loop checking for matching denovo sequences
for (bacName in row.names(attach_abund)){
  
  #store bacName abundances from each dataset as vectors a and b
  a = attach_abund[bacName,]
  b = free_abund[bacName,]
  tax = taxonomy_attach[bacName,]
  #store the average abundance from each dataset
  attach_avg = mean(as.numeric(attach_abund[bacName,]), na.rm = TRUE)
  free_avg = mean(as.numeric(free_abund[bacName,]), na.rm = TRUE)
  
  #do this if bacteria appears in both lists
  if (bacName %in% row.names(free_abund)){

    #Perform two-way t test on 2 data vectors. Put values into "ttest"
    ttest = t.test(a, b, alternative = "two.sided", mu=0, conf.level = 0.95)
    #store more abundant community
    if ((ttest$statistic > 0)&(ttest$p.value<=0.05)) {
      community = "Attached"
    } else if ((ttest$statistic<0)&(ttest$p.value<=0.05)) {
      community = "Free"
    } else{
      community = "Mixed"
    }
    #create a 1 row df with all info for bacName
    new_result <- list (OTU_ID = bacName, Lifestyle = community, pvalue = ttest$p.value, Avg_attach_abund = attach_avg, Avg_free_abund = free_avg, Taxonomy = tax)
  
  #do this if bacteria is only in 3.0 fraction  
  }else{
    #create a 1 row df with all info for bacName
    new_result <- list (OTU_ID = bacName, Lifestyle = "Attached Only", pvalue = "NA", Avg_attach_abund = attach_avg, Avg_free_abund = "NA", Taxonomy = tax)
  }
  #add the new info to the existing results dataframe
  Ttest_result = rbind(Ttest_result, new_result, stringsAsFactors = FALSE)
}

#do this if bacteria is only in 0.2 fraction
for (bacName in row.names(free_abund)){
  tax = taxonomy_free[bacName,]
  if (!(bacName %in% row.names(attach_abund))){
    free_avg = mean(as.numeric(free_abund[bacName,]), na.rm = TRUE)
    #create a 1 row df with all info for bacName
    new_result <- list (OTU_ID = bacName, Lifestyle = "Free Only", pvalue = "NA", Avg_attach_abund = "NA", Avg_free_abund = free_avg, Taxonomy = tax)
    #add the new info to the existing results dataframe
    Ttest_result = rbind(Ttest_result, new_result, stringsAsFactors = FALSE)
  }
}

#save df of results as .csv file
write.csv(Ttest_result, file = "Community t_test result 2k.csv")