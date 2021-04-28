# Applied Multivariate Data Analysis
# Final Project
# Emrecan Ozdogan and Kyle Naddeo

#--------------------------------------------#
# Linear and Quadratic Discriminant Analysis #
#--------------------------------------------#

# Import
library(dplyr)
library(magrittr)
library(ggplot2)
library(MASS)

# Read in data
#---------------------------------------------------------------------------------------#
data  = read.table("../data/plasma_data.txt", header=F)
names = read.table("../data/plasma_variable_names.txt", header=F, stringsAsFactors=F)
colnames(data) = names[1,]

# Remove Noncontinuous variables
data = data[, -which(names(data) %in% c("SEX", "SMOKSTAT", "VITUSE"))]

# Get Individual datasets Function
#------------------------------------------------------------------------------------------------------#
get_plasma_datasets = function(dataset, num_BETAPLASMA_bins, num_RETPLASMA_bins){
  data_w_bins = mutate(dataset, "Binned_BETAPLASMA" = cut_number(BETAPLASMA, num_BETAPLASMA_bins),
                   "Binned_RETPLASMA" = cut_number(RETPLASMA, num_RETPLASMA_bins))
  
  # Make individual datasets for BETA and RET
  BETA_data = dplyr::select(data_w_bins, -c(Binned_RETPLASMA, BETAPLASMA, RETPLASMA))%>%
                     rename(Bins = Binned_BETAPLASMA)
  RET_data  = dplyr::select(data_w_bins, -c(Binned_BETAPLASMA, BETAPLASMA, RETPLASMA))%>%
                     rename(Bins = Binned_RETPLASMA)
  
  # Make a reference table for bins
  BETA_category_data = cbind("Bins" = 1:num_BETAPLASMA_bins,
                              "Ranges" = levels(unique(data$Binned_BETAPLASMA)))
  RET_category_data = cbind("Bins" = 1:num_RETPLASMA_bins,
                             "Ranges" = levels(unique(data$Binned_RETPLASMA)))
  
 return(list("BETA_data" = BETA_data, "Beta_category_data" = BETA_category_data,
             "RET_data" = RET_data, "RET_category_data" = RET_category_data))
}

# LDA/QDA Analysis Function
#---------------------------------------------------------------------------------------------------------#
discrim_analysis = function(discrim_data){
  # LDA
  lda.discrim_data = lda(Bins ~ ., data = discrim_data)
  lda.predict = predict(lda.discrim_data)
  
  lda_prediction_table = table(lda.predict$class, discrim_data[ ,ncol(discrim_data)])
  lda_accuracy = sum(diag(lda_prediction_table))/sum(lda_prediction_table)
  
  # QDA
  qda.discrim_data = qda(Bins ~ ., data = discrim_data)
  qda.predict = predict(qda.discrim_data)
  
  qda_prediction_table = table(qda.predict$class, discrim_data[ ,ncol(discrim_data)])
  qda_accuracy = sum(diag(qda_prediction_table))/sum(qda_prediction_table)
  
  # Return Accuracies
  return(list(lda_accuracy, qda_accuracy))
}

# Analysis
# --------------------------------------------------------------------#
min_num_bins = 4
max_num_bins = 50
BETA_optimal_bin_acc = c(min_num_bins, 0, min_num_bins, 0)
RET_optimal_bin_acc = c(min_num_bins, 0, min_num_bins, 0)
for (i in seq(min_num_bins, max_num_bins)){
    # Get Individual data sets
    binned_data = get_plasma_datasets(data,
                                      num_BETAPLASMA_bins = i,
                                      num_RETPLASMA_bins = i)
    
    ## BETA ##
    # Check each bin has sufficient population
    if (min(table(binned_data$BETA_data$Bins)) > 10){
      # Preform Analysis
      BETA_results  = discrim_analysis(binned_data$BETA_data)
      
      # Save Best
      if (as.numeric(BETA_results[1]) > BETA_optimal_bin_acc[2]){
        BETA_optimal_bin_acc[1:2] = c(i, as.numeric(BETA_results[1]))
      }
      if (as.numeric(BETA_results[2]) > BETA_optimal_bin_acc[4]){
        BETA_optimal_bin_acc[3:4] = c(i, as.numeric(BETA_results[2]))
      }
      
    }
    
    ## RET ##
    # Check each bin has sufficient population
    if (min(table(binned_data$RET_data$Bins)) > 10){
      # Preform Analysis
      RET_results  = discrim_analysis(binned_data$RET_data)
      
      # Save Best
      if (as.numeric(RET_results[1]) > RET_optimal_bin_acc[2]){
        RET_optimal_bin_acc[1:2] = c(i, as.numeric(RET_results[1]))
      }
      if (as.numeric(RET_results[2]) > RET_optimal_bin_acc[4]){
        RET_optimal_bin_acc[3:4] = c(i, as.numeric(RET_results[2]))
      }
      
    }
    
}

# Display results
results = round(rbind(BETA_optimal_bin_acc,
          RET_optimal_bin_acc), digits = 2)
colnames(results) = c("LDA Bins", "LDA Acc", "QDA Bins" , "QDA Acc")
rownames(results)  = c("BETA", "RET")
print(results)

