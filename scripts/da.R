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
data  = read.table("../data/plasma_data.txt", header=F)
names = read.table("../data/plasma_variable_names.txt", header=F, stringsAsFactors=F)
colnames(data) = names[1,]

# Remove Noncontinuous variables
data = data[, -which(names(data) %in% c("SEX", "SMOKESTAT", "VITUSE"))]

# Add bins
num_BETAPLASMA_bins = 10
num_RETPLASMA_bins = 10
data %<>% mutate("Binned_BETAPLASMA" = cut_number(BETAPLASMA, num_BETAPLASMA_bins),
                 "Binned_RETPLASMA" = cut_number(RETPLASMA, num_RETPLASMA_bins))

# Make individual datasets for BETA and RET
BETA_data = dplyr::select(data, -c(Binned_RETPLASMA, BETAPLASMA, RETPLASMA))%>%
                   rename(Bins = Binned_BETAPLASMA)
RET_data  = dplyr::select(data, -c(Binned_BETAPLASMA, BETAPLASMA, RETPLASMA))%>%
                   rename(Bins = Binned_RETPLASMA)

# Make a reference table foe bins
BETA_category_data = cbind("Bins" = 1:num_BETAPLASMA_bins,
                            "Ranges" = levels(unique(data$Binned_BETAPLASMA)))
RET_category_data = cbind("Bins" = 1:num_RETPLASMA_bins,
                           "Ranges" = levels(unique(data$Binned_RETPLASMA)))

# Check Freq Table
table(BETA_data$Bins)
table(RET_data$Bins)

# LDA/QDA Analysis Function
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

# Preform Analysis
discrim_analysis(BETA_data)
discrim_analysis(RET_data)

