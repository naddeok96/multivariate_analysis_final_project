# Applied Multivariate Data Analysis
# Final Project
# Emrecan Ozdogan and Kyle Naddeo

#------------------------------#
# Linear Discriminant Analysis #
#------------------------------#

# Read in data
data  = read.table("../data/plasma_data.txt", header=F)
names = read.table("../data/plasma_variable_names.txt", header=F, stringsAsFactors=F)
colnames(data) = names[1,]

# Remove Noncontinuous variables
data = data[, -which(names(data) %in% c("SEX", "SMOKESTAT", "VITUSE"))]