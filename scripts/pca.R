# Applied Multivariate Data Analysis
# Final Project
# Emrecan Ozdogan and Kyle Naddeo

#------------------------------#
# Principal Component Analysis #
#------------------------------#

# Imports
library(corrplot)

# Read in data
data  = read.table("../data/plasma_data.txt", header=F)
names = read.table("../data/plasma_variable_names.txt", header=F, stringsAsFactors=F)
colnames(data) = names[1,]

# Remove Noncontinuous variables
data = data[, -which(names(data) %in% c("SEX", "SMOKESTAT", "VITUSE"))]

# Standardize
std.data = scale(data, center=T, scale=T)

# Preform PCA
data.pca = prcomp(std.data)
summary(data.pca)

## Visualizations ##
# Pairs
pairs(std.data, pch=19, col = "#CC6600", lower.panel=NULL, font.labels = 3)

# Correlation Plot
corrplot(cor(data), method="color", type = "upper", addCoef.col="black", outline=F, diag=F, 
         col=colorRampPalette(c("deepskyblue1","white","indianred3"))(200), 
         tl.cex = 1, number.cex = 1,  cl.cex = 1,  tl.col = "black")

# Scree Plot
screeplot(data.pca, type="line", main="Scree Plot") 

# Scatter of first two
plot(data.pca$x[, 1:2], pch=19)
text(data.pca$x[, 1:2],labels=seq(length(data.pca$x[,1])), pos=3)



