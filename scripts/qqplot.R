


# Applied Multivariate Data Analysis
# Final Project
# Emrecan Ozdogan and Kyle Naddeo

#---------------#
# Chi-2 QQ plot #
#---------------#

# Read in data
data  = read.table("../data/plasma_data.txt", header=F)
names = read.table("../data/plasma_variable_names.txt", header=F, stringsAsFactors=F)
colnames(data) = names[1,]

# Remove Noncontinuous variables
data = data[, -which(names(data) %in% c("SEX", "SMOKSTAT", "VITUSE", "BETAPLASMA", "RETPLASMA"))]

# Dimensions
feature_size = ncol(data)
sample_size  = nrow(data)

# Parameters
s = cov(data)
xbar = colMeans(data)

# Manahalanobid Distance
data.cen = scale(data, center=T, scale=T)
d2 = diag(data.cen%*%solve(s)%*%t(data.cen))

# Chi2 Quantiles
qchi = qchisq((1:sample_size - 0.5)/sample_size, df=feature_size)

# Sorted d2 value
sorted_d2 = sort(d2, index.return=TRUE)

# Plot
plot(qchi, sorted_d2$x, pch=19, xlab="Chi-2 Quantiles", 
     ylab="Mahalanobis squared distances", main="Chi-2 QQ Plot")

# Mark the outliers
num_outliers = 4
points(qchi[(sample_size-num_outliers+1):sample_size], 
       sorted_d2$x[(sample_size-num_outliers+1):sample_size], cex = 3,col='blue')


print("Outlier Indicies")
outliers = sorted_d2$ix[(sample_size-num_outliers+1):sample_size]
reduced_data = data[-outliers, ]


# Pairs
pairs(reduced_data, pch=19, col = "#CC6600", lower.panel=NULL, font.labels = 3)

# Correlation Plot
corrplot(cor(reduced_data), method="color", type = "upper", addCoef.col="black", outline=F, diag=F, 
         col=colorRampPalette(c("deepskyblue1","white","indianred3"))(200), 
         tl.cex = 1, number.cex = 1,  cl.cex = 1,  tl.col = "black")

