# Applied Multivariate Data Analysis
# Final Project
# Emrecan Ozdogan and Kyle Naddeo

#-----------------#
# Factor Analysis #
#-----------------#

# Read in data
data  = read.table("../data/plasma_data.txt", header=F)
names = read.table("../data/plasma_variable_names.txt", header=F, stringsAsFactors=F)
colnames(data) = names[1,]

# Remove Noncontinuous variables
data = data[, -which(names(data) %in% c("SEX", "SMOKESTAT", "VITUSE"))]

# Load Function Given to Us
eval.model = function(mod.fit) {
  resid.FA = round(mod.fit$correlation - (mod.fit$loadings[,] %*%
                                           t(mod.fit$loadings[,]) + diag(mod.fit$uniqueness)), 4)
  larger0.1  = sum(abs(resid.FA)>0.1)
  larger0.2  = sum(abs(resid.FA)>0.2)
  max.resid  = max(abs(resid.FA))
  mean.resid = colMeans(abs(resid.FA))
  list(LRT.pvalue = mod.fit$PVAL, resid.FA = resid.FA, larger0.1 = larger0.1,
       larger0.2 = larger0.2, max.resid = max.resid, mean.resid = mean.resid)
}

# Perform Factor Analysis 
num_factors = 1
mod.fit = factanal(x = data, factors = num_factors, rotation = "varimax")

while (mod.fit$PVAL < 0.05) {
  print(c(num_factors, mod.fit$PVAL))
  num_factors = num_factors + 1
  mod.fit = factanal(x = data, factors = num_factors, rotation = "varimax", scores = "regression")
}

eval.model(mod.fit)
print(c(num_factors, mod.fit$PVAL))
print(mod.fit)

# Examine Plots
FA3.positive<-mod.fit$scores[,3] - min(mod.fit$scores[,3])
common.limits<-c(min(mod.fit$scores[,1:2]), max(mod.fit$scores[,1:2]))
col.symbol<-ifelse(test = mod.fit$scores[,3]>0, yes = "red", no = "blue")
symbols(x = mod.fit$scores[,1], y = mod.fit$scores[,2], circles =
          FA3.positive, xlab = "Common factor #1", ylab = "Common factor #2", main
        =
          "Common factor scores", inches = 0.25, xlim = common.limits, ylim =
          common.limits, panel.first = grid(), fg = col.symbol)
text(x = mod.fit$scores[,1], y = mod.fit$scores[,2])
abline(h = 0)
abline(v = 0)

















