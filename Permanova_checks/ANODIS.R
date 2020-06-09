### R code to perform an Analysis of Distance using the first two Principal Components axes 

# A. read in or estimate PCA on data. ANODIS in this example will use dataframe with 3 columns: 
# 1. treatment indicator
# 2. PC1
# 3. PC2

# if data is already assembled into a matrix with the 1st two PC columns:

data.in = read.csv("Dragonfly assemblage data.csv") 

# otherwise, create data set from PCA output.  
# This example uses the PCA function in the FactoMineR package. 

pc.in=out$ind$coord[,1:2]

# where out = PCA(orig.data.in,scale.unit=T,graph=FALSE)
# data.in would then be

data.in = cbind(Treatment = pop.ind, PC1 = pc.in[,1], PC2 = pc.in[,2])

# where pop.ind is the group indicator variable 

head(data.in)

# Treatment      PC1        PC2
# 1   Natural -2.16337 -1.2300800
# 2   Natural -3.35751  0.1242660
# 3   Natural -1.72002  0.0719634
# 4   Natural -1.48639 -0.0646673
# 5   Natural -1.11385  0.3605570
# 6   Natural -0.86583  0.4979210

# B. perform appropriate linear model analysis on each principal component axis.
#    For this example, aov performs a one-way lm analysis
#    - Using summary on the aov function gives a 2 x 5 matrix in the 1st list item

test.1=summary(aov(pc1~treatment, data=data.in))
test.2=summary(aov(pc2~treatment, data=data.in))

## result of anovas
test.1
# Df Sum Sq Mean Sq F value Pr(>F)  
# Treatment    2  23.06  11.532   4.262  0.018 *
#   Residuals   68 184.00   2.706                 
# ---
#   Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1

test.2
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Treatment    2  54.37  27.185   46.32 2.03e-13 ***
#   Residuals   68  39.91   0.587                     
# ---
#   Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1


# C. Get Sums of Squares from each test to get overall results.

ssTreat.x = test.1[[1]][1,2]
ssError.x = test.1[[1]][2,2]
dfTreat.x = test.1[[1]][1,1]
dfError.x = test.1[[1]][2,1]

ssTreat.y = test.2[[1]][1,2]
ssError.y = test.2[[1]][2,2]
dfTreat.y = test.2[[1]][1,1]
dfError.y = test.2[[1]][2,1]

# D. Add Sums of Square for overall test
sumSS.trtmnt=(ssTreat.x)/(dfTreat.x)
sumSS.error=(ssError.x)/(dfError.x)

# F-test
F.all =sumSS.trtmnt/mean.error
F.df1 = dfTreat.x+dfTreat.y
F.df2 = dfError.x+dfError.y

1-pf(F.all,F.df1,F.df2)

[1] 5.992883e-05











