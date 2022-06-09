# Replace original mixture components with merged mixture components.
# To obtain the merged mixture component values, please follow the information given in the Supplementary Materials:
# Copy number feature encoding -> Feature distributions and mixture modelling (p. 17)

args = commandArgs(trailingOnly = TRUE)
INPUT=args[1]
print(paste("INPUT:", INPUT))
OUTPUT=args[2]
print(paste("OUTPUT:", OUTPUT))


x=readRDS(INPUT)

# Replace components in segment size
s=x$segsize

# Example on how the mean and sd of the merged component were obtained:
# Original components:
# Mean    SD      Weight
# 310,738	92,808	0.0537
# 544,169	159,909	0.0700
# New component:
# New mean is mean from original components
# mean(c(310738, 544169)) => 427454
# New sd is sd from random sample from the original components
# sd(c(rnorm(1e6, 310738, 92808), rnorm(1e6, 544169, 159909))) => 174237
# Weight is the sum of the original weights.
# This produces the new component with:
# Mean    SD      Weight
# 427454  174237  0.1237

s[4,] = c(427454, 174237, 0.1237)
s=s[-5,]

s[10,] = c(6371179,	985163,	0.0659)
s=s[-c(11,12),]

s[11,] = c(8780945, 1086234, 0.0506)
s=s[-c(12,13),]

s[12,]=c(13435586, 1439442, 0.0414)
s=s[-c(13,14),]

s[13,]=c(18665397, 1664701, 0.0331)
s=s[-14,]

x$segsize=s

# Replace components in changepoint
s=x$changepoint

s[2,]=c(0.3428350356, 0.1572276, 0.1180555609)
s=s[-3,]

s[3,]=c(1.008185883, 0.2528096, 0.37043868)
s=s[-c(4:8),]

s[4,]=c(1.407541335, 0.3988254, 0.06601043381)
s=s[-c(5:8),]

s[5,]=c(2.20266882, 0.4656683, 0.2117234981)
s=s[-c(6:9),]

x$changepoint=s

saveRDS(x, OUTPUT)
