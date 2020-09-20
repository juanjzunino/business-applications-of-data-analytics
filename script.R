###############################################
# Business Application of Data Analytics      #
# - Professor: Santiago Gallino               #
# - Student: Juan José Zunino                 #
#                                             #
###############################################

# 1) Setup environment
rm(list=ls())
setwd('/Users/juanjozunino/Documents/GitHub/business-applications-of-data-analytics/')


# 2) Import libraries
library(dplyr)
library(ggplot2)


# 3) Load data
data <- read.csv('loan_ds.csv')


# 4) Descriptive analytics
prop.table(table(data$Loan_Status[data$Gender == 'Male']))
prop.table(table(data$Loan_Status[data$Gender == 'Female']))


# 5) Matching
