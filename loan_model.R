###############################################
# Business Application of Data Analytics      #
# - Professor: Santiago Gallino               #
# - Student: Juan José Zunino                 #
#                                             #
###############################################

# Setup environment
rm(list=ls())
setwd('/Users/juanjozunino/Documents/GitHub/business-applications-of-data-analytics/')

# Import libraries
library(dplyr)
library(ggplot2)

# Load data
data <- read.csv('loan_ds.csv')

data$Loan_ID <- NULL

# Descriptive analytics
summary(data$Education)

trated <- data$Education

round(prop.table(table(data$Loan_Status[data$Education == 'Graduate'])), 2)
round(prop.table(table(data$Loan_Status[data$Education == 'Not Graduate'])), 2)

table(complete.cases(data))

unique(data$Married)
# Propensity score
propensity.score <- glm(Education ~ .,
                        family = binomial,
                        x = TRUE,
                        data=data)

summary(propensity_score)

logit.