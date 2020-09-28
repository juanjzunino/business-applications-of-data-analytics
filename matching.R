"
                          Business Application of Data Analytics      
                                
                              Professor: Santiago Gallino               
                               Student: Juan José Zunino                 
                                             
"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                     CONFIGURATION                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################
#        Environment        #
#############################

# Setup environment
rm(list=ls())
setwd('/Users/juanjozunino/Documents/GitHub/bada/')

# Import libraries
library(dplyr)
library(ggplot2)
library(coin)
library(optmatch)
library(exactRankTests)
library(MASS)
library(gap)



#############################
#           Data            #
#############################

# Load data
data <- read.csv('Data/dataset.csv', header = TRUE, na.strings = c(""))

# Load functions
source('functions.R')

# Parse data
dim(data)
summary(data) # Gender, Marrried and Self_Employed presents missing data

data$Loan_Status <- ifelse(data$Loan_Status == 'Y', 1, 0)
data$Dependents <- as.factor(ifelse(data$Dependents == 0, 0, 1))
data$Credit_History <- as.factor(data$Credit_History)
data$Loan_ID <- NULL


# Remove rows with missing data
table(complete.cases(data)) # 480 complete cases
data <- data[complete.cases(data),]
summary(data)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                 DESCRIPTIVE STATISTICS                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Proportions
table(data$Education)

round(prop.table(table(data$Loan_Status[data$Education == 'Graduate'])), 2)
round(prop.table(table(data$Loan_Status[data$Education == 'Not Graduate'])), 2)

# t-test
treated.y=data$Loan_Status[data$Education == 'Graduate'];
control.y=data$Loan_Status[data$Education == 'Not Graduate'];

t.test(treated.y, control.y)

summary(lm(Loan_Status ~ ., data=data))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                    PROPENSITY SCORE                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################
#     Propensity score      #
#############################

# Model
propscore.model <- glm(Education ~ .,
                        family = binomial,
                        x = TRUE,
                        y = TRUE,
                        data = data[,setdiff(colnames(data), "Loan_Status")])
summary(propscore.model)

# Set treatment variable
treated <- propscore.model$y



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                    MATCHING (BALANCE)                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################
#          Logit            #
#############################

# Find logit(propensity score)
logit.propscore <- predict(propscore.model)


#############################
#      Distance matrix      #
#############################

# Construct a distance matrix which gives the absolute 
# difference between the propensity scores 
# of the treated and control subjects

distmat.propensity=matrix(rep(0,sum(treated==1)*sum(treated==0)),nrow=sum(treated==1))

for(i in 1:sum(treated==1)){
  for(j in 1:sum(treated==0)){
    distmat.propensity[i,j]=abs((logit.propscore[treated==1])[i]-(logit.propscore[treated==0])[j])
  }
}



# Create a subject index and name the rows and columns of distance
# matrix by ### this subject index
subject.index=seq(1,length(treated),1)
rownames(distmat.propensity)=subject.index[treated==1]
colnames(distmat.propensity)=subject.index[treated==0]


#############################
#       Pair matching       #
#############################

# Generate matches
matchvec <- pairmatch(distmat.propensity)

# Create vectors of the subject indices of the treatment
treated.subject.index=rep(0,sum(treated==1))
matched.control.subject.index=rep(0,length(treated.subject.index))
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
subjects.match.order=as.numeric(names(matchvec))

# The subject indices in 
# the order of matchvec
for(i in 1:length(treated.subject.index)){
  matched.set.temp=which(matchedset.index.numeric==i)
  matched.set.temp.indices=subjects.match.order[matched.set.temp]
  if(treated[matched.set.temp.indices[1]]==1){
    treated.subject.index[i]=matched.set.temp.indices[1]
    matched.control.subject.index[i]=matched.set.temp.indices[2]
  }
  if(treated[matched.set.temp.indices[2]]==1) {
    treated.subject.index[i]=matched.set.temp.indices[2]
    matched.control.subject.index[i]=matched.set.temp.indices[1]
  }
}


#############################
# Standardized Differences  #
#############################

# Calculate standardized differences 
# Covariates used in propensity score model

Xmat=propscore.model$x[,-1] # Get rid of the intercept
treatedmat=Xmat[treated==1,];

# Standardized differences before matching
controlmat.before=Xmat[treated==0,];
controlmean.before=apply(controlmat.before,2,mean);
treatmean=apply(treatedmat,2,mean);
treatvar=apply(treatedmat,2,var);
controlvar=apply(controlmat.before,2,var);
stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);

# Standardized differences after matching
controlmat.after=Xmat[matched.control.subject.index,];
controlmean.after=apply(controlmat.after,2,mean);

# Standardized differences after matching
stand.diff.after=(treatmean-controlmean.after)/sqrt((treatvar+controlvar)/2);
cbind(stand.diff.before,stand.diff.after)


#############################
#           t-test          #
#############################

# Before
pval.t.test.b=rep(0,ncol(Xmat)-1);
varnames=(dimnames(Xmat)[[2]])[2:ncol(Xmat)];
for(j in 2:ncol(Xmat)){
  pval.t.test.b[j-1]=t.test(treatedmat[,j],controlmat.before[,j])$p.val;
}

cbind(varnames,pval.t.test.b);

# After
pval.t.test=rep(0,ncol(Xmat)-1);
varnames=(dimnames(Xmat)[[2]])[2:ncol(Xmat)];
for(j in 2:ncol(Xmat)){
  pval.t.test[j-1]=t.test(treatedmat[,j],controlmat.after[,j])$p.val;
}

cbind(varnames,pval.t.test);


#############################
#       Conclusions         #
#############################

treated.y=data$Loan_Status[treated.subject.index];
control.y=data$Loan_Status[matched.control.subject.index];

t.test(treated.y, control.y)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                   MATCHING (CALIPER)                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#############################
#   Matrix of covariates    #
#############################

# Matrix of covariates, excluding intercept
Xmat=propscore.model$x[,-1]

# Rank based Mahalanobis distance
distmat=smahal(treated,Xmat)


#############################
#          CALIPER          #
#############################

# Add caliper
logit.propscore=predict(propscore.model)
distmat2=addcaliper(distmat,treated,logit.propscore)

### Create a subject index and name the rows and columns
#of distance matrix by ### this subject index

subject.index=seq(1,length(treated),1)
rownames(distmat2)=subject.index[treated==1]
colnames(distmat2)=subject.index[treated==0]

#############################
#       Pair matching       #
#############################


# Generate matches
matchvec=pairmatch(distmat2)


# Create vectors of the subject indices of the treatment
treated.subject.index=rep(0,sum(treated==1))
matched.control.subject.index=rep(0,length(treated.subject.index))
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
subjects.match.order=as.numeric(names(matchvec))


# The subject indices in 
# the order of matchvec
for(i in 1:length(treated.subject.index)){
  matched.set.temp=which(matchedset.index.numeric==i)
  matched.set.temp.indices=subjects.match.order[matched.set.temp]
  if(treated[matched.set.temp.indices[1]]==1){
    treated.subject.index[i]=matched.set.temp.indices[1]
    matched.control.subject.index[i]=matched.set.temp.indices[2]
  }
  if(treated[matched.set.temp.indices[2]]==1){
    treated.subject.index[i]=matched.set.temp.indices[2]
    matched.control.subject.index[i]=matched.set.temp.indices[1]
  }
}



#############################
# Standardized Differences  #
#############################

# Calculate standardized differences 
# Covariates used in propensity score model
Xmat=propscore.model$x[,-1]
treatedmat=Xmat[treated==1,];

# Standardized differences before matching
controlmat.before=Xmat[treated==0,];
controlmean.before=apply(controlmat.before,2,mean);
treatmean=apply(treatedmat,2,mean);
treatvar=apply(treatedmat,2,var);
controlvar=apply(controlmat.before,2,var);
stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);

# Standardized differences after matching
controlmat.after=Xmat[matched.control.subject.index,];
controlmean.after=apply(controlmat.after,2,mean);

# Standardized differences after matching
stand.diff.after=(treatmean-controlmean.after)/sqrt((treatvar+controlvar)/2);

cbind(stand.diff.before,stand.diff.after)



#############################
#           t-test          #
#############################

covnames=names(stand.diff.before)[-1]
t.test.pval.vec=rep(0,length(covnames))
for(i in 2:ncol(treatedmat)){
  t.test.pval.vec[i-1]=t.test(treatedmat[,i],controlmat.after[,i])$p.value
}
cbind(covnames,t.test.pval.vec)



#############################
#       Conclusions         #
#############################

treated.y=data$Loan_Status[treated.subject.index];
control.y=data$Loan_Status[matched.control.subject.index];

t.test(treated.y, control.y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#