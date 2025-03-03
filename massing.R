# Coursework MA321-6-SP: Team Project Task - R Code to Get Started
# Version:February 2025
rm(list=ls())
# --- Setup ---
# Set your working directory to the folder where your data is stored
# Example: setwd("path/to/your/directory")
# If you're using a University lab computer, ensure you save your work on your network drive 
# or back it up using cloud storage (e.g., Apple iCloud, Google Drive) or a USB stick.
# Always keep multiple backups of your work to prevent data loss.

# --- Load Data ---
# Copy the file "gene-expression-invasive-vs-noninvasive-cancer.csv" from Moodle to your working directory

InitialData <- read.csv(file = "gene-expression-invasive-vs-noninvasive-cancer T.csv")

# --- Check the Data ---
# Use the following commands to understand the structure and dimensions of the dataset
str(InitialData)
# Output Example:
# 'data.frame': 78 obs. of 4773 variables
# $ X             : int  1 2 3 4 5 6 7 8 9 10 ...
# $ J00129        : num  -0.448 -0.48 -0.568 -0.819 ...
# $ Contig29982_RC: num  -0.296 -0.512 -0.411 -0.267 ...
# $ Contig42854   : num  -0.1 -0.031 -0.398 0.023 ...

dim(InitialData)  # Returns dataset dimensions (rows and columns)
# Example Output:
# [1] 78 4773

dimnames(InitialData)[[2]][4770:4773]  # View the names of the last columns
# Example Output:
# [1] "NM_000895" "NM_000898" "AF067420" "Class"

# Summary of the dataset:
# - 78 rows (patients)
# - 4773 columns: 4772 columns represent gene expression measurements, 
#   and column 4773 contains the "Class" variable (values: 1 or 2).

# Check the distribution of the "Class" variable
table(InitialData[, 'Class'])
# Example Output:
# Class
#   1   2 
#  34  44 



# --- Randomization Setup ---

# The script assigns a subset of variables to each team.
# In the file 'teamsubsets.csv', each team number is associated with 500 columns (variables).

# Load the file 'teamsubsets.csv', which contains the numbers of the teams 
# and their associated variable subsets.
teamsubsets <- read.csv('teamsubsets.csv')

# Specify the number of your team to identify your teams's subset of variables.
# Replace 50 by the number of your team.

your_team <- 2

my_team_subset2 <- teamsubsets[your_team,]
my_team_subset2<- as.numeric(strsplit(my_team_subset2, " ")[[1]])
my_team_subset2 <- my_team_subset2[-1]

str(my_team_subset2)
# Extract the subset of variables for the number of your team.
# The result is a vector of 500 variables associated with the number of your team.

print(my_team_subset2) # Print your team's subset of variables.

# Assume that InitialData is the preloaded dataset containing the original variables.

Class <- InitialData$Class # Extract the "Class" column, which represents the labels or targets.

# Select only the columns (variables) specified in the subset (my_subset).
X <- InitialData[, as.integer(my_team_subset2)]

# Combine the "Class" column with the selected variables to create the final dataset.

MyTeam_DataSet2 <- cbind(Class, X)

# The dataset 'MyTeam_DataSet' contains:
# - The "Class" column as the first column.
# - The 500 variables associated with your registration number.

str(MyTeam_DataSet2)

#

dim(MyTeam_DataSet2)

# The data set has 78 rows (observations/patients) and 501 variables (class variable and 500 features/genes)

dimnames(MyTeam_DataSet2)[[2]]

# For example with team number 50 you get
# > dimnames(MyTeam_DataSet50)[[2]]
# [1] "Class"          "X98260"         "Contig26706_RC" "NM_002923"      "AL137736"    
# [6] "NM_006086"      "AK000345"       "NM_014668"      "Contig50184_RC" "Contig25659_RC"
# ...
# [496] "NM_001797"      "NM_015623"      "AL117661"       "NM_006762"      "NM_001360"     
# [501] "NM_002729"
# 

MyTeam_DataSet2[1:5,1:6]

# For example with team number 50 you get
# > MyTeam_DataSet50[1:5,1:6]
#   Class X98260 Contig26706_RC NM_002923 AL137736 NM_006086
# 1     2  0.114         -0.181    -0.079   -0.099    -0.709
# 2     2  0.015         -0.131     0.143   -0.067    -0.646
# 3     2 -0.439         -0.103     0.020    0.319    -0.052
# 4     2  0.263          0.128    -0.267   -0.357     0.578
# 5     2  0.215          0.139    -0.188   -0.329    -0.493

# All analysis of your group has to use the data set MyTeam_DataSet
# defined by the 500 variables identified by your team number
# 
# Avoid plagiarism by choosing for your team variable and object names,
# and by using comments to explain each line of code in your team's words


print(my_team_subset2)
sum(is.na(MyTeam_DataSet2))

library(caret)       # For feature selection & ML
library(tidyverse)   # For data manipulation
library(randomForest)# Random Forest for feature importance
library(MASS)        # LDA for dimension reduction
library(ggplot2)

#################################################################



Gene_data = MyTeam_DataSet2[,2:501]
Classes = MyTeam_DataSet2$Class

### t-test?
p_valuesT = apply(Gene_data, 2, function(gene) t.test(gene ~Classes)$p.value)
sort(p_valuesT)

top50_ttest = order(p_valuesT)[1:50]
sub_50_ttest = MyTeam_DataSet2[,top50_ttest]



#????????????????????????????????????????????????????????????????????????????????v
### LDA? # THIS DOES NOT WORK AT THE MOMENT 
lda_model = lda(MyTeam_DataSet2$Class~., data = MyTeam_DataSet2)
lda_pred = predict(lda_model)

plot(lda_pred$x[,1], col=Classes, pch=19, 
     xlab="LDA1", ylab="LDA2", main="LDA for Cancer Classification")

prcomp = prcomp(MyTeam_DataSet2, )
plot(prcomp, main = "PCA results", xlab = "Component")
summary(prcomp)

PC1 = prcomp$x[,1]
PC1
PC2 = prcomp$x[,2]
PC2
Classes = as.factor(Classes)

ggplot(MyTeam_DataSet2, aes(x = PC1, y = PC2, color = Classes)) +
  geom_point(size = 2.5) +
  labs(title = "PCA of Gene Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))
###?????????????????????????????????????????????????????????????????????????


########

#using random forest to reduce ( recursive feature eliminaton)
# supervised feature selection

#setting class as a factor
MyTeam_DataSet2$Class = as.factor(MyTeam_DataSet2$Class)
#making a control section
control = rfeControl(functions = rfFuncs, method = "cv", number = 5)
#using Random forest elimination to determine the most important genes
set.seed(22427)
RFelim_results = rfe(MyTeam_DataSet2[],MyTeam_DataSet2$Class,sizes = seq(from = 1, to = 500, by = 10),rfeControl = control)
print(RFelim_results$optVariables)
plot(RFelim_results, type = c("g","o"))

# looking at the plot, it seems to be that around 50 variables is ideal for accuracy, 
# any less than this and you reduce that data too much. has good balance
set.seed(22427)
RFelim_results = rfe(MyTeam_DataSet2[],MyTeam_DataSet2$Class,sizes = 50,rfeControl = control)
reduced_genes = print(RFelim_results$optVariables)
print(RFelim_results$results)

Randomforest_subset = MyTeam_DataSet2[,reduced_genes]


# note, it treats class as a variable, this is important for any classification problem
###################

common_cols <- intersect(colnames(Randomforest_subset), colnames(sub_50_ttest))
