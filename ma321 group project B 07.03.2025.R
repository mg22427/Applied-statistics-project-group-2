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
library(randomForest)# Random Forest for feature importance
library(MASS)        # LDA for dimension reduction
library(ggplot2)     # for graphing
library(tidyverse)   # For data manipulation, is a must have just in case 
library(factoextra)  # extra plotting of PCA
library(viridis)     # better visualisation
#################################################################



Gene_data = MyTeam_DataSet2[,2:501]
Classes = MyTeam_DataSet2$Class

### t-test?
p_valuesT = apply(Gene_data, 2, function(gene) t.test(gene ~Classes)$p.value)
sort(p_valuesT)

top50_ttest = order(p_valuesT)[1:50]
sub_50_ttest = MyTeam_DataSet2[,top50_ttest]

# i would love to have tried chi-square test, however i dont know anything about this dataset
# and so dont think it would be applicable as i dont know what would be "expected"


#????????????????????????????????????????????????????????????????????????????????v
### LDA? # THIS DOES NOT WORK AT THE MOMENT, results are colinear  
lda_model = lda(MyTeam_DataSet2$Class~., data = MyTeam_DataSet2)
lda_pred = predict(lda_model)

plot(lda_pred$x[,1], col=Class, pch=19, 
     xlab="LDA1", ylab="LDA2", main="LDA for Cancer Classification")

# trying to use PCA instead to reduce dimensions.
prcomp = prcomp(MyTeam_DataSet2, center = T, scale. = F) # Making a pca
plot(prcomp, main = "PCA results", xlab = "Component")   # plotting the pca

# i wasnt happy with the plotting so found a better package 
### https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/eigenvalue

fviz_eig(prcomp, choice = "eigenvalue") # plots eigenvalues of PCA
fviz_eig(prcomp, choice = "variance") # Plots explained vaience of PCA
fviz_eig(prcomp, addlabels = T)       # adding labels 
fviz_eig(prcomp, addlabels = F, geom = "bar", barcolor = "black", barfill = heat.colors(n=78), ncp = 78) # plotting all PCS
fviz_eig(prcomp, addlabels = F, geom = "bar", barcolor = "black", barfill = heat.colors(25,rev = F), ncp = 25) # thinning to 25 
fviz_eig(prcomp, addlabels = T, geom = , barcolor = "black", barfill = heat.colors(10,rev = F), ncp = 10) # thinnig to 10
# looking at both the graph and the summary data, we can see that after the first 9 or 10 components, there is a diminishing return of the
# explained variance of each PC. It takes another 68 PC's to make up another 50 percent of the variance explained
# 

PC1 = prcomp$x[,1] # making PC1 a vector
PC2 = prcomp$x[,2] # making PC2 a vector

# Plotting PC1 and PC2 against eachother with CLASS
ggplot(MyTeam_DataSet2, aes(x = PC1, y = PC2, color = as.factor(Class))) +
  geom_point(size = 2.5) +
  labs(title = "PCA of Gene Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))

# this seems to not work very well, theres lots of clustering with all of the components in PC1 and PC2


###?????????????????????????????????????????????????????????????????????????


########
# * * * 1
#using random forest to reduce ( recursive feature eliminaton)
# supervised feature selection

#setting class as a factor
MyTeam_DataSet2$Class = as.factor(MyTeam_DataSet2$Class)
#making a control section
control = rfeControl(functions = rfFuncs, method = "cv", number = 5)
#using Random forest elimination to determine the most important genes
set.seed(321)
RFelim_results = rfe(MyTeam_DataSet2[],MyTeam_DataSet2$Class,sizes = seq(from = 1, to = 500, by = 10),rfeControl = control)
print(RFelim_results$optVariables)
plot(RFelim_results, type = c("g","o"))

# looking at the plot, it seems to be that around 50 variables is ideal for accuracy, 
# any less than this and you reduce that data too much. has good balance
set.seed(321)
RFelim_results = rfe(MyTeam_DataSet2[],MyTeam_DataSet2$Class,sizes = 51,rfeControl = control)
reduced_genes = print(RFelim_results$optVariables)
print(RFelim_results$results)

Randomforest_subset = MyTeam_DataSet2[,reduced_genes]



####
#references i have read so far
 # https://www.geeksforgeeks.org/feature-selection-techniques-in-machine-learning/ * * * 1

# to surmise, there the t- test, PCA and random forest selection 
# the reduced PCA set is 
  PCA25 = prcomp$x[,1:25] # PCA scores for the first 25 PC
  PCA10 = prcomp$x[,1:10] # PCA scores for the first 10 PC

# for Random forest the top 50 genes are 
  Randomforest_subset

  
  
  #33333333333
  library(factoextra)  # k-means, hierarchical clustering
  library(cluster)      # to help analyse clustering
  library(dbscan)       # DBSCAN clustering
  library(NbClust)      # to find k
  
  # prepare data (without Class)
  Gene_data <- MyTeam_DataSet2[, -1]  # Class 열 제거 (비지도 학습이므로)
  scaled_Gene_data <- scale(Gene_data)  # 스케일링 (k-means에 필수)
  
  #K-Means Clustering
  set.seed(321)
  # finding k (Elbow Method)
  fviz_nbclust(scaled_Gene_data, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    ggtitle("Optimal K using Elbow Method")
  # for the elbow method, you find the point where a sharp decline slows down so it is 4.
  
  # finding k (Silhouette)
  fviz_nbclust(scaled_Gene_data, kmeans, method = "silhouette") +
    ggtitle("Optimal K using Silhouette Method")
  # for the silhouette method, you find the point where the silhouette score is highest. so it is 2.
  
  
  sum(is.na(scaled_Gene_data))
  
  
  zero_var_cols <- apply(scaled_Gene_data, 2, var) == 0
  sum(zero_var_cols)
  scaled_Gene_data_fixed <- scale(scaled_Gene_data)
  
  
  set.seed(321) 
  gap_stat <- clusGap(scaled_Gene_data, FUN = kmeans, K.max = 10, B = 50)
  
  # Gap Statistic
  fviz_gap_stat(gap_stat)
  
  set.seed(321)
  # Gap Statistic to find k 
  fviz_nbclust(scaled_Gene_data, kmeans, nstart = 25, method = "gap_stat", nboot = 50) +
    labs(subtitle = "Gap Statistic Method for Optimal K")
  
  
  #K-Means comparing
  library(ggplot2)
  # k = 2
  set.seed(321)
  kmeans_2 <- kmeans(scaled_Gene_data, centers = 2, nstart = 25)
  fviz_cluster(kmeans_2, data = scaled_Gene_data) + ggtitle("K-Means Clustering (k=2)")
  # k = 4
  set.seed(321)
  kmeans_4 <- kmeans(scaled_Gene_data, centers = 4, nstart = 25)
  fviz_cluster(kmeans_4, data = scaled_Gene_data) + ggtitle("K-Means Clustering (k=4)")

  
  # checking k=2 
  sil_k2 <- silhouette(kmeans_2$cluster, dist(scaled_Gene_data))
  fviz_silhouette(sil_k2) + ggtitle("Silhouette Plot for k=2")
  # checking k=4 
  sil_k4 <- silhouette(kmeans_4$cluster, dist(scaled_Gene_data))
  fviz_silhouette(sil_k4) + ggtitle("Silhouette Plot for k=4")
  

  #gap stat(k=1) shows that the data is clustered to a single group so it probably saying clustering is unnecessary
  #elbow method (k=4) it suggest that four clusters may be appropriate.
  #silhouette method(k=2) shows the highest internal cohesion occurs when the data is divided to two groups
  #so I think it is appropriate to use either 2 or 4. 
  # as the data naturally divides to two distinct group like invasive and non invasive, so I think it is good to chose 2.
  # if the cluster shows too simple then lets try k=4.
  # or we can delete the outlier and do the clustering again.

  # k=2 K-Means 
  set.seed(321)
  kmeans_res <- kmeans(scaled_Gene_data, centers = 2, nstart = 25)
  
  # K-Means visulisation
  fviz_cluster(kmeans_res, data = scaled_Gene_data, ellipse.type = "norm") +
    ggtitle("K-Means Clustering")
  
  #Hierarchical Clustering (euclidean)
  dist_matrix <- dist(scaled_Gene_data, method = "euclidean")
  
  # using ward.D
  hclust_res <- hclust(dist_matrix, method = "ward.D")
  
  # Hierarchical Clustering Dendrogram 
  plot(hclust_res, labels = FALSE, main = "Hierarchical Clustering Dendrogram")
  rect.hclust(hclust_res, k = 2, border = "red")  # k=2
  
  #DBSCAN Clustering
  set.seed(321)
  # DBSCAN 
  dbscan_res <- dbscan(scaled_Gene_data, eps = 1.5, minPts = 5)
  # visualise DBSCAN
  fviz_cluster(dbscan_res, data = scaled_Gene_data, stand = FALSE) +
    ggtitle("DBSCAN Clustering")
  #(Silhouette Score)
  silhouette_kmeans <- silhouette(kmeans_res$cluster, dist_matrix)
  fviz_silhouette(silhouette_kmeans) +
    ggtitle("Silhouette Score for K-Means Clustering")
  #https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
  #https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
  
  