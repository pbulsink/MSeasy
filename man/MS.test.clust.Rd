\name{MS.test.clust}

\alias{MS.test.clust}

\title{Test for the best clustering method}

\description{
This function tests the efficiency of several unsupervised clustering methods to group similar mass spectra from mass spectrometry (MS) data. 
Using a dataset where molecules are already well-identified and represented by several samples/individuals mass spectra, the clustering algorithms are tested for their ability to find the correct structure of the dataset (correctly assign the different mass spectra to the pre-defined molecules). }

\usage{MS.test.clust(data_tot, nclust)}

\arguments{
  \item{data_tot}{data matrix with the name of the molecule in the first column, the name of the sample in the second column, the retention time (or retention index) in the third column and the relative mass spectrum displayed in the following columns.}
  \item{nclust}{number of molecules in the dataset}
}
\details{
This function tests the efficiency of several unsupervised clustering methods to group similar mass spectra from mass spectrometry data. Using a dataset where molecules are already well-identified and represented by several samples/individuals mass-spectra, the clustering algorithms are tested for their ability to correctly assign the different mass spectra to the pre-defined molecules. 


Since the total number of true molecules is usually unknown in complex biological substance, the use of unsupervised clustering algorithms is required. These include partitional and hierarchical algorithms. Partitional algorithms, such as K-means or Partitioning Around Medoids (PAM), determine all clusters at once and do not consider any hierarchic-al/neighboring relations among clusters. For these algorithms, the number of clusters should be specified beforehand. Unlike partitional methods, hierarchical algorithms are iterative methods for clustering datasets (hierarchical clustering analysis HCA), based on the neighboring relations among clusters. Two types of algorithms exist: agglomerative or divisive. Agglomerative algorithms begin with each element as a separate cluster and merge them into successively larger clusters. Divisive algorithms begin with the whole set and proceed to divide it into successively smaller clusters. 
An important step in hierarchical methods is the choice of a distance metric and a link method. Several options are implemented for the link method (average, single, complete, centroid and ward). The p-order Minkowski distance and correlations are the most commonly used measures of dissimilarity and similarity (Jobson, 1991). Minkowski distance is typically used with p being 2 (Euclidean distance) or 1 (Manhattan distance). A study found that clustering large dimension data was more efficient using p indices of Minkowski distances smaller than 1 (Aggarwal, et al., 2000; Hinneburg, et al., 2000). For that reason, we also implemented two values of p indices below 1 (p= 1/2; p=1/3
The clustering algorithms tested are Partition Around Medoid (PAM), hierarchical divisive clustering (Diana), hierarchical agglomerative clustering (hclust), with various combinations of distance metrics and link methods.

The results of clustering algorithms are evaluated with three quality indices that assess which clustering scheme best fits the data.
The matching coefficient computes for correct assignment of each mass spectrum to the expected molecules. 
When one cluster groups the mass spectra corresponding to the same molecule, then 1 is attributed and when one cluster contains mass spectra of different molecules, then 0 is attributed. The sum is then divided by the total expected number of molecules/clusters. The value of the matching coefficient varies from 0 to 1 and 1 indicates perfect clustering. 
Matching coefficient= Number of clusters grouping mass spectra of the same molecule divided by the total number of clusters.

The second cluster validity index is called silhouette width and was described by Rosseeuw (1987).
This index is based on two criteria: cluster compactness and isolation.  

Silhouette width s(i) is defined as:
s(i)=(b-a)/max(a,b)  

where a is the average distance of a point from the other points of the same cluster (variation intracluster / compactness)
and b represents the minimum of the average distances of the point from the points of the other clusters (cluster separation)


Another quality index, the Dunn index D, is defined as:

D=[min{k,l-numbers of clusters}dist(Ck, Cl)]/[max{m-cluster number}diam(Cm)]  

k,l,m - numbers of clusters which come from the same partitioning,
dist(Ck,Cl) - inter cluster distance between clusters Ck and Cl,
diam(Cm) - intra cluster diameter computed for cluster Cm. 
}

\value{
This function will return three matrices with the distance metric in column and the clustering algorithms in row.
  \item{Dunn.test}{           display the Dunn index }
  \item{silhouette.test}{           display the Silhouette Width}
  \item{matching.coef}{           display the matching coefficient}
This function produces a pdf file \emph{Graph_MStestClust.pdf} displaying graphics with matching coefficient and silhouette width in the folder \emph{output_MStestclust_resultDate_time} to help identifying the best clustering method.
}
\author{Elodie Courtois, Yann Guitton, Florence Nicole}

\examples{
#Not Run
\dontrun{
data(Data_testclust)
MS.test.clust(Data_testclust,10)
}
}

