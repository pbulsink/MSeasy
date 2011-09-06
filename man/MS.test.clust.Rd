\name{MS.test.clust}

\alias{MS.test.clust}

\title{Test for the best clustering method}

\description{
This function tests the efficiency of several unsupervised clustering methods to group similar mass spectra from mass spectrometry (MS) data. 
Using a dataset where molecules are already well-identified and represented by several samples/individuals' mass-spectra, the clustering algorithms are tested for their ability to find the correct structure of the dataset (correctly assign the different mass spectra to the pre-defined number of molecules). }

\usage{MS.test.clust(data_tot, nclust)}

\arguments{
  \item{data_tot}{data matrix with the name of the molecule in the first column, the name of the sample in the second column, the retention time (or retention index) in the third column and the mass spectrum displayed in the other column.}
  \item{nclust}{number of molecules in the dataset}
}
\details{
This function tests the efficiency of several unsupervised clustering methods to group similar mass spectra from mass spectrometry data. Using a dataset where molecules are already well-identified and represented by several samples/individuals mass-spectra, the clustering algorithms are tested for their ability to correctly assign the different mass spectra to the pre-defined number of molecules. 
The clustering algorithms tested are partition around medoid (pam), hierarchical divisive clustering (Diana), hierarchical agglomerative clustering (hclust), with various combinations of distance metrics and link methods. Distance metrics include euclidean, correlation and manhattan. Link methods include single, average, complete, centroid and ward.

The results of clustering algorithms are evaluated with three quality indices that assess which clustering scheme best fits the data.
The matching coefficient computes for correct assignment of each mass spectrum to the expected molecules. 
When one cluster groups the mass spectra corresponding to the same molecule, then 1 is attributed and when one cluster contains mass spectra of different molecules, then 0 is attributed. The sum is then divided by the total expected number of molecules/clusters. The value of the matching coefficient varies from 0 to 1 and 1 indicates perfect clustering. 
Matching coefficient= Number of clusters grouping mass spectra of the same molecule divided by the total number of clusters.

The second cluster validity index is called silhouette width and described by Rosseeuw (Rousseeuw, 1987).
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
  \item{silhouette.test}{           display the silhouette index}
  \item{matching.coef}{           display the matching coefficient}
This function produces a pdf file \emph{Graph_MStestClust.pdf} displaying graphics in the folder \emph{Output_Date_time} to help identifying the best clustering method.
}
\author{Elodie Courtois, Yann Guitton, Florence Nicole}

\examples{
\dontrun{
data(Data_testclust)
MS.test.clust(Data_testclust,10)
}
}

