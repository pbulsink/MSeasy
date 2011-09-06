\name{MS.clust}

\alias{MS.clust}

\title{Mass spectra clustering and creation of a fingerprinting or profiling matrix}

\description{
MS.clust runs unsupervised clustering methods on  mass spectra. It can identify the optimal number of clusters using a cluster validity index (silhouette width), produces different files for facilitating the identification of molecules within a complex dataset of numerous mass spectra, and returns a fingerprinting or profiling matrix for homogeneous clusters (see details below for the definition of homogeneous clusters).}

\usage{MS.clust(data_tot, quant=FALSE, clV, ncmin, ncmax, Nbc, varRT = 0.1, 
disMeth, linkMeth, clustMeth)}
\format{
\itemize{
	\item{\code{header line  }}{  the first row must contains columns' names}
    \item{\code{first column  }}{  name of the sample/analysis}
    \item{\code{second column  }}{  retention time (RT) of the peak (or retention index (RI)) }
    \item{\emph{optionally} \code{third column  }}{  corrected peak area }
    \item{\emph{optionally} \code{fourth column  }}{  percent of the total corrected area }
    \item{\code{following columns  }}{  relative mass spectrum of the peak (mass spectrum at the apex or obtained by averaging 5 percent of the mass spectra surrounding the apex; Each mass fragment is transformed to a relative percentage of the highest mass fragment per spectrum); The intensity of one mass fragment (m/z) per column}
 }}

\arguments{
  \item{data_tot}{R object data frame as returned by \emph{MS.DataCreation} (\code{initial_DATA.txt}), or a \emph{user made file} (.txt, .csv...) with the first row containing columns' names; first column contains sample/analysis name; second column contains retention time of the peak (or retention index); optionally third and fourth columns may contain respectively corrected  peak area and percent of the total corrected area, and following columns contains the mean relative mass spectrum (the  intensity of one mass fragment (m/z) per column; each mass fragment is transformed to a relative percentage of the highest mass fragment per spectrum)}

  \item{clV}{TRUE indicates that the function clValid will be used to identify the optimal number of clusters. FALSE, when the number of clusters is already known, escapes the clValid step and goes directly to the clustering. }

  \item{ncmin}{If clV = TRUE, a numeric value giving the minimum number of clusters to be evaluated.}

  \item{ncmax}{If clV = TRUE, a numeric value giving the maximum number of clusters to be evaluated.}

  \item{Nbc}{If clV = FALSE, a numeric vector giving the number(s) of clusters to be evaluated. e.g., Nbc=c(20,25) would evaluate the number of clusters 20 and 25. }

  \item{varRT}{RT or RI's range to define homogeneous clusters, i.e. the accepted range of variation of RT/RI for a given molecule. Default value is set to 0.1.
				If the varRT is baseless (analyses from different GC columns for example), set the varRT to a high value (1000 for example).}

  \item{clustMeth}{A character vector giving the clustering methods. Available options are \emph{hierarchical}, \emph{diana}, \emph{kmeans} and \emph{pam}}
  
  \item{disMeth}{The metric used to determine the distance matrix. Possible choices are \emph{euclidean}, \emph{manhattan}, and \emph{correlation}. For pam and diana, only euclidean and manhattan are available.}

  \item{linkMeth}{For hierarchical clustering, the agglomeration method used. Available choices are \emph{ward}, \emph{single}, \emph{complete}, \emph{centroid} and \emph{average}. For all others \code{clustMeth}, \code{linkMeth=NULL}}

  \item{quant}{\code{TRUE} if columns 3 and 4 of the input file contains corrected peak area (CorrArea) and percent of the total corrected area (PercTot). CorrArea is used for absolute quantification when associated with the use of external and/or internal standards. PercTot is used for relative quantification (no external or internal standard needed). This option generates two distinct profiling matrices in outfiles, one with CorrArea and one with PercTot. \code{FALSE} if these two columns are absent. Then, a fingerprinting matrix (absence or presence of each molecule) is generated}
  }
\details{
MS.clust runs several unsupervised clustering methods on a dataset composed of numerous mass spectra from different samples/analyses. 
When the total number of molecules in the dataset is unknown, MS.clust can first identify the optimal number of clusters with a cluster validity index (silhouette width) after running the clustering on a range of numbers of clusters (clValid procedure, \code{clV=TRUE}).  

A graphic window displays the mean silhouette width as a function of the number of clusters. A red line indicates the optimal number of clusters with the highest silhouette width. The values of silhouette width for the different numbers of clusters are summarized in a table named \emph{res_clValid.txt} and saved in a folder called \emph{Output_resultdate_time}.  

Following the graphic display, the user is asked \emph{How many clustering separations ?}, i.e. how many times  should the dataset be cut into clusters. The answer is an integer. If the graph indicates a unique and clear optimal at the apex of a curve, only one cut at the optimal number of clusters is expected. If the graph display an optimal value located on a plateau, the user might be interested to perform different cut, one at the optimal number of clusters together with one at the minimum and one at the maximum numbers of clusters delimiting the plateau.  

Afterward, the user is asked \emph{How many clusters?} The answer is an integer. If several values, each integer should be entered and followed by Enter key.
When the number(s) of clusters to be analyzed is defined, the clustering is performed. The arguments of clustMeth includes hierarchical, diana, kmeans and pam. For disMeth and linkMeth, the arguments are similar to those of the clValid package. See arguments below and the documentation of this package for more details.
Following the clustering, the function identifies homogeneous and inhomogeneous clusters. The homogeneous clusters are defined by a variation in retention time lower than varRT (0.1 by default). Homogeneous clusters may correspond to a well-defined molecule, with clear mass spectra. Inhomogeneous clusters usually need manual investigations to be further classified as molecules.
}

\value{
MS.clust produces different files in folder \emph{Output_MSclust_resultdate_time} for facilitating the identification of molecules within a dataset composed of numerous mass spectra:
\item{Output_cluster.txt}{in column, the number of the cluster, quality of the cluster based on the variation of retention time (0 if inhomogeneous, 1 if homogeneous), number of distinct individuals within the cluster and total number of peaks in the cluster (to check for unique occurence of each given analysis in the cluster), mean retention time (RT), range of retention time (max(RT)-min(RT)), mean silhouette width. Follow the 8 highest mass fragments (m/z) and the complete mean relative mass spectrum.}
\item{Output_peaks.txt}{in column, the number of the cluster, the sample name, the retention time, the silhouette width, the neighbor cluster, \emph{optionally} if \code{quant=TRUE} corrArea and PercTotal, the 8 highest mass fragments and the complete mean relative mass spectrum.}
\item{Hist_cluster_ok_RT.pdf}{a pdf file displaying the histogram of the distribution of retention times for each homogeneous cluster.}
\item{Hist_cluster_ok_silhouette.pdf}{a pdf file displaying the histrogram of the distribution of silhouette width for each homogeneous cluster.}
\item{Hist_cluster_problem_RT.pdf}{a pdf file displaying the histrogram of the distribution of retention times for each inhomogeneous cluster.}
\item{Hist_cluster_problem_silhouette.pdf}{a pdf file displaying the histrogram of the distribution of silhouette width for each inhomogeneous cluster.}

Depending on the quant option

\item{Output_fingerprintingmatrix.txt}{a fingerprinting matrix (0 for absence, 1 for presence) with samples' names in the first column, retention time in the second column and presence or absence for homogeneous clusters in the following columns.}
 
or 

\item{Output_profilingmatrix_CorrArea.txt}{a profiling matrix (0 for absence, corrected area if present) with samples' names in the first column, retention time in the second column and corrected area for homogeneous clusters in the following columns.}

\item{Output_profilingmatrix_PercTot.txt}{a profiling matrix (0 for absence, percent of the total corrected area if present) with samples' names in the first column, retention time in the second column and percent of the total corrected area for homogeneous clusters in the following columns.}
}

\author{Elodie Courtois, Yann Guitton, Florence Nicole}

\seealso{
cluster, kohonen, class, mclust, amap, ClValid, fpc, flexmix
}

\examples{

data(Agilent_quantT_MSclust)
MS.clust(Agilent_quantT_MSclust, quant=TRUE, clV=TRUE, ncmin=10, ncmax=50,
 varRT = 0.1, disMeth="euclidean", linkMeth="ward", clustMeth="hierarchical")
1
21
## 21 clusters have been determined as the optimal number of cluters.
##with the option quant=TRUE, generate profiling matrices in output

data(Agilent_quantF_MSclust)
MS.clust(Agilent_quantF_MSclust, quant=FALSE, clV=FALSE, Nbc=21, 
varRT = 0.1, disMeth="euclidean", linkMeth="ward", clustMeth="hierarchical") 
##with clV=FALSE, if you already know the number of molecules in the dataset
##with the option quant=FALSE, generate a fingerprinting matrix in output

data(ASCII_MSclust)
MS.clust(ASCII_MSclust, quant=FALSE, clV=TRUE, ncmin=10, ncmax=50, 
varRT = 0.1, disMeth="euclidean", linkMeth=NULL, clustMeth="kmeans") 
3 	
26
28
30
## output files are generated for three different numbers of clusters.
## with 3 as the number of clustering separations
## 26 # First number of clusters 
## 28 # Second number of clusters
## 30 # Third number of clusters

}

