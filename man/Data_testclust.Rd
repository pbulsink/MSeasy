\name{Data_testclust}

\alias{Data_testclust}

\docType{data}

\title{
Demonstration dataset for MS.test.clust }

\description{
To test for the best unsupervised clustering method, a dataset where molecules are already identified is created.  Each molecule is represented by several samples' mass spectra. Here, the dataset contains 10 molecules obtained in different samples (84 Lavandula GC-MS analyses). In the function MS.test.clust, different clustering methods are tested for their abilities to find the correct structure of the dataset. Three different cluster validity indices are calculated to evaluate the results: the matching coefficient, the silhouette width and the Dunn index (see MS.test.clust for details)}

\usage{data(Data_testclust)}

\format{A data frame with 10 molecules from 84 GC-MS analyses.
  \itemize{
	\item{\code{header line  }}{  the first row must contains the columns' names}
    \item{\code{first column }}{ name of the molecule}  
    \item{\code{second column }}{ sample name}  
    \item{\code{third column }}{ retention time}  
    \item{\code{following columns }}{ mean relative mass spectrum of the molecule (the intensity of one mass fragment (m/z) per column; Mean mass spectrum calculated by averaging 5 percent of the mass spectra surrounding the apex; The intensity of each mass fragment is transformed to a relative percentage of the highest mass fragment per spectrum)}
}
}

\examples{
data(Data_testclust)
}

\keyword{datasets}


