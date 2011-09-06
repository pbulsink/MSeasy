\name{ASCII_MSclust}

\alias{ASCII_MSclust}

\docType{data}
\title{
Demonstration dataset for MS.clust 
}
\description{
This demonstration dataset includes 2 GC-MS analyses of Petrel, representing a total of 67 chromatogram's peaks. It can be used with the function MS.clust:  

(i) to identify the optimal number of clusters.

(ii) to obtain the fingerprinting matrix (absence or presence of peaks for all samples)
}

\usage{data(ASCII_MSclust)}

\format{ A data frame with 67 chromatogram's peaks from 2 GC-MS analyses.
 \itemize{
	\item{\code{header line  }}{  the first row contains columns' names}
    \item{\code{first column  }}{  name of the sample/analysis}
    \item{\code{second column  }}{  retention time of the peak }
        \item{\code{following columns  }}{  mean relative mass spectrum of the peak (the intensity of one mass fragment (m/z) per column; Mean mass spectrum calculated by averaging 5 percent of the mass spectra surrounding the apex; The intensity of each mass fragment is transformed to a relative percentage of the highest mass fragment per spectrum)}
 }}

\examples{
data(ASCII_MSclust)
}

\keyword{datasets}
