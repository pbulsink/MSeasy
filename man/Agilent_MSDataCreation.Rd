\name{Agilent_MSDataCreation}

\alias{Agilent_MSDataCreation}

\docType{data}
\title{
Demonstration folder for MS.DataCreation with option DataType=Agilent
}

\description{
This demonstration folder includes 2 GC-MS analyses of Lavandula obtained from Agilent. The two analyses represent a total of 54 chromatogram's peaks. The folder can be used with the function MS.DataCreation that collects and assembles the information from chromatograms and mass spectra of the two samples in a initial data matrix with peaks in row and mass spectrum in columns.  
}

\usage{Agilent_MSDataCreation}

\format{ A folder with two different sub-folders, each corresponding to one GC-MS analysis. Each sub-folder contains an netCDF file (mass spectra) and a rteres file (chromatogram). 
 }

\examples{
data(Agilent_MSDataCreation)
}

\keyword{datasets}
