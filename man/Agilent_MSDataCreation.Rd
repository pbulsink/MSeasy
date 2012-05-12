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

\format{ A folder with two different sub-folders (.D folders) , each corresponding to one GC-MS analysis. Each sub-folder contains a rteres file (chromatogram). To be used in \code{MS.DataCreation}, a complementary folder containing the netCDF files (mass spectra) should be downloaded from MSeasy web site: url{http://sites.google.com/site/rpackagemseasy/downloads/Agilent_example.zip}. 
 }

\examples{
pathAgilent<-system.file("doc/AGILENT_MSDataCreation",package="MSeasy")

\dontrun{
#command lines to download netCDF files
url1<-"http://sites.google.com/site/rpackagemseasy/downloads/Agilent_example.zip"
download.file(url=url1, destfile="AgilentCDF.zip")
unzip(zipfile="AgilentCDF.zip", exdir=".") 
print(paste("Two example folders have been downloaded in ", getwd(), sep=" "))
unlink("AgilentCDF.zip") 
}
}
\keyword{datasets}
