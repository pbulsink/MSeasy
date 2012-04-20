\name{SearchNIST}

\alias{SearchNIST}

\title{Identification of putative molecules using the NIST mass spectral library search tool}

\description{
SearchNIST performs identification of putative molecules using the NIST mass spectral library search tool. The input file is a MSP file. It can be obtained directly from MSeasy output files by using the conversion function \code{MSeasyToMSP}. Tcltk package is needed. Warning: this function only works on Windows plateforms ! }

\usage{SearchNIST(mspfile=NULL)}

\arguments{
  \item{mspfile}{file returned by \code{MSeasyToMSP} or any MSP file compatible with NIST. Default is NULL (or if left empty), launched a popup window to specify the location of the MSP file}}

\value{ A text file called \emph{ResultsFromNIST.txt}with the results of the NIST mass spectral library search tool is created in the folder \emph{output_SearchNIST_result_time}.}
	
\author{Yann Guitton}

\seealso{
tcltk, MSeasyToMSP
}

\examples{
#Not run
\dontrun{ 
pathexample<-system.file("doc/Output_examples",
package="MSeasy")
MSeasyToMSP(file.path(pathexample,"Output_peak21.txt"), cluster=5)

SearchNIST(mspfile=NULL) # select the MSP file created by the code below
#or
SearchNIST(mspfile=file.path(pathexample,"ForNIST.msp")) 
}
}

