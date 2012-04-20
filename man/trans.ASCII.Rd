\name{trans.ASCII}

\alias{trans.ASCII}

\title{Transform GC-MS data in ASCII format to suitable data for MS.DataCreation}

\description{
This function transform each ASCII file (i.e. each GC-MS analysis in ASCII format) into a new file compatible with MS.DataCreation.
}

\usage{trans.ASCII(path, mz)}

\arguments{
  \item{path}{Name of the folder containing all the GC-MS analyses in ASCII format. If left empty, a popup window opens to browse your computer \emph{tcltk package required}}
  \item{mz}{Range of mass fragments delimiting the mass spectrum (each mass fragment is characterized by its mass-to-charge ratio m/z)}
}

\details{
When your raw GC-MS data cannot be exported to netCDF (or mzXML) but only to  the international ASCII format (.txt). The data in ASCII format have to be transformed with the function trans.ASCII for further analyses with MS.DataCreation (option \emph{DataType="ASCII"}). }

\value{
trans.ASCII creates a folder named \emph{output_transASCII_Date_Hour} which contains the same number of files than \code{path}.
}

\author{Elodie Courtois, Yann Guitton, Florence Nicole}

\examples{
\dontrun{
##not run
##For ASCII GC-MS files 
path<-system.file("doc/ASCII_TransASCII",package="MSeasy")
trans.ASCII(path=path,mz=30:250)
}
}
