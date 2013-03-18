\name{MSeasy-package}

\alias{MSeasy-package}

\alias{MSeasy}

\docType{package}

\title{Unsupervised and untargeted processing of Gas Chromatography-Mass Spectrometry (GC-MS) data }

\description{Package for the detection of molecules in complex mixtures of compounds.
It creates an initial_DATA matrix from several GC-MS analyses by collecting and assembling the information from chromatograms and mass spectra (\emph{MS.DataCreation}).
It tests for the best unsupervised clustering method to group similar mass spectra into molecules (\emph{MS.test.clust}).
It runs the optimal unsupervised clustering method on the initial_DATA matrix, identifies the optimal number of clusters, produces different files for facilitating the quality control and identification of putative molecules, and returns fingerprinting or profiling matrices (\emph{MS.clust}).
It converts output files from \emph{MS.clust} for NIST mass spectral library search and ARISTO webtool search}

\details{
\tabular{ll}{
Package: \tab MSeasy\cr
Type: \tab Package\cr
Version: \tab 5.3.2 \cr
Date: \tab 2012-08-17\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

}
\author{Elodie Courtois, Yann Guitton, Florence Nicole
Maintainer: Yann Guitton <yann.guitton@gmail.com>
} 


