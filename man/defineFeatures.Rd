% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runPrunedXcmsSet.R
\name{defineFeatures}
\alias{defineFeatures}
\title{define features for individual pipeline}
\usage{
defineFeatures(
  files,
  mtab_files,
  rtParam = NULL,
  pgParam = PeakDensityParam(rep(1, length(files)))
)
}
\arguments{
\item{files}{the file paths of XCMSnExp}

\item{mtab_files}{the file path of meta tables, often returned by function \code{peakPicking}}

\item{rtParam}{the retention time adjustment parameters, passed to \code{adjustRtime}}

\item{pgParam}{the peak grouping parameter, passed to \code{groupChromPeaks}}
}
\value{
It returns an object of class \code{XCMSnExp}
}
\description{
define features for individual pipeline
}
