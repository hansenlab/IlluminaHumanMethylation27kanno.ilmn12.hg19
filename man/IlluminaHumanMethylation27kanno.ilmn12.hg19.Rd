\name{IlluminaHumanMethylation27kanno.ilmn12.hg19}
\alias{IlluminaHumanMethylation27kanno.ilmn12.hg19}
\alias{Locations}
\alias{Other}
\alias{Manifest}
\docType{data}
\title{
  Annotation data for the \sQuote{IlluminaHumanMethylation27k} micro array.
}
\description{
 This package is based on the file
  \sQuote{HumanMethylation27_v-1-2.csv} from the zip archive
  from \url{www.illumina.com}.

  Additional SNP annotation is generated by the authors (described below).

  A script for creating the data object in this package is contained in
  \code{scripts/manifest.R}.
}
\format{
  An object of class \code{IlluminaMethylationAnnotation} for
  \code{IlluminaHumanMethylation27kanno.ilmn12.hg19}. All others are of
  class \code{DataFrame}
}
\details{
  The following changes/ modifications / addition has been made to the
  source material.

  The creation of this object based on the Illumina annotation and
  additional SNP information (see below) is contained in the
  \code{createAnno.R} script in the \code{scripts} directory.

  Contents of the columns are the responsibility of Illumina; see their
  documentation for details.
}

\source{
  See description.
}

\examples{
data(IlluminaHumanMethylation27kanno.ilmn12.hg19)
data(Locations)
data(Other)
data(Manifest)
}
\keyword{datasets}
