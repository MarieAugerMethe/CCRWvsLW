\name{movFD}
\alias{movFD}
\docType{package}
\title{
Transform locations into appropriate measures for movement analysis.
}
\description{
This functiontakes the locations creates with the simulation loops and calculates the step length and turning angle as weel as assessing the missing locations.
}
\usage{
movFD(xy, tl)	
}
\arguments{
  \item{xy}{
    matrix: locations of the movement path
  }
  \item{tl}{
    integer: length of the individual movement path (needed for the rediscretised path)
  }
}
