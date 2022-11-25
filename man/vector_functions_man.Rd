\name{multivarvector}
\alias{vector_functions}
\title{
Calculate angles between vectors and mean vector lengths
}
\description{
This package contains two main components that allow calculating the angles
between vectors as well as mean vector lengths to quantify parallelism of microbial community change.
Angles can be calculated between pairs of host populations based on ordination files (PCA/PCoA/NMDS) with varying numbers of axes.
}
\usage{
multivariate.vector(x)
}
\arguments{
  \item{x}{
A dataframe with information on host populations and axes scores obtained from some ordination method
}
}
\details{
The first column of the dataframe must specifiy the sample id, the second column must specifiy information
on the identity of a host population pair (e.g., lake), the third column must specificy the host groups between
which vectors are calculated (e.g., ecotype), the first column must combine the information of the second & third column separated by an underscore (i.e., lake_ecotype). From the fifth column onwards, all columns must contain scores for the different axes.
}
\value{
multivariate.vector returns angles between vectors and mean vector lengths for host population pairs.
}
\references{
}
\author{
Andreas Härer <andreas_haerer@gmx.de>, Diana Rennison <drennison@ucsd.edu>
}
\note{
}
\seealso{
}
\examples{
library("dplyr")
library("plyr")
library("tibble")
library("data.table")

load("~/multivariatevector/data/stickleback_example_data.rda")
multivariate.vector(stickleback_example_data)
}
