\name{sample}
\alias{sample}
\alias{coord}
\title{
Sample Package
}
\description{
The sample-function computes the exact amount of a high resolution grid to a low resolution grid.
The algorithm selects step by step each grid cell of the dataset with the lower resolution and sums up the area, or optional the cell value based on the dataset with the higher resolution. High resolution cells, which belong to multiple low grid cells, will be split up proportionally. The algorithm is programmed for lat-long WGS84 and uses the coord-function for computing the coordinates of each cell. 
}
\usage{
sample(r1,r2,cell.constraint=0,area=TRUE,crop=TRUE)
coord(r1, vector=length(r1[]))
}
\arguments{

  \item{r1}{The coarse raster which gives the final extent and resolution}
  \item{r2}{The fine raster which will be re-sampled}
  \item{cell.constraint}{Consider only cells with a higher value as the cell.constraint}
  \item{area}{If TRUE the amount of the area of the high resolution grid cells will be computed for each low resolution grid cell. If FALSE the grid cell values are summed up}
  \item{crop}{If TRUE r2 will be crop to the extent of r1. If FALSE the lapping cells will be proportionately considered. (Exception: The last lapping cell in the first row and the first lapping cell in the last row are not considered, due to a structural issues.)}
  \item{vector}{Gives the cell id of the raster r1 for which the coordinates should be computed. Default is all.}
}



\details{
keine Ahnung
}
\value{
Was ist Value?
}

\references{
[1] Not yet
}
\author{
Steffen Kriewald
}



\keyword{ utils }

\examples{
r1 <- raster(nrows=2, ncols=3, xmn=10, xmx=28, ymn=-2, ymx=10)
r1[] <- 1:6

r2 <- raster(nrows=3, ncols=4, xmn=11, xmx=19, ymn=3, ymx=9)
r2[] <- 0
r2[c(1,2,3,7,11)] <- 1:5

plot(r1)
plot(r2, add=T)

sample_area <- sample(r1,r2)
plot(sample_area[[2]])

sample <- sample(r1,r2,area=FALSE)
plot(sample[[2]])
}