\name{raster.to.png}

\title{
Sample Package
}
\description{
Plot raster as png
  png.path <- getOption("png.path",default=getwd())}

\usage{
raster.to.png(raster, name="test", mypal=c("wheat1", "red3","green"),n=10, style="equal", fixed=NULL,trans=1,legend=FALSE)
}
\arguments{

  \item{r1}{The coarse raster which gives the final extent and resolution}
  \item{r2}{The fine raster which will be re-sampled}

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
raster <- raster(nrows=360,ncols=720,xmn=-180, xmx=180, ymn=-90, ymx=90)
raster[] <- 1:(360*720)
raster.to.png(raster)
}