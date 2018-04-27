\name{srtm}
\alias{srtm.download}
\alias{srtm.list}
\alias{srtm.region}
\alias{srtm.slr}
\title{
SRTM Package
}
\description{
Tool for handling srtm-data including cca [1]
}
\usage{
srtm.download(url=3, all=TRUE, dest=getwd(),n,s,e,w)
srtm.list()
srtm.region(e,w,s,n,path=getwd(),del.tmp=TRUE,plot=FALSE,out="output",type="")
srtm.slr(srtm,slr,pos=FALSE,type="",out="output",cca.s=2,mode=3,plot=FALSE)
}
\arguments{

  \item{all}{
 If TRUE all Files which are available, will be downloaded. If FALSE you can specify a single region using n,s,e,w.
}
  \item{cca.s}{
 The radius/shell size of the burning procedure in the cca
(i.e. how tolerant to small gaps the algorithm is)
}
  \item{del.tmp}{
 If TRUE the created temporary Folder, which contains the unzipped srtm files will be deleted.
}
  \item{dest}{
 Defining a download folder.
}
  \item{e}{
 east coordinate in decimal degree (e=12.2)
}
  \item{mode}{
The algorithm comes in three versions which affect which close cells are 
included to the considered cluster:
(mode=1) nearest neighbors
(mode=2) cells within a shell (i.e. squares of certain size)
(mode=3) cells within a radius
Whereas (mode=1) is equivalent to (mode=3) with r=1 and 
(mode=2) with r=1 is equivalent to (mode=3) with r=2.
}
  \item{n}{
 north coordinate in decimal degree (n=54.2)
}
  \item{out}{
 Defining the name of the output file. Only active if you have chosen a type.
}
  \item{path}{
 Defining the path to the folder with downloaded srtm zip-files.
}
  \item{plot}{
 If TRUE a plot of the extracted srtm will be generated.
}
  \item{pos}{
 A vector containing the starting position for the clustering. Will automatically reseted, if the starting point is not defined as ocean.
}
  \item{s}{
 south coordinate in decimal degree (s=53.2)
}
  \item{slr}{
 Defining a Sea-level in meter for flooding.
}
  \item{srtm}{
 A R raster object. All internal operations will be done on such a raster, excluding the cca (ordinary matrix).
}
  \item{type}{
 Defining the type of the output file. You can choose between "shp" and "tif". If type="", as per default, no output file will be generated.
}

  \item{url}{
 Defining alternative download server as String. Per default you can choose between:
 (url=1) "ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/"
 (url=2) "http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/"
 (url=3) "http://droppr.org/srtm/v4.1/6_5x5_TIFs/"
 (url=4) "ftp://xftp.jrc.it/pub/srtmV4/tiff/"
}
\item{w}{
 west coordinate in decimal degree (w=13.2)
}
}
\details{
keine Ahnung
}
\value{
Was ist Value?
}

\references{
[1] Rozenfeld HD, et al. (2008) Laws of population growth. Proc Nat Acad Sci USA 105:18702-18707.
}
\author{
Steffen Kriewald
}



\keyword{ utils }

\examples{

n<-32.9
s<-32.6
w<--17.3
e<--16.6

srtm.download(all=FALSE,n=n,s=s,w=w,e=e)
madeira <- srtm.region(n=n,s=s,w=w,e=e,plot=TRUE)
slr <- srtm.slr(madeira,slr=100,cca.s=1,mode=2,plot=TRUE)
}
