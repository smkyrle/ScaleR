\name{ScaleR}
\alias{ScaleR}
\title{ScaleR}
\usage{
ScaleR(x,y,inter, method, plot, seed, k)
}
\description{
Scales data by desired scaling factor (inter), and performs either Ridge regression or PLS.
}
\examples{
ScaleR(x,y,method='PLS',inter=0.1, plot=TRUE, seed=111, k=3)
}
