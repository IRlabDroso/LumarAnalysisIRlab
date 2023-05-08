% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{Finding_pics}
\alias{Finding_pics}
\title{Finding_pics}
\usage{
Finding_pics(csv_input, exp_odorant = NULL)
}
\arguments{
\item{csv_input}{rolling_mean table output of \code{\link{format_csv}}.}

\item{exp_odorant}{An odorant name present in the csv_input.}
}
\description{
Find the pics in the data and superpose them for 1 odorant
}
\examples{
csv_df = format_csv(csv)
Finding_pics(csv_df$rolling_mean, exp_odorant="water 1")

}