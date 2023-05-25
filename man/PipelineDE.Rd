\name{PipelineDE}
\alias{PipelineDE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Differentially Expression analysis with a R Shiny graphic-user interface
}
\description{
Launches in the default browser an R Shiny app that allows to identify differentially expressed genes with edgeR and DESeq2, immediately after the alignment of the reads
}
\usage{
PipelineDE(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{}
}
}
\value{
The results of the analysis can be downloaded from the app in the tsv file format. The output is a data frame, as described with the function \code{\link[=DEinr]{DEinr()}}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=DEinr]{DEinr()}}: the function called by the app.
}
\examples{
PipelineDE()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
