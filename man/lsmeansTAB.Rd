\name{lsmeansTAB}
\alias{lsmeansTAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculates Least Squares Means and Confidence Intervals for the factors of a fixed part of mixed effects model of lmer object.
}
\description{
 Produces a table which resembles to what SAS software gives in proc mixed statement. 
}
\usage{
lsmeansTAB(model, data, test.effs=NULL, plot=FALSE ,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
   linear mixed effects model (lmer object from lme4 package).
}
\item{data}{
    
     data frame containing the variables named in formula in model.
}

\item{test.effs}{
    
     charachter vector specyfying the names of terms to be tested. If NULL all the terms are tested.
}

\item{plot}{
    
     logical, if TRUE plots are performed: two-way interaction (if present), LSMEANS with confidence intervals
}

  \item{\dots}{
     other potential arguments.
}
}
\details{
 The value of t-statistics is calculated using Satterthwaite's approximation of degrees of freedom.
}
\value{
 Produces Least Squares Means table with p-values and Confidence intervals. 
 }
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{totalAnalysis}}, \code{\link{randTAB}}, \code{\link{anovaTAB}}
}
\examples{
library(lme4)
data(ham)
m1<-lmer(Informed.liking ~ Kjonn*Information +(1|Consumer), data=ham)
lsmeansTAB(m1, ham, plot=TRUE, test.effs="Kjonn:Information")

data(TVbo)
m<-lmer(Coloursaturation ~ TVset*Picture + (1|Assessor), data=TVbo)
lsmeansTAB(m, TVbo, plot=TRUE)
lsmeansTAB(m, TVbo, test.effs="TVset")
}
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
