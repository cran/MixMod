\name{randTAB}
\alias{randTAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs Chi square test on random effects of linear mixed effects model. 
}
\description{
Returns a table with  value of Chi square test, p-value, elimination number for each random effect
}
\usage{
randTAB(model, data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{model}{
    
     linear mixed effects model (lmer object from lme4 package).
}

\item{data}{
    
     data frame containing the variables named in formula in model.
}

  \item{\dots}{
     other potential arguments.
}
}

\value{
Produces a table which tests whether the random terms are significant.
}

\details{
 The columns of the table are:
 
 
 Chisq: The value of the chi square test
 
 Chi Df: The degrees of freedom for the Chi square test
 
 elim.num: The order in which effect has been eliminated from the model. 0 means it has not been eliminated (either due to it's significance or due to the fact that the reduction was not required). Random effects that have 0 variance or correlation equal 1, -1 are always eliminated.
 
  p.value: The p-value of the Chi square test for the effect

}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff
}
\note{
If the effect has random slopes, then first the correlations between itercept and slopes are checked for significance
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{totalAnalysis}}, \code{\link{anovaTAB}}, \code{\link{lsmeansTAB}}
}
\examples{
options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
library(lme4)
data(carrots)

#lmer model with correlation between intercept and slopes
#in the random part
m<-lmer(Preference~sens1+sens2+Homesize+Homesize:sens1+(1+sens2+sens1|Consumer), data=carrots)

# table with  p-values for the random effects
randTAB(m, carrots)
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}