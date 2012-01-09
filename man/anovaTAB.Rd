\name{anovaTAB}
\alias{anovaTAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculates Analysis of Variance Table for mixed model of object of class lmer.
}
\description{
 Performs type 3 hypothesis (definition from SAS Theory) for fixed effects part of a linear mixed effects model.
 The F statistics is calculated based on Satterthwaite's approximation for denominator degrees of freedom. 
}
\usage{
anovaTAB(model, data, ...)
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
\details{
 The definition of type 3 hypothesis is taken from SAS software, the contrast matrix is calculated according to SAS theory Proc GLM. The options for the contrasts should be "contr.SAS".
}
\value{
 Produces a table which tests whether the model fixed terms are significant. The test is F test with Satterthwaite's approximation of degrees of freedom.  
}
\references{
Theory SAS Proc GLM
}
\author{
 Alexandra Kuznetsova, Per Bruun Brockhoff
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{totalAnalysis}}, \code{\link{randTAB}}, \code{\link{lsmeansTAB}}
}
\examples{

library(lme4)
data(ham)
m1<-lmer(Informed.liking ~ Kjonn*Information*Product + (1|Product:Consumer) + 
(1|Information:Consumer) + (1|Consumer), data=ham)
# anova table of lme4 package
anova(m1)
# anova table produced by MixMod package
anovaTAB(m1, ham)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
