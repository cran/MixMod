\name{totalAnalysis}
\alias{totalAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs complete analysis of linear mixed effects model:
}
\description{
calculates p-values for the fixed effects (F test based on Sattethwaite's approximation), p-values for the random effects (based on Chi square test),  LSMEANS for the fixed part of the model, performs automatic elimination of the 
random part, automatic elimination of the fixed part. All analysis may be performed on lmer object of lme4 package. 
}
\usage{
totalAnalysis(model, data,  alpha.rand = 0.05, 
alpha.fix = 0.05, isFixReduce = FALSE, 
isRandReduce = FALSE, test.effs=NULL, plot=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{model}{
    
     linear mixed effects model (lmer object from lme4 package).
}

\item{data}{
    
     data frame containing the variables named in formula in model.
}

\item{alpha.rand}{
    
     significance level for elimination of the random part
}

\item{alpha.fix}{
    
     significance level for elimination of the fixed part
}

\item{isFixReduce}{
    
    logical for whether the reduction of the fixed part is required
}

\item{isRandReduce}{
    
    logical for whether the reduction of the random part is required
}

\item{test.effs}{
    
    charachter vector specyfying the names of terms to be tested in LSMEANS. If NULL all the terms are tested.
}

\item{plot}{
    
     logical, if TRUE plots are performed: two-way interaction (if present), LSMEANS with confidence intervals
}

  \item{\dots}{
     other potential arguments.
}
}
\details{
 First the random part is reduced (if possible and required), based on Chi square test and REML=TRUE, then the fixed
 part is eliminated (if possible and required) in a way: the highest order interactions are tested first: if they are significant, the main effects are not tested for significance.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{rand.table}{table with value of Chi square test, p-values e t.c. for the random effects}
  \item{anova.table }{table which tests whether the model fixed terms are significant (Analysis of Variance)}
  \item{lsmeans.table}{Least Squares Means table with p-values and Confidence intervals}
  \item{model}{Final model - object of class lmer or gls (after all the required reduction has been performed)}
  }
%% ...


\author{
Alexandra Kuznetsova, Per Bruun Brockhoff
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{anovaTAB}}, \code{\link{randTAB}}, \code{\link{lsmeansTAB}}
}
\examples{

library(lme4)
data(carrots)

#construction of lmer model, 
#with correlation between slope and
#intercept
m2<-lmer(Preference~sens2+sens1 + BITTER + Homesize + Age+sens1:Homesize
+ sens2:Homesize + sens2:sens1 + Homesize:Age+sens2:Homesize:Age
+ sens1:Homesize:Age + (1+sens2|Consumer), data=carrots)


#reduction only of the random part
totalAnalysis(m2, carrots, isRandReduce=TRUE, isFixReduce=FALSE, plot=TRUE)




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}