\name{MixMod-package}
\alias{MixMod-package}
\alias{MixMod}
\docType{package}
\title{
 Analysis of Mixed Effect Models
}
\description{
This package performs different kinds of analysis on lmer object(of lme4 package), such as automatic elimination of the random part (based on Chi square test), automatic elimination of the fixed part, Analysis of variance (with p-values calculated based on Satterthwaite's approximation of denominator degrees of freedom, Least squares Means table  with confidence intervals for the fixed part of the model).
}
\details{
\tabular{ll}{
Package: \tab MixMod\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-10\cr
License: \tab GPL\cr
}

The calculation of statistics for the fixed part was developed according to SAS Proc Mixed Theory. Therefore 
the contrasts matrix for the fixed effects should be contr.SAS. The option for the contrasts is automatically changed, and the model is updated. 
In the elimination of the random part the models with REML=TRUE are compared, therefore the logical argument REML of lmer object should be the default one REML=TRUE.

}
\author{
Alexandra Kuznetsova <alku@imm.dtu.dk>, Per Bruun Brockhoff 
}
\references{
Theory SAS proc mixed (for calculation of denominator degrees of freedom)
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%%~~ directory ~~
%%\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
\examples{
options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
library(lme4)
data(ham)
m<-lmer(Informed.liking ~ Kjonn*Information*Product +(1|Consumer), data=ham)
totalAnalysis(m, ham,isRandReduce=TRUE, isFixReduce=TRUE)

data(TVbo)
library(lme4)
m<-lmer(Cutting ~ TVset*Picture+(1|Assessor)
+ (1|Assessor:Picture) + (1|Assessor:TVset), data=TVbo)
totalAnalysis(m, TVbo, test.effs="TVset", isFixReduce=TRUE, isRandReduce=TRUE)
anovaTAB(m, TVbo)
randTAB(m, TVbo)
lsmeansTAB(m, TVbo, plot=TRUE)

library(lme4)
data(carrots)
m2<-lmer(Preference~sens2+sens1+BITTER+Homesize+
sens2:Homesize+(1+sens2|Consumer), data=carrots)
randTAB(m2, carrots)
an<-anovaTAB(m2, carrots)
totalAnalysis(m2, carrots, isRandReduce=TRUE, isFixReduce=TRUE)

}
