\name{TVbo}
\alias{TVbo}
\title{
  TV dataset
}
\description{
  The \code{TVbo} dataset comes from Bang and Olufsen company. The main purpose was to test products, specified by two attributes Picture and TVset. 15 different response variables (characteristics of the product) were assessed by trained panel list.
  
}
\usage{
TVbo
}
\format{
\describe{
  \item{\code{Assessor}}{
    Numbering identifying assessors
  }
  \item{\code{TVset}}{
    factor - attribute of the product
  }
  \item{\code{Picture}}{
    factor: attribute of the product
  }
  \item{\code{15 Characteristics of the product}}{
    numeric variables: Coloursaturation, Colourbalance, Noise, Depth, Sharpness, Lightlevel, Contrast, Sharpnessofmovement, Flickeringstationary, Flickeringmovement,
 Distortion, Dimglasseffect, Cutting, Flossyedges, Elasticeffect   
  }
  
}}
\source{
  Bang and Olufsen company
}

\examples{

data(TVbo)
library(lme4)
library(MixMod)

m<-lmer(Cutting ~ TVset*Picture+(1|Assessor)
+ (1|Assessor:Picture) + (1|Assessor:TVset), data=TVbo)
totalAnalysis(m, TVbo, test.effs="TVset", isFixReduce=TRUE, isRandReduce=TRUE)
anovaTAB(m, TVbo)
randTAB(m, TVbo)
lsmeansTAB(m, TVbo, plot=TRUE)

}

\keyword{datasets}
