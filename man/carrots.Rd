\name{carrots}
\alias{carrots}
\title{
  Consumer preference mapping of carrots
}
\description{
  In a consumer study 103 consumers scored their preference of 12 danish carrot types on a scale from 1 to 7. Moreover the consumers scored the degree of sweetness, bitterness and crispiness in the products. The carrots were harvested in autumn 1996 and tested in march 1997. 
In addition to the consumer survey, the carrot products were evaluated by a trained panel of tasters, the sensory panel, with respect to a number of sensory (taste, odour and texture) properties. Since usually a high number of (correlated) properties(variables) are used, in this case 14, it is a common procedure to use a few, often 2, combined variables that contain as much of the information in the sensory variables as possible. This is achieved by extracting the first two principal components in a principal components analysis(PCA) on the product-by-property panel average data matrix. In this data set the values of the first two principal components are provided (\code{sens1} and \code{sens2}).
}
\usage{
carrots
}
\format{
\describe{
  \item{\code{Consumer}}{
    Numbering identifying consumers
  }
  \item{\code{Frequency}}{
    factor with 5 levels; "How often do you eat carrots?"
    1: once a week or more, 2: once every two weeks, 3: once every three weeks,
    4: at least once month, 5: less than once a month
  }
  \item{\code{Gender}}{
    gender: 1 male, 2 female 
  }
  \item{\code{Age}}{
    factor with 4 levels: 1: less than 25 years, 2: 26-40 years, 3: 41-60 years, 4  more than 61 years
  }
  \item{\code{Homesize}}{
    factor with two levels. Number of persons in the household. 1: 1 or 2 persons, 2: 3 or more persons
  }
  \item{\code{Work}}{
    7 different types of employment. 1: unskilled worker(no education), 2: skilled worker(with education), 3: office worker, 4: housewife (or man), 5: independent businessman/ self-employment, 6: student, 7: retired 
  }
  \item{\code{Income}}{
   1: <150000, 2: 150000-300000, 3: 300000-500000, 4: >500000
  }
}}
\source{
  Per Bruun Brockhoff, The Royal Veterinary and Agricultural University, Denmark.
}

\examples{

data(carrots)
library(lme4)
library(MixMod)

m<-lmer(Preference ~ sens1+sens2+Homesize + Homesize:sens1+Homesize:sens2
+(1+sens1+sens2|Consumer), data=carrots)

#total analysis with 
#automatic elimination of the random and fixed parts
totalAnalysis(m, carrots, isRandReduce=TRUE, isFixReduce=TRUE)
}

\keyword{datasets}
