### examples form the manual, comparisons with SAS
library(MixMod)
library(lme4)

#### EXAMPLES MixMod package, comparisons with SAS (testMixMod.sas)

## example MixMod-package 
data(ham)
m<-lmer(Informed.liking ~ Product*Information*Kjonn + (1|Consumer), data=ham)
t<-totalAnalysis(m, ham, isRandReduce=TRUE, isFixReduce=TRUE)
#OK with SAS
t<-totalAnalysis(m, ham, plot=TRUE)
#OK with SAS


data(TVbo)
m<-lmer(Cutting ~ TVset*Picture +  (1|Assessor) + (1|Assessor:Picture) + (1|Assessor:TVset), data=TVbo)
t<-totalAnalysis(m, TVbo, test.effs="TVset", isFixReduce=TRUE, isRandReduce=TRUE)
#OK with SAS
anovaTAB(m, TVbo)
# OK with SAS
randTAB(m, TVbo)
# OK
lsmeansTAB(m, TVbo, plot=TRUE)
#OK with SAS


data(carrots)
m<-lmer(Preference~sens2+sens1 + BITTER + Homesize + sens2:Homesize + (1+sens2|Consumer), data=carrots)
randTAB(m, carrots)
#OK
an<-anovaTAB(m, carrots)
#OK with SAS
totalAnalysis(m, carrots, isRandReduce=TRUE, isFixReduce=TRUE)
#OK with SAS
totalAnalysis(m, carrots, isRandReduce=TRUE, isFixReduce=FALSE)
#OK with SAS

## example anovaTAB
library(lme4)
data(ham)
m<-lmer(Informed.liking ~ Product*Information*Kjonn + (1|Product:Consumer) + (1|Information:Consumer) + (1|Consumer), data=ham)
anova(m)
anovaTAB(m, ham)
# OK with SAS


## example lsmeansTAB
library(lme4)
data(TVbo)
m<-lmer(Coloursaturation ~ TVset*Picture + (1|Assessor),  data=TVbo)
lsmeansTAB(m, TVbo, plot=TRUE)
lsmeansTAB(m, TVbo, test.effs="TVset")
#OK with SAS


### example randTAB
library(lme4)
data(carrots)
m<-lmer(Preference~sens1+sens2+Homesize+Homesize:sens1+(1+sens2+sens1|Consumer), data=carrots)
randTAB(m, carrots)
totalAnalysis(m, carrots, isRandReduce=TRUE)$rand.table

### example totalAnalysis
library(lme4)
data(ham)
m<-lmer(Informed.liking ~ Information*Kjonn + (1|Product:Consumer) + (1|Information:Consumer) + (1|Consumer), data=ham)
totalAnalysis(m, ham, isRandReduce=TRUE, isFixReduce=TRUE, plot=TRUE)



data(carrots)
m<-lmer(Preference~sens2+sens1+BITTER+Homesize+Age+sens1:Homesize + sens2:Homesize+sens2:sens1+Homesize:Age+sens2:Homesize:Age+sens1:Homesize:Age+(1+sens2|Consumer), data=carrots)
t<-totalAnalysis(m, carrots, isRandReduce=TRUE, isFixReduce=TRUE, plot=TRUE)
# OK with SAS
t<-totalAnalysis(m, carrots,test.effs="Homesize:Age" , isRandReduce=FALSE, isFixReduce=FALSE, plot=TRUE)
#OK with SAS


#### example missing values
## example MixMod-package 
data(ham)
ham2<-ham
ham2[2,6]<-NA
m<-lmer(Informed.liking ~ Product*Information*Kjonn + (1|Consumer), data=ham2)
t<-totalAnalysis(m, ham2)
ham2[2,2]<-NA
m<-lmer(Informed.liking ~ Product*Information*Kjonn + (1|Consumer), data=ham2)
t<-totalAnalysis(m, ham2)
#OK with SAS