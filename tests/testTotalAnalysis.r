library(MixMod)
library(lme4)


data(ham)
m<-lmer(Informed.liking ~ Age*Product*Information*Kjonn + (1+Age|Consumer), data=ham)
t<-totalAnalysis(m, ham, plot=TRUE, isFixReduce=TRUE, isRandReduce=TRUE)
# OK with SAS
t<-totalAnalysis(m, ham, plot=TRUE)
# OK with SAS
# difference in df for Age and Age:Kjonn


data(TVbo)
m<-lmer(Cutting ~ TVset*Picture +  (1|Repeat) + (1|Assessor) + (1|Repeat:Assessor) + (1|Assessor:Picture) + (1|Assessor:TVset) + (1|Repeat:Assessor:TVset) + (1|Repeat:TVset) + (1|Repeat:Picture), data=TVbo)
t<-totalAnalysis(m, TVbo, plot=TRUE)
t<-totalAnalysis(m, TVbo, test.effs=c("TVset"), isFixReduce=TRUE, isRandReduce=TRUE)
# OK with SAS


m<-lmer(Cutting ~ TVset*Picture - TVset +  (1|Repeat) + (1|Assessor) + (1|Repeat:Assessor) + (1|Assessor:Picture) + (1|Assessor:TVset) + (1|Repeat:Assessor:TVset) + (1|Repeat:TVset) + (1|Repeat:Picture), data=TVbo)
t<-totalAnalysis(m, TVbo, plot=TRUE)
t<-totalAnalysis(m, TVbo, test.effs=c("Picture","Picture:TVset"), isFixReduce=TRUE, isRandReduce=TRUE)
# OK with SAS

for (i in 5:19){
res1=lmer(TVbo[,i]~TVset*Picture +(1|Assessor) + (1|Repeat) + (1|Assessor:Repeat),data=TVbo)
print(totalAnalysis(res1, TVbo,isRandReduce=TRUE, isFixReduce=FALSE, plot=TRUE))
}

