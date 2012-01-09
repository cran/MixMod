totalAnovaRandLsmeans<-function(model, data,  alpha.rand = 0.05, alpha.fix = 0.05, isFixReduce = FALSE, isRandReduce = FALSE, isTotal=FALSE, isAnova=FALSE, isRand=FALSE, isLSMEANS=FALSE, test.effs=NULL, plot=FALSE)
{
  
  
  #not to show the warnings  
  options(warn=-1)    
    
  result<-NULL
  anova.table<-NULL
  
  
  
  #model<-lmer(formula=formula, data=data)
  
  #update model
  # change unordered contrasts to contr.SAS
  # change REML to TRUE
  options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
  
  #model<-update(model, REML=TRUE)
  mf.final<-update.formula(formula(model),formula(model))
  
   
  
  #update data 
  #eliminate missing values
  if(!paste(mf.final)[2] %in% names(data))
  {
     mframe<-model.frame(mf.final, data=data, na.action=na.pass)
     data$response<-mframe[,1]
     fm<-paste(mf.final)
     fm[2]<-"response"
     mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
     mf.final<-update.formula(mf.final,mf.final)
  }
  
  data<-get_all_vars(mf.final, data)
  data<-data[complete.cases(data),]
  
  
 
  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  
  
  
    # check if there are no fixed effects
  #if(length(attr(delete.response(terms(model)),"term.labels"))==0)
 # {
  #  fm<-getFormula(model, withRand=TRUE)
  #  fm<-as.formula(paste(fm[2],fm[1],paste("1+", fm[3], sep=""), sep=""))
  #  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  #}
   
  
   
  # save the call of the model              
  result$call<-model@call
  
  #check if there are correlations between intercept and slope
  result$corr.intsl<-checkCorr(model)
  
  
  
  
  #Reduction of random effects
  #if(isRand || isTotal )
  #{     
    # perform reduction of random effects
    #if(!isRandReduce)
    #  result.rand<-elimRandLmer(model, data, 1)
    #else
    #  result.rand<-elimRandLmer(model, data, alpha.rand)
    if(!isRandReduce)
      result.rand<-elimRandEffs(model, data, 1)
    else
      result.rand<-elimRandEffs(model, data, alpha.rand)  
   
   
    model<-result.rand$model
    result$rand.table<-result.rand$TAB.rand   
    if(isRand)
      return(result)
    
  #}
  
  
  
  #perform reduction of fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
  {
    result$anova.table<-anova(model)
    result$model<-model
    lsmeans.summ<- matrix(ncol=7,nrow=0)
    colnames(lsmeans.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
    result$lsmeans.table<-lsmeans.summ
    return(result)
  }
    
  #perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova<-TRUE
  is.first.sign<-TRUE
  
  
  while(!stop)
  {
  
      
    
      # if there are no fixed terms
      if(nrow(anova(model))==0)
      {
        if(is.null(anova.table))
        {
          
          if(isLSMEANS)
          {
             lsmeans.summ<- matrix(ncol=7,nrow=0)
             colnames(lsmeans.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
             result$lsmeans.table<-lsmeans.summ
             return(result)
          }
          result$model<-model
          result$anova.table<-anova(model)
          return(result)        
          
        }          
        break
      }
        
          
      
      # save lmer outcome in rho environmental variable
      rho<-rhoInit(model)     
      
      # calculate asymptotic covariance matrix A
      h <- hessian(function(x) Dev(rho,x), rho$param$vec.matr)
      rho$A <- 2*solve(h)
      
      #Check if A is positive-definite
      isposA<-all(eigen(rho$A)$values>0)
      if(!isposA)
      {
        print("ERROR: asymptotic covariance matrix A is not positive!")
        result$model<-model
        TABs<-emptyAnovaLsmeansTAB()
        result$anova.table<-TABs$TAB.fixed
        result$lsmeans.table<-TABs$TAB.lsmeans
        return(result)
      }
      
      #calculate lsmeans of the final model
      if(isLSMEANS)
      {
        lsmeans.tab<-calcLSMEANS(model, data, rho, alpha.fix, test.effs=test.effs, plot=plot)
        result$lsmeans.table<-lsmeans.tab
        return(result)
      }
      
          
      # Calculate  F-test with Satterthwaite's approximation
      # create X design matrix for fixed effects
      X.design<-createDesignMat(model,data)
       
      
      # define full coefficients for model
      #coefs.real<-fullCoefs(model, data)
      #coefs.real<-getDummyCoefs(model, data)
      
      
      #save full coefficients in rho
      #rho$s.test<-coefs.real
      nums.dummy.coefs<-getNumsDummyCoefs(model, data)
      rho$nums.zeroCoefs<-nums.dummy.coefs$nums.zeroCoefs
      rho$nums.Coefs<-nums.dummy.coefs$nums.Coefs
      
     
      
      #define the terms that are to be tested
      test.terms<-attr(terms(model),"term.labels")
      
      #initialize anova table
      if(is.first.anova)
      {
        anova.table<-initAnovaTable(model, isFixReduce)
        is.first.anova<-FALSE
        elim.num<-1
      }
      
      # calculate general set matrix for type 3 hypothesis
      L<-calcGeneralSetForHypothesis(X.design, rho)
      
      for(i in 1:length(test.terms))
      {
        
         result.fstat<-calcFpvalue(test.terms[i], L, model, rho)   
         if(!is.null(result.fstat))
         {
            anova.table[which(rownames(anova.table) == test.terms[i]),2]<-round(result.fstat$denom,2)
            anova.table[which(rownames(anova.table) == test.terms[i]),3]<-round(result.fstat$Fstat,2)
            anova.table[which(rownames(anova.table) == test.terms[i]),which(colnames(anova.table)=="p.value")]<-round(result.fstat$pvalue,4)
            #if(!is.first.anova)
            #{
            #  anova.model<-anova(model)
            #  anova.table[which(rownames(anova.table) == test.terms[i]),1:4]<-anova.model[which(rownames(anova.model) == test.terms[i]),]
            #}
              
         }
       
      }
     
      if(!isFixReduce)
        break
      else
      {
        resNSelim<-elimNSFixedTerm(model, anova.table, data, alpha.fix, elim.num)
        if(is.null(resNSelim))
          break
        else
        {
          model<-resNSelim$model
          mf.final<-update.formula(formula(model),formula(model))
          model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
          anova.table<-resNSelim$anova.table
          elim.num<-elim.num+1
        }        
      }      
    
  }
  
  
  if(isTotal || isAnova)
  {
    result$anova.table<-anova.table
    if(isAnova)
      return(result)
  }
  
  lsmeans.tab<-calcLSMEANS(model, data, rho, alpha.fix, test.effs=test.effs, plot=plot)
  result$lsmeans.table<-lsmeans.tab
   
  #update model
  mf.final<-update.formula(formula(model),formula(model))
  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  
  #save model
  result$model<-model
  return(result)
}


# generic functions for total analysis on mixed models
#totalAnalysis <- function(model, data,  alpha.rand = 0.05, alpha.fix = 0.05, isFixReduce = FALSE, isRandReduce = FALSE, test.effs=NULL, plot=FALSE, ...) UseMethod("totalAnalysis")


totalAnalysis<-function(model, data,  alpha.rand = 0.05, alpha.fix = 0.05, isFixReduce = FALSE, isRandReduce = FALSE, test.effs=NULL, plot=FALSE, ...)
{  
  result<-totalAnovaRandLsmeans(model, data,  alpha.rand, alpha.fix, isFixReduce, isRandReduce, isTotal=TRUE, test.effs=test.effs, plot=plot)
  class(result)<-"totalAnalysis"
  result
}

### UNUSED function
#totalAnalysis.formula <- function(formula, data, ...)
#{
#  model <- lmer(formula=formula, data=data)
#  resAnalysis<-totalAnalysis.default(model, data, ...)
#  resAnalysis$call<-match.call()
#  resAnalysis
#}

print.totalAnalysis <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nRandom effects:\n")
  printCoefmat(data.matrix(x$rand.table), digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)        
  if(nrow(x$anova.table)!=0)
  {
    if(class(x$model) == "lm" | class(x$model) == "gls")
    {
      cat("\nFixed effects:\n")
      print(x$anova.table)
      cat("\nLeast squares means:\n")
      print(x$lsmeans.table) 
      cat("\nFinal model:\n")
      print(x$model)
      return()
    }
    else
    {
      cat("\nFixed effects:\n")
      printCoefmat(data.matrix(x$anova.table), dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      cat("\nLeast squares means:\n")
      printCoefmat(data.matrix(x$lsmeans.table), dig.tst=1  ,tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
    }    
  }
  else
    print(x$anova.table)
  cat("\nFinal model:\n")
  print(x$model@call)
  
  
}


# generic functions for ANOVA
#anovaTAB <- function(model, data, ...) UseMethod("anovaTAB")


#anovaTAB.default<-function(model, data, ...)
anovaTAB<-function(model, data, ...)
{
  result<-totalAnovaRandLsmeans(model, data, isAnova=TRUE)  
  res <- list(anova.table=result$anova.table)
  class(res)<-"anovaTAB"
  res
}

print.anovaTAB <- function(x, ...)
{

  if(nrow(x$anova.table)!=0)
  {
    cat("Analysis of Variance Table:\n")
    printCoefmat(data.matrix(x$anova.table), dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 , P.values=TRUE, has.Pvalue=TRUE)
  }
  else
    print(x$anova.table)
}


# generic functions for random effects
#randTAB <- function(model, data, ...) UseMethod("randTAB")


#randTAB.default<-function(model, data, ...)
randTAB<-function(model, data, ...)
{
  result<-totalAnovaRandLsmeans(model, data, isRand=TRUE)  
  res <- list(rand.table=result$rand.table, isCorr = result$corr.intsl)
  class(res)<-"randTAB"
  res
}

print.randTAB <- function(x, ...)
{

  cat("Analysis of Random effects Table:\n")
  #if(x$isCorr)
  #  print(x$rand.table)
  #else
    printCoefmat(data.matrix(x$rand.table), digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)        
}



# generic functions for LSMEANS
#lsmeansTAB <- function(model, data, test.effs=NULL, plot=FALSE , ...) UseMethod("lsmeansTAB")


#lsmeansTAB.default<-function(model, data, test.effs=NULL, plot=FALSE , ...)
lsmeansTAB<-function(model, data, test.effs=NULL, plot=FALSE , ...)
{
  result<-totalAnovaRandLsmeans(model, data, isLSMEANS=TRUE, test.effs=test.effs, plot=plot)  
  res <- list(lsmeans.table=result$lsmeans.table)
  class(res)<-"lsmeansTAB"
  res
}

print.lsmeansTAB <- function(x, ...)
{

  cat("Least Squares Means table:\n")
  printCoefmat(data.matrix(x$lsmeans.table),dig.tst=1  ,tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 , P.values=TRUE, has.Pvalue=TRUE)
        
}
