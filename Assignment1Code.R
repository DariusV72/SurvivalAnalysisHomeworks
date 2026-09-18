set.seed(123)#reproduce results, can be removed for a more random simulation
library(survival)

findCmaxAnyPrc <- function( targetPrc, nCount, lambda){
  if(targetPrc == 0){
    return(1.8*(10^308)) #fun fact, this is about the biggest number R can handle, after that it's INF
  }
  
  timesValues <- rexp(nCount, lambda)#generates a fresh batch of values each run, idk if it's the best choice though
  cmaxValues <- seq(1, 5*max(timesValues), by = 0.1)#make cmax based on the values simulated
  empiricalCensoring <- sapply(cmaxValues, function(cmax){#basically, for each cmax value
    #calculates what percentage of the simulated values would be censored, i.e. event time 
    #smaller than censored time
    cnes <- runif(nCount, 0, cmax)
    mean(timesValues > cnes)
  })
  #based on the percentages calculated, choose the cmax value with the closes percenetage to the target from the list
  #as in, if the target is 50%, the cmax entry where the empirical percentage is 49% or 51% is chosen, no safety for higher or lower
  bestCmax<-cmaxValues[which.min(abs(empiricalCensoring - targetPrc))]
  
  list(cmax = bestCmax)
}


generalSimulationFunction <- function(simulationR = 2000, nCount, targetPrc, targetCmax, lambda){
  
  naiveEstimator <- 0
  mleEstimator <- 0
  
  naiveEstimatorsResults <- vector(mode = "logical", length = simulationR)
  mleEstimatorsResults <- vector(mode = "logical", length = simulationR)
  
  censoringProportions <- vector(mode = "logical", length = simulationR)
  
  rZeroNAcounter <- vector(mode = "logical", length = simulationR)
  
  for (run in 1:simulationR){
    #these are taken straight from the lecture slides, minus the if's
    eventsT <- rexp(nCount, lambda)
    
    if(targetPrc == 0){
      censored <- rep(Inf, nCount)#no censoring, so to make sure there's nothing censored, just use R's biggest value
    }else{
      censored <- runif(nCount, 0, targetCmax)
    }
    
    
    Y <- pmin(eventsT, censored)
    
    delta <- as.integer(eventsT <= censored)
    
    naiveEstimator <- 1/mean(Y)
    
    censoringProportions[run] <-1- mean(delta)
    
    if (sum(delta) == 0){
      mleEstimator <- NA
      rZeroNAcounter[run] <- TRUE
    } else {
      mleEstimator <- sum(delta)/sum(Y)
    }
    
    #fit.exp<-survreg(Surv(Y, delta)~1, dist = "exponential")
    #added this for fun, to see the output, slows down the program a lot so I'd avoid uncommenting it 
    #print(fit.exp)
    
    naiveEstimatorsResults[run] <- naiveEstimator
    mleEstimatorsResults[run] <- mleEstimator
  }
  
  mseNaive <- mean((naiveEstimatorsResults - lambda)^2)
  mseMle <- mean((mleEstimatorsResults - lambda)^2, na.rm = TRUE)
  
  naiveBias <- mean(naiveEstimatorsResults) - lambda
  mleBias <- mean(mleEstimatorsResults, na.rm = TRUE) - lambda
  
  naiveVariance <- mseNaive - naiveBias^2
  mleVariance <- mseMle - mleBias^2
  
  list(censoringPropsAvg = mean(censoringProportions), naiveBias = naiveBias, mleBias = mleBias,
       naiveVariance = naiveVariance, mleVariance = mleVariance, mseNaive = mseNaive,
       mseMle = mseMle, rZeroCases = mean(rZeroNAcounter))
}


censorPercentage1 <- 0.1
censorPercentage2 <- 0.5
censorPercentage3 <- 0


n1 <- 200
n2<- 500

lambda <- 0.05


#censoring dist family = uniform (for now)

times <- rexp(n1, lambda)
times2 <- rexp(n2, lambda)

#find cmaxes, assignment says calibration on 1000 values, so I'll just hardcode that

calibrationCmax10Prc200Entries <- findCmaxAnyPrc(censorPercentage1, 1000, lambda)
calibrationCmax50Prc200Entries <- findCmaxAnyPrc(censorPercentage2, 1000, lambda)
calibrationCmax0Prc200Entries <- findCmaxAnyPrc(censorPercentage3, 1000, lambda)


cmax10PRC<- calibrationCmax10Prc200Entries$cmax #cmax for 10 percentage censoring, 200 values
cmax50PRC<- calibrationCmax50Prc200Entries$cmax #cmax for 50 percentage censoring, 200 values
cmax0PRC <- calibrationCmax0Prc200Entries
#do simulations
simulationR <- 2000

simulation10Prc200Entries <- generalSimulationFunction(simulationR =simulationR, n1, censorPercentage1, 
                            cmax10PRC, lambda)



simulation50Prc200Entries <- generalSimulationFunction(simulationR =simulationR, n1, censorPercentage2, 
                                                       cmax50PRC, lambda)

simulation10Prc500Entries <- generalSimulationFunction(simulationR =simulationR, n2, censorPercentage1, 
                                                       cmax10PRC, lambda)



simulation50Prc500Entries <- generalSimulationFunction(simulationR =simulationR, n2, censorPercentage2, 
                                                       cmax50PRC, lambda)

simulation0Prc200Entries <- generalSimulationFunction(simulationR =simulationR, n1, censorPercentage3, 
                                                       cmax0PRC, lambda)

simulation0Prc500Entries <- generalSimulationFunction(simulationR =simulationR, n2, censorPercentage3, 
                                                       cmax0PRC, lambda)





simulation10Prc200Entries #10%, n = 200
simulation50Prc200Entries #50%, n = 200
simulation0Prc200Entries#0%, n = 200

simulation10Prc500Entries#10%, n = 500
simulation50Prc500Entries#50%, n = 500
simulation0Prc500Entries#0%, n = 500

