estimatePatterns <- function(patternCounts, 
							 epsilon = 0,
							 eta = 0, 
							 column = -1, 
							 fast = TRUE, 
							 steps = 20000, 
							 plot = TRUE,
							 yLimit1 = -1, 
							 yLimit2 = -1)							 
{	
	if(column == -1){
		columns <- 1:(ncol(patternCounts)-1)
		} else {
		columns <- column
		}

	compareData <- list()
	for(i in columns){
			compareData[[i]] <- estimatePatternsS(patternCounts, epsilon, eta, column = i, fast, steps)
	}		
	
	if(plot){	
		pdf("patternDistribution.pdf")
		for (i in columns){
		 	plotGraph(compareData[[i]], yLimit1, yLimit2)
		 	}
		dev.off()		
		} 
	
	if(length(columns)==1){
		compareData <- compareData[[1]]
	}
		
	return(compareData)
		
}

estimatePatternsS<- function(patternCounts, epsilon, eta, column, fast, steps)
{

    # Arguments checks

    if (!is.numeric(epsilon)) {
        stop("epsilon must be numeric\n")
    }

    if (epsilon < 0 || epsilon >= 1) {
        stop("epsilon must be between 0 and 1\n")
    }

    if (!is.element(column, 1:(ncol(patternCounts) - 1))) {
        stop("column beyond allowed range\n")
    }

    options(scipen=999)

    # Read the data

    patternCounts <- subset(patternCounts, select=c(1, column+1))
    names(patternCounts) <- c("Patterns","Counts")
    patternCounts <- patternCounts[patternCounts$Counts!=0,]
    ### FIXME implement proper checks of the patterns format
    trim <- function (x) gsub("[#%$*()&a-zA-Z ]", "", x) #removes all extra symbols.
    patternCounts$Patterns <- trim(as.character(patternCounts$Patterns))
    nCpGsites <- nchar(patternCounts$Patterns[1])

    # Eta input

    if (length(eta) == 1) {
        if (eta < 1 && eta >= 0) {
            eta <- rep(eta, nCpGsites)
        } else {
            stop("eta must be between 0 and 1\n")
        }
    } else {
        if (length(eta) != nCpGsites) {
            stop("length of eta is not equal to the number of CpG sites\n")
        } else if (any(eta < 0 || eta >= 1)) {
            stop("eta must be between 0 and 1\n")
        }
    }

    # Create the vector of patterns

    if (fast) {
        methData <- patternCounts
    } else {
            binary <- function(x) if (all(x < 2)) x else cbind(binary(x %/% 2), x %% 2)
            cytosineBinary <- binary(0:(2 ^ nCpGsites - 1))
            mPattern <- array(dim=c(2^nCpGsites,1))
            for(i in 1:(2 ^ nCpGsites)){
                mPattern[i, 1] <- paste(cytosineBinary[i, ], collapse="")
            }
            counts <- array(0, dim=c(2^nCpGsites,1))
            methData <- data.frame(Patterns=mPattern, Counts=counts)
            methData$Patterns <- as.character(methData$Patterns)
            for(i in 1:(2 ^ nCpGsites)){
                patternMatches <- patternCounts$Patterns == mPattern[i, 1]
                if(any(patternMatches)) {
                    methData[i, 2] <- patternCounts[patternMatches, 2]
                }
            }
    }

    yPatterns <- methData$Counts
    totalPatterns <- sum(yPatterns)
    readDistribution <- yPatterns / totalPatterns
    patternsMax <- which.max(yPatterns)
    yWithoutMax <- yPatterns[-patternsMax]

    # Create the pattern array

    size <- nrow(methData)

     stringToVector <- function(x) {
         patternVector <- array(dim=nchar(x))
         for (i in 1:nchar(x)) {
             patternVector[i] <- as.numeric(substr(x, i, i))
         }
         return(patternVector)
     }

     patternArray <- array(dim=c(size, nCpGsites))

     for (i in 1:size) {
      patternArray[i, ] <- stringToVector(methData$Patterns[i])
     }

     # Construct the conversion matrix

    if (fast) {
        conversionMatrix <- array(0, dim=c(size, size))
        for (i in 1:size) {
            for (j in 1:size) {
                fromIndex <- patternArray[i, ]
                toIndex <- patternArray[j, ]
                conversionMatrix[i, j] <- prod((1 - epsilon - eta + 2 * epsilon * eta) ^ ((fromIndex + toIndex == 0) * 1) *
                                               (epsilon + eta - 2 * epsilon * eta) ^ ((fromIndex - toIndex == -1) * 1) *
                                               eta ^ ((fromIndex - toIndex == 1) * 1) *
                                               (1 - eta) ^ ((fromIndex + toIndex == 2) * 1))
             }
        }
    } else {
        conversionRule <- list(array(0, dim=c(2, 2)))
        for (i in 1:nCpGsites) {
          conversionRule[[i]] <- array(c(1 - epsilon - eta[i] + 2 * epsilon * eta[i],
                                         eta[i],
                                         epsilon + eta[i]- 2 * epsilon * eta[i],
                                         1 - eta[i]),
                                       dim=c(2,2))
        }
        conversionMatrix <- 1
        for (i in 1:nCpGsites) {
            conversionMatrix <- kronecker(conversionMatrix, conversionRule[[i]])
        }
    }


    # Define the likelihood function
    # Note that likelihoodOpt's argument has one less entry than the function likelihood

    likelihood <- function(theta) {
        phi <- theta %*% conversionMatrix
        if (fast) {
            likelihood <- -sum(yPatterns * log(as.vector(phi)))
        } else {
            likelihood <- -sum(yPatterns[yPatterns != 0] * log(phi[yPatterns != 0]))
        }
        return(likelihood)
    }

    expand <- function(theta, patternsMax) {
        expanded <- c(theta[(1:length(theta)) < patternsMax],
                      1 - sum(theta),
                      theta[(1:length(theta)) >= patternsMax])
        return(expanded)
    }

    likelihoodOpt <- function(theta){
        likelihoodOpted <- likelihood(expand(theta, patternsMax))
        return(likelihoodOpted)
    }

    # Optimisation

    startingVector <- yWithoutMax / totalPatterns
    yZeros <- which(startingVector %in% c(0))
    startingVector[yZeros]<- yPatterns[patternsMax] / (100000 * totalPatterns)

    constraintMatrix <- rbind(diag(size - 1), rep(-1, size - 1))
    constraintVector <- append(rep(0, size - 1), -1)
    opt <- constrOptim(startingVector, likelihoodOpt, grad=NULL,
                        ui=constraintMatrix,
                        ci=constraintVector,
                        method="Nelder-Mead",
                        control=list(maxit=steps))

    recovered <- expand(opt$par, patternsMax)

    # Decide the "0's"

    minZeros <- recovered
    minLikelihood <- likelihood(recovered)
    listReduced <- 1:size

    for (i in listReduced[listReduced != patternsMax]) {
        minTemp <- minZeros
        minTemp[patternsMax] <- minTemp[patternsMax] + minTemp[i]
        minTemp[i] <- 0
        if(likelihood(minTemp) <= minLikelihood) {
            minZeros <- minTemp
            minLikelihood <- likelihood(minTemp)
        }
    }

    # Generate output

    compareData <- data.frame(Pattern=methData$Patterns,
                              Coverage=yPatterns,
                              observedDistribution=readDistribution,
                              estimatedDistribution=minZeros)

    if(!fast) {
        compareData <- compareData[compareData$observedDistribution != 0 | compareData$estimatedDistribution != 0, ]
        rownames(compareData) <- 1:nrow(compareData)
    }

    if(opt$convergence!=0) {
        warning('Constrained optimisation did not converge.\n')
    }

    compareData$spurious <- compareData$observedDistribution !=0 & compareData$estimatedDistribution == 0	

    return(compareData)
}

# Plot graphs.

plotGraph <- function(compareData, yLimit1, yLimit2)
{
      if(yLimit1 ==-1 ){
        yLimit1 <- ceiling(max(compareData$observedDistribution, compareData$estimatedDistribution) * 10.2) / 10
      }
      if(yLimit2 ==-1 ){
        yLimit2 <- quantile(sort(compareData$observedDistribution), 0.9)
      }

      layout(matrix(c(1, 2, 3), 3, 1, byrow=TRUE), heights=c(5, 1, 5))

      par(mar=c(5.4, 5, 2, 5))
      plot(compareData$estimatedDistribution,
           pch=4,
           col='white',
           xlab="pattern",
           ylab="proportion",
           cex.lab=1.5,
           ylim=c(0, yLimit1))
      points(compareData$estimatedDistribution, pch=4, col='red')
      par(new=TRUE)
      plot(compareData$Coverage,
           col='blue',
           pch=3,
           xlab="",
           ylab="",
           ylim=c(0, sum(compareData$Coverage) * yLimit1),
           axes=FALSE)
      axis(side=4)
      mtext("coverage", side=4, line=2.5)
      abline(0, 0, col="grey")

      par(mar=c(0, 0, 0, 5))
      plot.new()
      legend('right', c("observed distribution","estimated distribution"), pch=c(3,4), col=c("BLUE", "RED"))

      par(mar=c(5.4, 5, 2, 5))
      plot(compareData$estimatedDistribution,
           pch=4,
           col='white',
           xlab="pattern",
           ylab="proportion",
           cex.lab=1.5,
           ylim=c(0, yLimit2))
      points(compareData$estimatedDistribution, pch=4, col='red')
      par(new=TRUE)
      plot(compareData$Coverage,
           col='blue',
           pch=3,
           xlab="",
           ylab="",
           ylim=c(0,sum(compareData$Coverage)*yLimit2),
           axes=FALSE)
      axis(side=4, at=c(compareData$Coverage, 0, 10))
      mtext("coverage", side=4, line=2.5)
      abline(0, 0, col="grey")
}

# vim:ft=r:ts=2:sw=2:sts=2:expandtab:
