# This is modified from AUCRFcv to output the ROC curve data and to correct
# a bug in the calculation of the AUC

get_spec <- function(sens, model_sens, model_spec){
	max(model_spec[sens <= model_sens])
}


cross_validate <- function(x, nCV=10, M=10){

	switch(class(x), AUCRF = {
	  callRF <- x$call
	  data <- x$data
	  callRF$data <- as.name("newData")
	  yname <- as.character(eval(x$call$formula)[[2]])
	}, stop("x must be a AUCRF object."))

	cvAUC <- rep(0, M)

	sens <- seq(0,1,0.01)
	spec <- data.frame(sens=sens)

	for (m in 1:M) {

	  mpredict <- NULL
	  indPermuted <- matrix(c(sample(rownames(data)), rep(NA,
	      nCV - nrow(data)%%nCV)), ncol = nCV, byrow = TRUE)

	  for (k in 1:nCV) {
	    indTest <- indPermuted[, k]
	    indTest <- indTest[!is.na(indTest)]
	    indTrain <- rownames(data)[!(rownames(data) %in% indTest)]
	    newData <- data[indTrain, ]
	    kaucRF <- eval(callRF)
	    mpredict <- rbind(mpredict, predict(kaucRF$RFopt, newdata = data[indTest, ],
												type = "vote"))
	  }

	  mvotes <- data.frame(y = data[, yname], mpredict[rownames(data), ])
	  class(mvotes) <- c("votes", "data.frame")
	  colnames(mvotes) <- c("y", "0", "1")

		cv_roc <- roc(mvotes$y,mvotes[,"1"])

		spec[,as.character(m)] <- sapply(sens, get_spec, cv_roc$sensitivities, cv_roc$specificities)
		cvAUC[m] <- cv_roc$auc
	}

	mean_sens_spec <- data.frame(sens=spec$sens, spec=apply(spec[,-1], 1, mean))

	return(list(cv_est = mean(cvAUC), sens_spec = mean_sens_spec))
}
