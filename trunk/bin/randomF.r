#03-16-06 functions called by randomF.py


randomF_with_test <- function(data, correctList, mty=2, ntree=50)
{
	correctL = factor(correctList)
	ind <- sample(2, nrow(data), replace = TRUE, prob=c(0.8, 0.2))
	print (nrow(data))
	#data_rf = randomForest(data, y=correctL, ntree=50)
	data_rf = randomForest(data[ind==1,], y=correctL[ind==1], xtest=data[ind==2, ], ytest=correctL[ind==2], ntree=ntree, mty=mty)
	# pred = predict(data_rf, data[ind==2,])
	#, xtest=data[ind==2,c(1,2,3,4,5)], ytest=data[ind==2, 6])
	return (data_rf)
}


randomF <- function(data, correctList, mty=2, ntree=50)
{
	correctL = factor(correctList)
	print (nrow(data))
	data_rf = randomForest(data, y=correctL, ntree=ntree, mty=mty)
	# pred = predict(data_rf, data[ind==2,])
	#, xtest=data[ind==2,c(1,2,3,4,5)], ytest=data[ind==2, 6])
	return (data_rf)
}

predictRandomF <- function(rf_model, new_data)
{
	pred = predict(rf_model, new_data)
	return (pred)
}

spliceData <- function(data)
{
	ind <- sample(2, nrow(data), replace = TRUE, prob=c(0.8, 0.2))
	return (data[ind==2,])
}
