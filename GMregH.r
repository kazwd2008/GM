###################################################################
###################################################################
# Original code is of function ÅgbmregÅh included in the WRS package 
# provided at Github [https://github.com/nicebread/WRS] by Dr. Felix
# Schonbrodt for Wilcox (2012) Introduction to Robust Estimation and
# Hypothesis Testing, 3rd ed., Elsevier. 
###################################################################
# bmreg: compute a bounded M regression using Huber Psi and Schweppe
# weights. The predictors are assumed to be stored in the n by p matrix # x.
  #----------------------------------------------------------------
# GMregH is prepared based on bmreg for weighted estimation K. Wada.
#----------------------------------------------------------------

GMregH <- function(x, y, g=rep(1, length(y)), itr=150, bend=2*sqrt((ncol(x)+1)/nrow(x))){

  x		<-as.matrix(x)
  init	<-lsfit(x,y,g)
  resid	<-init$residuals
  x1	<-cbind(1,x)
  nu	<-sqrt(1-hat(x1))
  low	<-ncol(x)+1

  for(it in 1:itr){
	ev    <- sort(abs(resid))
	scale <- median(ev[c(low:length(y))])/qnorm(.75)
	rov   <- (resid/scale)/nu
	psi   <- ifelse(abs(rov) <= bend,rov,bend*sign(rov))  # Huber Psi
	wt    <- nu*psi/(resid/scale)
	new   <- lsfit(x,y,wt*g)

	if(max(abs(new$coef-init$coef)) < .0001) break
		init$coef <- new$coef
		resid     <- new$residuals
  }

  resid<- y-x1%*%new$coef
  if (max(abs(new$coef-init$coef)) >= .0001)
	paste("failed to converge in",itr,"steps")
	list(coef=new$coef,resid=resid,w=wt, it=it)
}
###################################################################
