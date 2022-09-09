import numpy as np

def normalPDF(x,mu,sigma):
	p = 1.0/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2.0*sigma**2.0))
	return(p)

def mixture_logLikelihood(x,rho,mu_l,sigma_l,mu_h,sigma_h):
	pG_l = normalPDF(x,mu_l,sigma_l)
	pG_h = normalPDF(x,mu_h,sigma_h)
	logLik = sum(np.log(rho*pG_l + (1-rho)*pG_h))
	return(logLik)

def EM_gaussian_mixture_model(x,rho,mu_l,sigma_l,mu_h,sigma_h,eps=0.0001):
	logLik = []
	while True:
		# calculate likelihood of x
		pG_l = normalPDF(x,mu_l,sigma_l)
		pG_h = normalPDF(x,mu_h,sigma_h)
		logLik.append(sum(np.log(rho*pG_l + (1-rho)*pG_h)))
		# test for convergence
		if len(logLik)>1 and (logLik[-1] - logLik[-2]) < eps:
			break
		# calculate probabilities
		ppG_l = rho*pG_l/(rho*pG_l + (1-rho)*pG_h)
		ppG_h = (1-rho)*pG_h/(rho*pG_l + (1-rho)*pG_h)
		# update parameter
		rho = sum(ppG_l)/len(x)
		mu_l = sum(x*ppG_l)/sum(ppG_l)
		sigma_l = np.sqrt(sum(ppG_l*(x-mu_l)**2)/sum(ppG_l))
		mu_l = sum(x*ppG_l)/sum(ppG_l)
		sigma_l = np.sqrt(sum(ppG_l*(x-mu_l)**2)/sum(ppG_l))

	logLik = np.array(logLik)
	return(logLik,rho,mu_l,sigma_l,mu_h,sigma_h)

