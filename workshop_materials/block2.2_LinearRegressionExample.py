""" LinearRegressionExample.py

Synthetic data example in 2 dimensions. """
import sys
sys.path.append("..\\")

## Standard stuff
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import methods
np.random.seed(23)

def axes_setup(axes):
	axes.spines["left"].set_position(("axes",-0.025))
	axes.spines["top"].set_visible(False)
	axes.spines["right"].set_visible(False)
	return

def LinearRegression(X,Y):

	""" Linear regression, this implementation computes
	the estimator and covariance matrix based on marginalized variance. 

	NB: X is assumed to be an array with shape = (num_data points, num_features) """

	## Get the dimensions
	num_data_points = X.shape[0]
	try:
		num_features = X.shape[1]
	except:
		num_features = 1
		X = X.reshape((num_data_points,1))

	## Compute the required matrix inversion
	## i.e. inv(x.T*x), which comes from minimizing
	## the residual sum of squares (RSS) and solving for
	## the optimum coefficients.
	xTx_inv = np.linalg.inv(np.dot(X.T,X))

	## Now use that matrix to compute the optimum coefficients
	## and their uncertainty.
	beta_hat = np.dot(xTx_inv,np.dot(X.T,Y))

	## Compute the estimated variance in the data points
	residual = Y - np.dot(X,beta_hat)
	RSS = (residual)**2
	var = RSS.sum(axis=0)/(num_data_points - 1)

	## Then the uncertainty (covariance matrix) is simply a 
	## reapplication of the inv(x.T*x):
	beta_var = var*xTx_inv

	## Reshape the outputs
	if num_features > 1:
		beta_hat = beta_hat.reshape((num_features,))
	else:
		beta_hat = beta_hat[0]
		beta_var = beta_var[0,0]

	return beta_hat, beta_var, residual

if __name__ == "__main__":

	## Make a simple dataset
	x = np.linspace(0,10,30)
	sigma = 2.
	y = 0.8*x + 7. 
	y = y + np.random.normal(0,sigma,size=x.shape) ## Measurement uncertainty

	## Fit the regression model
	X = np.ones((len(x),2))
	X[:,1] = x
	beta_hat, beta_var, _ = LinearRegression(X,y)

	## Propogate uncertainty to the fit
	y_hat = np.dot(X,beta_hat)
	y_cov = np.dot(X,np.dot(beta_var,X.T))
	y_samples = np.random.multivariate_normal(y_hat,y_cov,
											  size=(1000,))
	y_low = np.percentile(y_samples,2.5,axis=0)
	y_high = np.percentile(y_samples,97.5,axis=0)

	## Plot the result
	fig, axes = plt.subplots(figsize=(9,7))
	axes_setup(axes)
	axes.plot(x,y,
			  ls="None",
			  marker="o",markersize=10,
			  color="k",
			  zorder=5)
	axes.fill_between(x,y_low,y_high,
					  edgecolor="None",facecolor="grey",
					  alpha=0.4)
	axes.plot(x,y_hat,color="grey",lw=3)
	axes.set(xlabel=r"Input variable, $X$",
			 ylabel=r"Output variable, $Y$")
	fig.tight_layout()
	fig.savefig("_plots\\lr_example.png")
	plt.show()

