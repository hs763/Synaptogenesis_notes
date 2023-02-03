"""
Additional functions for IB MCB lectures

This module contains functions used for IB MCB lectures on mixed models. They mostly focus on
diagnostic plotting.

The functions are based on a script originally written by Jason Sadowski for OLS data
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
# import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy as np
import scipy.stats as stats
import statistics
# imports to check for object type
from pymer4.models import Lmer
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.regression.linear_model import OLS

def _resid_fitted(self, ax = None, col = None):
  """Residual vs Fitted plot

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots from another function.
  col: string
    String with the name of the random effect

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """
  # get the fitted values
  fits = _get_fits(self)
  # extract the residuals
  residuals = _get_resids(self)
  # and find largest 3 residuals
  top3 = abs(residuals).sort_values(ascending=False).iloc[:3]
  # the smooth to overlay on the plot
  smoothed = lowess(residuals, fits)
  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()
  if col is not None:
    col = _get_data(self)[col]
  # now create a scatter plot
  sns.scatterplot(x=fits, y=residuals,hue= col, legend = None, ax=ax)
  # add the smooth
  ax.plot(smoothed[:, 0], smoothed[:, 1], color='r')
  # axes and labels
  ax.set_ylabel('Residuals')
  ax.set_xlabel('Fitted Values')
  ax.set_title('Residuals vs. Fitted')
  # add the dashed line at 0
  ax.plot([min(fits), max(fits)], [0, 0], color='k', linestyle=':')
  # label the outliers
  for i in top3.index:
    ax.annotate(i, xy=(fits[i], residuals[i]))
  return (ax)
    
def _qqplot(self, ax = None, col = None):
  """qqplot of residuals

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots from another function.
  col: string
    String with the name of the random effect

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """

  # make a copy of the data dataframe
  plot_data = _get_data(self)
  # add a standard residuals column
  plot_data['residuals'] = _get_resids(self)
  plot_data['standard_residuals'] = plot_data['residuals'].div(statistics.stdev(plot_data['residuals']))
  # if there is a col variable, add it to the dataframe
  if col is not None:
    plot_data[col] = plot_data[col].astype('category')
  # sort by the residuals
  plot_data.sort_values(by=['residuals'], inplace=True)
  # add theortical quantiles
  plot_data['theoretical_quantiles'] = stats.probplot(plot_data['standard_residuals'], \
                                              dist = 'norm', fit = False)[0]
  # TODO, it would be better to rank based on the resid of the standard on theoretical quantiles.
  # rather than directly on the standard residuals
  rankings = abs(plot_data['standard_residuals']).sort_values(ascending = False)
  # get the top three values
  top3 = rankings.iloc[:3]
  if isinstance(ax,type(None)):
    fig, ax = plt.subplots()

  x = plot_data['theoretical_quantiles']
  y = plot_data['standard_residuals']
  ax.plot([np.min([x,y]),np.max([x,y])],[np.min([x,y]),np.max([x,y])], \
           color = 'r', ls = '--')
  # we want a 1:1 line
#  min_xy = np.min([x,y])
#  max_xy = np.max([x,y])
#  ax.plot([min_xy,max_xy],[min_xy,max_xy], \
#           color = 'r', ls = '--')
  # now create a scatter plot
  sns.scatterplot(data=plot_data, x='theoretical_quantiles', y="standard_residuals",hue= col, legend = None, ax=ax)
  # axes and labels
  ax.set_title('Normal Q-Q')
  ax.set_ylabel('Standardized Residuals')
  ax.set_xlabel('Theoretical Quantiles')
  # label the outliers
  for val in top3.index:
    ax.annotate(val,xy=(plot_data['theoretical_quantiles'].loc[val], \
                        plot_data['standard_residuals'].loc[val]))

  return(ax)
   
def _scale_location(self,  ax = None, col = None):
  """scale location plot

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots from another function.
  col: string
    String with the name of the random effect

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """

  # get the residuals, standardise them and sqrt
  residuals = _get_resids(self)
  standard_residuals = residuals / statistics.stdev(residuals)
  sqrt_standard_resid = pd.Series(np.sqrt(np.abs(standard_residuals)))
  # and find largest 3 residuals
  sqrt_standard_resid.index = sqrt_standard_resid.index
  top3 = abs(pd.Series(sqrt_standard_resid)).sort_values(ascending=False).iloc[:3]
  # get the fitted values
  fits = _get_fits(self)
  # the smooth to overlay on the plot
  smoothed = lowess(sqrt_standard_resid, fits)

  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()
  # make a copy of the data dataframe
  plot_data = _get_data(self)
  # now create a scatter plot
  sns.scatterplot(data=plot_data, x=fits, y=sqrt_standard_resid, hue=col, legend=None, ax=ax)
  # add the smooth
  ax.plot(smoothed[:, 0], smoothed[:, 1], color='r')
  # axes and labels
  ax.set_ylabel('$\sqrt{|standardized \ Residuals|}$')
  ax.set_xlabel('Fitted Values')
  ax.set_title('Scale-Location')
  # label the outliers
  for i in top3.index:
    ax.annotate(i, xy=(fits[i], sqrt_standard_resid[i]))


def _cooks_plot(self, ax=None, col=None):
  """plot of Cooks distances

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots
  col: string
    String with the name of the random effect

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """

  # make a copy of the data dataframe
  plot_data = _get_data(self)
  # get the cook distance
  cooks_dist = _get_cooks(self)
  # get the top 3 values
  top3 = abs(pd.Series(cooks_dist)).sort_values(ascending=False).iloc[:3]
  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()
  # now create a scatter plot
  sns.scatterplot(data = plot_data, x = plot_data.index, y = cooks_dist, hue = col, legend = None, ax = ax)
  # add horizontal lines at 0.5 and 1
#  plt.axhline(y=1,color="r",linestyle="-")
 # plt.axhline(y=0.5,color="r",linestyle=":")
  # add horizontal line at 4N
  plt.axhline(y=4/len(plot_data),color="r",linestyle="-")
  # axes and labels
  ax.set_ylabel("Cook's distance")
  ax.set_xlabel('Index')
  ax.set_title("Cook's distances")
  # set the y limits of the plot to the range of the data plus 5% (or they will be extended to at least 1
  ax.set_ylim(0, 1.05*max(cooks_dist))
  # label the outliers
  for i in top3.index:
    ax.annotate(i, xy=(plot_data.index[i], cooks_dist[i]))


def _lmer_cooks_plot(self, ax=None, col=None):
  """plot of Cooks distances

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots
  col: string
    String with the name of the random effect

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """

  # make a copy of the data dataframe
  plot_data = self.data.copy()
  # get the cook distance
  cooks_dist = self.cooks_distance()
  # get the top 3 values
  top3 = abs(pd.Series(cooks_dist)).sort_values(ascending=False).iloc[:3]
  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()
  # now create a scatter plot
  sns.scatterplot(data = plot_data, x = plot_data.index, y = cooks_dist, hue = col, legend = None, ax = ax)
  # add horizontal lines at 0.5 and 1
#  plt.axhline(y=1,color="r",linestyle="-")
 # plt.axhline(y=0.5,color="r",linestyle=":")
  # add horizontal line at 4N
  plt.axhline(y=4/len(self.data),color="r",linestyle="-")
  # axes and labels
  ax.set_ylabel("Cook's distance")
  ax.set_xlabel('Index')
  ax.set_title("Cook's distances")
  # set the y limits of the plot to the range of the data plus 5% (or they will be extended to at least 1
  ax.set_ylim(0, 1.05*max(cooks_dist))
  # label the outliers
  for i in top3.index:
    ax.annotate(i, xy=(plot_data.index[i], cooks_dist[i]))

def lmer_plot_diag_resids(self, col = None):
  """diagnostic plots for Lmer models

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model from the pymer4 module
  
  Returns
  ---------------------------------------------------------
  axs: axis object
    A matplotlib axis object with 4 subplots. The order of the subplots is
       as follows:
       [0,0]: Residuals vs. Fitted values
       [0,1]: Normal QQ Plot
       [1,0]: Studentized Residuals vs. Fitted Values
       [1,1]: A cooks distance plot

  """

  # check object type
  if not (isinstance(self, Lmer)):
    raise Exception("this function is designed for Lmer models from pymer4")
  # check that if we have a variable to colour by, it is one of the random factors
  if col in self.ranef_var.index.values is False:
    raise Exception("error, the col variable is not a random effect")
  if self.family == 'gaussian':
    # set up the subplots
    fig, axs = plt.subplots(2,2)
    # and now create the plots
    _resid_fitted(self,ax = axs[0,0], col = col)
    _qqplot(self,ax = axs[0,1], col = col)
    _scale_location(self,ax = axs[1,0], col = col)
    _cooks_plot(self,ax = axs[1,1], col = col)
  else:
    # set up the subplots
    fig, axs = plt.subplots(1,2)
    # and now create the plots
    _resid_fitted(self,ax = axs[0], col = col)
    _cooks_plot(self,ax = axs[1], col = col)
  # sort out the layout
  fig.tight_layout()
  return(axs)



def lm_plot_diag_resids(self, col=None, R_leverage= False):
  """diagnostic plots for linear models (ols)

  Parameters
  ----------
  self: RegressionResultsWrapper
    A fitted OLS model from the statsmodels module

  Returns
  ---------------------------------------------------------
  axs: axis object
    A matplotlib axis object with 4 subplots. The order of the subplots is
       as follows:
       [0,0]: Residuals vs. Fitted values
       [0,1]: Normal QQ Plot
       [1,0]: Studentized Residuals vs. Fitted Values
       [1,1]: A leverage plot with contours of Cooks distances

  """
  # check object type
  if not (isinstance(self.model, OLS) or isinstance(self.model, GLM)):
    raise Exception("this function is designed for OLS models from statsmodels")
  # check that if we have a variable to colour by, it is one of the random factors
#  if col in self.ranef_var.index.values is False:
#    raise Exception("error, the col variable is not a random effect")
  # set up the subplots
  fig, axs = plt.subplots(2, 2)
  # and now create the plots
  _resid_fitted(self, ax=axs[0, 0], col=col)
  _qqplot(self, ax=axs[0, 1], col=col)
  _scale_location(self, ax=axs[1, 0], col=col)
  if (R_leverage):
    _lm_leverage(self, ax=axs[1, 1], col=col)
  else:
    _cooks_plot(self, ax=axs[1, 1], col=col)
  # sort out the layout
  fig.tight_layout()
  return (axs)



def glm_plot_diag_resids(self, col=None, R_leverage= False):
  """diagnostic plots for generalised linear models (ols)

  Parameters
  ----------
  self: GLMResultsWrapper
    A fitted GLM model from the statsmodels module

  Returns
  ---------------------------------------------------------
  axs: axis object
    A matplotlib axis object with 4 subplots. The order of the subplots is
       as follows:
       [0,0]: Residuals vs. Fitted values
       [0,1]: Normal QQ Plot
       [1,0]: Studentized Residuals vs. Fitted Values
       [1,1]: A leverage plot with contours of Cooks distances

  """
  # check object type
  if not (isinstance(self.model, GLM)):
    raise Exception("this function is designed for OLS models from statsmodels")
  # set up the subplots
  fig, axs = plt.subplots(1, 2)
  # and now create the plots
  _resid_fitted(self, ax=axs[0], col=col)
  _cooks_plot(self, ax=axs[1], col=col)
  # sort out the layout
  fig.tight_layout()
  return (axs)


#from pymer4.models import Lmer
#Lmer.plot_diag_resids = lmer_plot_diag_resids

def _get_fits(self):
  if (isinstance(self, Lmer)):
    if self.family == 'gaussian':
      return(self.fits)
    else:
      return(self.predict(data=self.data, verify_predictions=False,pred_type="link"))
  if isinstance(self.model, OLS) or isinstance(self.model, GLM):
    return(self.fittedvalues)
#  if type(self).__name__=="Lmer":
  raise Exception("error, unsupported model type")

def _get_resids(self):
  if (isinstance(self, Lmer)):
    return(pd.Series(self.residuals))
  if isinstance(self.model, OLS):
    return(self.resid)
  if isinstance(self.model, GLM):
    return(self.resid_deviance)
  raise Exception("error, unsupported model type")

def _get_data(self):
  if (isinstance(self, Lmer)):
    return (self.data.copy())
  if isinstance(self.model, OLS) or isinstance(self.model, GLM):
    return (self.model.data.frame.copy())
  raise Exception("error, unsupported model type")

def _get_cooks(self):
  if (isinstance(self, Lmer)):
    return (self.cooks_distance())
  if isinstance(self.model, OLS) or isinstance(self.model, GLM):
    influence = self.get_influence()
    cooks = influence.cooks_distance
    return (cooks[0])
  raise Exception("error, unsupported model type")


def _lm_leverage(fitted_model, col=None, student_residuals=None, \
             leverage=None, ax=None):
  """
  Parameters
  ---------------------------------------------------------
  fitted_model: A fitted linear regression model from the statsmodels package.
                Class: <statsmodels.regression.linear_model.OLS>
  student_residuals: A pandas series of the internally studentized residuals.
  ax: A specific matplotlib axis. Used if creating subplots

  Returns
  ---------------------------------------------------------
  ax: A matplotlib axis object

  The approach for coding the Cook's D lines comes from:
  https://emredjan.github.io/blog/2017/07/11/emulating-r-plots-in-python/

  By: Jason Sadowski
  Date: 2019-11-19
  """
  if isinstance(student_residuals, type(None)):
    student_residuals = pd.Series(fitted_model \
                                  .get_influence().resid_studentized_internal)
    student_residuals.index = fitted_model.resid.index
  if isinstance(leverage, type(None)):
    leverage = fitted_model.get_influence().hat_matrix_diag
  df = pd.DataFrame(student_residuals)
  df.columns = ['student_residuals']
  df['leverage'] = leverage
  sorted_student_resid = abs(df['student_residuals']) \
    .sort_values(ascending=False)
  top3 = sorted_student_resid.iloc[:3]
  smoothed = lowess(student_residuals, leverage)
  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()
  x = df['leverage']
  y = df['student_residuals']
  plot_data = _get_data(fitted_model)
  xpos = max(x) + max(x) * 0.05
  ax.scatter(x, y, edgecolors='k', facecolors='none')
  # now create a scatter plot
  sns.scatterplot(data = plot_data, x =x, y = y, hue = col, legend = None, ax = ax)

  ax.plot(smoothed[:, 0], smoothed[:, 1], color='r')
  ax.set_ylabel('Studentized Residuals')
  ax.set_xlabel('Leverage')
  ax.set_title('Residuals vs. Leverage')
  ax.set_ylim(min(y) - min(y) * 0.15, max(y) + max(y) * 0.15)
  ax.set_xlim(min(x) - max(x) * 0.05, max(x) + max(x) * 0.05)

  cooksx = np.linspace(min(x), xpos, 50)
  p = len(fitted_model.params)
  poscooks1y = np.sqrt((p * (1 - cooksx)) / cooksx)
  poscooks05y = np.sqrt(0.5 * (p * (1 - cooksx)) / cooksx)
  negcooks1y = -np.sqrt((p * (1 - cooksx)) / cooksx)
  negcooks05y = -np.sqrt(0.5 * (p * (1 - cooksx)) / cooksx)

  ax.plot(cooksx, poscooks1y, label="Cook's Distance", ls=':', color='r')
  ax.plot(cooksx, poscooks05y, ls=':', color='r')
  ax.plot(cooksx, negcooks1y, ls=':', color='r')
  ax.plot(cooksx, negcooks05y, ls=':', color='r')
  ax.plot([0, 0], ax.get_ylim(), ls=":", color='k', alpha=0.3)
  ax.plot(ax.get_xlim(), [0, 0], ls=":", color='k', alpha=0.3)
  ax.annotate('1.0', xy=(xpos, poscooks1y[-1]), color='r')
  ax.annotate('0.5', xy=(xpos, poscooks05y[-1]), color='r')
  ax.annotate('1.0', xy=(xpos, negcooks1y[-1]), color='r')
  ax.annotate('0.5', xy=(xpos, negcooks05y[-1]), color='r')
  ax.legend()
  for val in top3.index:
    ax.annotate(val, xy=(x.loc[val], y.loc[val]))
  return (ax)



def lmer_plot_diag_ranef(self):
  """diagnostic plots of randeom effects for Lmer models

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model from the pymer4 module

  Returns
  ---------------------------------------------------------
  axs: axis object
    A matplotlib axis object with 4 subplots. The order of the subplots is
       as follows:
       [0,0]: Residuals vs. Fitted values
       [0,1]: Normal QQ Plot
       [1,0]: Studentized Residuals vs. Fitted Values
       [1,1]: A cooks distance plot

  """
  # check object type
  if not (isinstance(self, Lmer)):
    raise Exception("this function is designed for Lmer models from pymer4")
  # coerce ranef to be a list (pymer4 returns a list one when multiple grouping factors are present
  if isinstance(self.ranef, list):
    ranef_list = self.ranef
  else:
    ranef_list = list()
    ranef_list.append(self.ranef)
  n_rows = self.ranef_var.shape[0]
  if self.family == "gaussian":
    n_rows = n_rows- 1
  groups = list(self.grps.keys())
  n_groups = len(groups)

  fig, axs = plt.subplots(n_rows, 2)

  row_counter = -1  # keep track of which row we are filling in
  for i_grp in range(n_groups):
    ranef_grp = ranef_list[i_grp].copy()  # the random effect dataframe for this group
    n_ranef_grp = len(ranef_grp.columns)  # number of columns in the dataframe
    for i_col in range(n_ranef_grp):
      row_counter = row_counter + 1
      id = ranef_grp.columns[i_col]
      this_ranef_pd = ranef_grp.iloc[:, 0]
      if n_rows > 1:
        _hist_ranef(ranef_grp, id=id, grp=groups[i_grp], ax=axs[row_counter, 0])
        _qqplot_ranef(ranef_grp, id=id, grp=groups[i_grp], ax=axs[row_counter, 1])
      else:
        _hist_ranef(ranef_grp, id=id, grp=groups[i_grp], ax=axs[0])
        _qqplot_ranef(ranef_grp, id=id, grp=groups[i_grp], ax=axs[1])
  # sort out the layout
  fig.tight_layout()
  return (axs)


def _qqplot_ranef(plot_data, id, grp, ax=None):
  """qqplot of random effects

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  id: string
    the id for the random effect
  grp: string
    the name of the grouping (i.e.random) factor
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots from another function.

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """

  # get the random effect
  # plot_data = self.ranef[id].to_frame()
  # add a standard residuals column
  plot_data['standard_coeffs'] = plot_data[id].div(statistics.stdev(plot_data[id]))
  # sort by the residuals
  plot_data.sort_values(by=['standard_coeffs'], inplace=True)
  # add theortical quantiles
  plot_data['theoretical_quantiles'] = stats.probplot(plot_data['standard_coeffs'], \
                                                      dist='norm', fit=False)[0]
  # TODO, it would be better to rank based on the resid of the standard on theoretical quantiles.
  # rather than directly on the standard residuals
  rankings = abs(plot_data['standard_coeffs']).sort_values(ascending=False)
  # get the top three values
  top3 = rankings.iloc[:3]
  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()

  x = plot_data['theoretical_quantiles']
  y = plot_data['standard_coeffs']
  ax.plot([np.min([x, y]), np.max([x, y])], [np.min([x, y]), np.max([x, y])], \
          color='r', ls='--')
  # we want a 1:1 line
  #  min_xy = np.min([x,y])
  #  max_xy = np.max([x,y])
  #  ax.plot([min_xy,max_xy],[min_xy,max_xy], \
  #           color = 'r', ls = '--')
  # now create a scatter plot
  sns.scatterplot(data=plot_data, x='theoretical_quantiles', y='standard_coeffs', legend=None, ax=ax)
  # axes and labels
  ax.set_title('Normal Q-Q')
  ax.set_ylabel('Std Random effect')
  ax.set_xlabel('Theoretical Quantiles')
  # label the outliers
  for val in top3.index:
    ax.annotate(val, xy=(plot_data['theoretical_quantiles'].loc[val], \
                         plot_data['standard_coeffs'].loc[val]))
  ax.set_aspect('equal', 'box')
  return (ax)



def _hist_ranef(plot_data, id, grp, ax=None):
  """histogram of random effects

  Parameters
  ----------
  self: Lmer
    A fitted Lmer model
  plot_data: a pandas data frame
    The coefficients for the random effect of interest
  id: string
    the id for the random effect
  grp: string
    the name of the grouping (i.e.random) factor
  ax: matplotlib.axis
    A specific matplotlib axis. Used if creating subplots from another function.

  Returns
  ----------
  ax: matplotlib.axis
    A matplotlib axis object
  """
  # get the random effect
  # plot_data = self.ranef[id].to_frame()
  if isinstance(ax, type(None)):
    fig, ax = plt.subplots()
  sns.histplot(data=plot_data, x=id, ax=ax)
  # axes and labels
  plot_label = id +" for "+ grp
  ax.set_title(plot_label)
  ax.set_xlabel('coefficients')
  ax.set_aspect('auto')
  return (ax)