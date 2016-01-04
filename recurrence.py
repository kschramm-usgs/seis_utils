"""recurrence.py

This package contains methods for different recurrence and earthquake frequency
statistics.

Austin Holland
Oklahoma Geological Survey
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import *

def bin_mag(m,dmag=0.1,mc=0.0,cum=True):
  """bin_mag(dmag=0.1,mc=0.0,cum=True)
  dmag is the bin size
  mc is the minimum magnitude to consider 
  number of events are the number of events greater than the magnitude in bin
  If cum is False then the number of events are differenced.
  """
  
  bins=np.arange(mc,np.max(m)+dmag,dmag)
  hist=[]
  # This loop may be slow, but it works correctly - strange bug for m_max
  for b in bins[:-1]:
    hist.append(len(np.where(m>=b)[0]))
  hist.append(len(np.where(m>=np.max(m))[0])) # bug fix to catch the largest quake
  hist=np.array(hist)
  if cum==False:
    hist=hist[::-1]
    hist=np.diff(hist)
    hist=hist[::-1]
  #print len(bins),len(hist)
  return bins,hist,m
  
def b_MLE(m,uncert=True,mc=0.0):
  """b_value[,b_sigma]=b_MLE(m,uncert=True)
  Return the b-value given an numpy array of magnitudes using the MLE determination of 
  Aki (1965) if uncert=True then return the 1-sigma uncertainty in the estimate.
  """
  m=np.extract(m>=mc,m)
  #print mc,'<=',np.min(m)
  n=len(m)
  b=np.log10(np.e)/(np.mean(m)-np.min(m))
  if uncert:
    return b,0.98*b/np.sqrt(n)
  else:
    return b

def gr_lsqr(mags,nev,mc=0.0):
  """gr_lsqr(mags,nev,mc=0.0)
  Calculate the Gutenberg Richter values using least squares.
  """
  from scipy import polyfit, polyval
  bins=[]
  count=[]
  for i,mag in enumerate(mags):
    if mag>=mc:
      bins.append(mag)
      count.append(nev[i])
  (b,a,V)=polyfit(bins,np.log10(count),1,cov=True)
  print V
  
#   For whatever reason a and b are switched
  y=polyval([b,a],bins)
  fit=10**y
  return fit,bins,a,np.abs(b) 
  
def gr_mle(m,dmag=0.1,mc=0.0,method='Aki',scale=1.,cum=True):
  """gr_lsqr(mags,nev,mc=0.0)
  Calculate the Gutenberg Richter values using the maximum likelihood estimator.
  method can be Aki (1965) or Bender (1983) difference should be minimal
  method='Aki','Bender'
  """  
  #print mc
  bins,neq,m=bin_mag(np.array(m),dmag=dmag,mc=mc,cum=cum)
  if method=='Aki':
    neq=neq/scale
    b,b_std = b_MLE(m,uncert=True,mc=mc)
#    print "Aki (1965)",b,b_std
  elif method=='Bender':
    neq=neq/scale
    b,b_std = b_max_likelihood(bins,neq,dmag=dmag,m_c=mc)
#    print "Bender (1983)",b,b_std
  else:
    print "Unrecognized method"
    
  a,a_std=a_value(b,bins,neq,mc=mc)
  return b,b_std,a,a_std,len(m)
        
def a_value(b,mval,number_obs,mc=0.0):
  """ Calculation of a-value given a b-value and a histogram of events
  
  a_value(b,mval,number_obs,mc=0.0)
  """
  nlog=[]
  m=[]
  for i in np.arange(0,len(number_obs)):
    if (number_obs[i]>0.) and (mval[i]>=mc):
      nlog.append(np.log10(number_obs[i]))
      m.append(mval[i])
  nlog=np.array(nlog)
  m=np.array(m)
  a=nlog+b*m
  #print a
  return np.mean(a),np.std(a)

def b_max_likelihood(mval, number_obs, dmag=0.1, m_c=0.0):
    """
    Calculation of b-value and its uncertainty for a given catalogue,
    using the maximum likelihood method of Aki (1965), with a correction
    for discrete bin width (Bender, 1983).

    :param mval: array of reference magnitudes
                 (column 0 from recurrence table)
    :type mval: numpy.ndarray
    :param number_obs: number of observations in magnitude bin
                       (column 1 from recurrence table)
    :type number_obs: numpy.ndarray
    :keyword dmag: magnitude interval
    :type dmag: positive float
    :keyword m_c: completeness magnitude
    :type m_c: float
    :returns: bvalue and sigma_b
    :rtype: float
    """

    # Exclude data below Mc
    id0 = mval >= m_c
    mval = mval[id0]
    number_obs = number_obs[id0]
    # Get Number of events, minimum magnitude and mean magnitude
    neq = np.sum(number_obs)
    m_min = np.min(mval)
    m_ave = np.sum(mval * number_obs) / neq
    # Calculate b-value
    bval = np.log10(np.exp(1.0)) / (m_ave - m_min + (dmag / 2.))
    # Calculate sigma b from Bender estimator
    sigma_b = np.sum(number_obs * ((mval - m_ave) ** 2.0)) / (neq * (neq - 1))
    sigma_b = 2.3 * (bval ** 2.0) * np.sqrt(sigma_b)
    return bval, sigma_b  
    
def simple_probability(b,a,M_m):
  """ Remove the consideration of Omori aftershock decay from the probability of Reasenberg & Jones
  (1989,1990,1994) from Bachmann et al. (2011)
  """
  P=1-np.exp(-10**(a-b*M_m))
  return P
  
class TimeFrame():
  """ Class to support time frames fro recurrence statistics
  """
   
  def __init__(self,dt,desc=''):
    """ TimeFrame(dt,desc='')
    dt is a timedelta object in days
    desc is how one wants to refer to the time frame
    """
    self.dt=dt
    self.desc=desc
    
  def __str__(self):
    return "%s dt=%d days" % (self.desc,self.dt.total_seconds()/(24*60*60))
    
  def prob2html(self):
    str="<tr><td style=\"font-weight: bold\">%s</td>" % (self.desc)
    for p in self.P:
      str+="<td>%.4f</td>" % (p)
    str+="<\tr>\n"
    return str
    
  def gr2html(self):
    str="<tr><td style=\"font-weight: bold\">%s</td>" % (self.desc)
    str+="<td> %.4f &plusmn; %.4f</td>" % (self.b,self.b_std)
    str+="<td> %.4f &plusmn; %.4f</td>" % (self.a,self.a_std)
    str+="<td> %d </td>" % (self.neq)
    str+="<\tr>\n"
    return str

def probheader(M_m):
  str="<tr><th>&nbsp;</th><th colspan=\"%d\">Magnitude (m)</th></tr>\n" % (len(M_m))
  str+="<tr><th>Duration</th>"
  for m in M_m:
    str+="<th>%.1f</th>" % (m)
  str+="</tr>\n"
  return str  

def grheader():
  str="<tr><th>Duration</th><th>b-value</th><th>a-value</th><th>&#35; Earthquakes</tr>\n"
  return str
  
def daily_probabilities(dur,starttime,m_m,m_c=2.5,dt=timedelta(days=1),endtime=None):
  """prob_time(dt,starttime,m_m,m_c=2.5,endtime=None)
     return the daily probability for earthquakes 
     with magnitudes m_m=[] is a list of magnitudes to calculate the probability for
     dur is the duration to consider
  """
  day=timedelta(days=1)