import numpy as np
import logging

class Response():
  """ class Response 
  This class offers a number of helper functions that allow one to look at and evaluate 
  instrument response.  This includes converting from poles and zeros in Hz and converting
  to the SEED standard of radians/s.  It also has helper functions to compare different
  responses on the same plot and read dataless seed files and plot the associated response.
  """

  def __init__(self,desc='',poles=[],zeros=[],a0=1.,units='',gain=1.):
    """ __init__(desc='',poles=[],zeros=[],a0=1.,units='',gain=1.)
    Private method to build the response where desc is a description of the data
    poles is an list of poles, zeros is a list of zeros a0 is the normalization constant 
    gain is the seismometer gain and units can indicate (Hz, or radians/s)
    """
    self.desc=desc
    self.gain=gain
    self.poles=poles
    self.zeros=zeros
    self.a0=a0
    self.units=units

  def __str__(self):
    """ __str__()
    Private method to create a user readable output from a response object
    """
    ostr="Response %s - gain: %g\n" % (self.desc,self.gain)
    ostr+="Units %s\n" % (self.units)
    ostr+="Normalization Factor A0: %g\n" % (self.a0)
    ostr+="Poles (Number %d):\n" % (len(self.poles))
    for p in self.poles:
      ostr+="\t%.4e\t\t%.4e i\n" % (p.real,p.imag)
    ostr+="Zeros (Number %d):\n" % (len(self.zeros))
    for z in self.zeros:
        ostr+="\t%.4e\t\t%.4e i\n" % (z.real,z.imag)
    ostr+="\n"
    return ostr
    
  def plot(self,nfft=16**2,t_sample=0.005,show=False):
    """ plot(self,nfft=16**2,t_sample=0.005,show=False)
    Method to plot the response of an object where nfft and t_sample are documented in
    obspy.signal.pazToFreqResp.  Currenty all figures are generated in PDF format and saved
    by the unique description _freqresp.pdf.
    
    Requires:
    matplotlib
    obspy
    """
    import matplotlib.pyplot as plt
    from obspy.signal import pazToFreqResp
    h,f=pazToFreqResp(self.poles,self.zeros,self.a0,t_sample,nfft,freq=True)
    plt.figure()
    plt.loglog(f,np.abs(h))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.grid(which='both')
    plt.savefig(self.desc+'_freqresp.pdf',bbox='tight',transparent=True)
    if show==True:
      plt.show()

  def fromHz2radians(self):
    if self.units=='Hz':
      self.poles=2*np.pi*np.array(self.poles)
      self.zeros=2*np.pi*np.array(self.zeros)
      logging.debug("#poles - #zeros",len(self.poles)-len(self.zeros))
      self.a0=self.a0*np.power(2*np.pi,(len(self.poles)-len(self.zeros)))
      self.units='radians/s'
    else:
      logging.warn("Units must be in Hz")

  def fromPAZdict():
    pass