""" seisutils.response 
Work with seismic instrument response.  Certainly not all inclusive at this point.  Should
add phase plots as well.

Austin Holland <holland.austin@gmail.com>
"""
import numpy as np
import logging

class Response():
  """ class Response 
  This class offers a number of helper functions that allow one to look at and evaluate 
  instrument response.  This includes converting from poles and zeros in Hz and converting
  to the SEED standard of radians/s.  It also has helper functions to compare different
  responses on the same plot and read dataless seed files and plot the associated response.
  """

  def __init__(self,desc='',poles=[],zeros=[],a0=1.,units='',gain=1.,countsperV=1.):
    """ __init__(desc='',poles=[],zeros=[],a0=1.,units='',gain=1.,,countsperV=1.)
    Private method to build the response where desc is a description of the data
    poles is an list of poles, zeros is a list of zeros a0 is the normalization constant 
    gain is the seismometer gain and units can indicate (Hz, or radians/s)
    """
    self.desc=desc
    self.gain=gain
    self.poles=poles
    self.zeros=zeros
    self.countsperV=countsperV
    self.a0=a0
    self.sensitivity=self.gain*self.countsperV
    self.units=units

  def __str__(self):
    """ __str__()
    Private method to create a user readable output from a response object
    """
    ostr="Response %s - gain: %g\n" % (self.desc,self.gain)
    ostr+="Digitizer counts/V %g\n" % (self.countsperV)
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
    
  def plot(self,nfft=16**2,t_sample=0.005,show=False,include_seismometer=True):
    """ plot(self,nfft=16**2,t_sample=0.005,show=False)
    Method to plot the response of an object where nfft and t_sample are documented in
    obspy.signal.pazToFreqResp.  Currenty all figures are generated in PDF format and saved
    by the unique description _freqresp.pdf.
    
    Requires:
    matplotlib
    obspy.signal.pazToFreqResp
    """
    import matplotlib.pyplot as plt
    from obspy.signal import pazToFreqResp
    h,f=pazToFreqResp(self.poles,self.zeros,self.a0,t_sample,nfft,freq=True)
    plt.figure()
    if include_seismometer:
      plt.loglog(f,np.abs(h)*self.gain)
      plt.ylabel('V/m/s')
    else:
      plt.loglog(f,np.abs(h)*self.gain)
      plt.ylabel('Amplitude')   
    plt.xlabel('Frequency (Hz)')    
    plt.grid(which='both')
    plt.savefig(self.desc+'_freqresp.pdf',bbox='tight',transparent=True)
    if show==True:
      plt.show()

  def fromHz2radians(self):
    """ fromHz2radians()
    Method to convert a response in Hz to one in radians per second.  This method
    requires the response units to be 'Hz'
    """
    if self.units=='Hz':
      self.poles=2*np.pi*np.array(self.poles)
      self.zeros=2*np.pi*np.array(self.zeros)
      logging.debug("#poles - #zeros",len(self.poles)-len(self.zeros))
      self.a0=self.a0*np.power(2*np.pi,(len(self.poles)-len(self.zeros)))
      self.units='radians/s'
    else:
      logging.warn("Units must be in Hz")

  def fromPAZdict(paz,desc='',units='radians/s'):
    """ fromPAZdict(paz,desc='',units='radians/s')
    Method to populate a response object from the PAZ dict output from obspy.xseed.Parser.
    Description (desc) is recomended as a good identifier.
    """
    self.gain=paz['seismometer_gain']
    self.poles=paz['poles']
    self.zeros=paz['zeros']
    self.a0=paz['gain']
    self.sensitivity=paz['sensitivity']
    self.units=units
    self.countsperV=self.sensitivity/self.gain
    
  def check_normalization(freq=1.0):
    """ check_normalization(freq=1.0)
    This method checks the normalization at or near a given frequency. This method 
    returns the normalization factor and the frequency at which the discrete determination 
    was evaluated.

    Requires:
    obspy.signal.pazToFreqResp
    """
    from obspy.signal import pazToFreqResp
    h,f=pazToFreqResp(self.poles,self.zeros,self.a0,t_sample,nfft,freq=True)
    # Find the index of the frequency closest to freq
    i_f=np.argmin(f-freq)
    norm=np.abs(h[i])
    return norm,f[i_f]
  
  def to_delimited(self,delimeter=",",freq=1.0):
    """ to_csv(self,delimeter=",",freq=1.0)
    This method creates a one line string for parameters delimited by the designated 
    delimeter the fields include:
    desc
    gain
    norm #Normalization factor at or near freq see check_normalization
    f_norm #Frequency for normalization factor
    countsperV
    sensitivity
    """
    norm,f_norm=self.check_normalization(freq=freq)
    out=[self.desc,self.gain,norm,f_norm,self.countsperV,self.sensitivity]
    if isinstance(delimeter, str):
      o_str=delimeter.join(map(str,out))
    else:
      logging.warn("Response.to_deimited warning: delimeter must be a string")
      o_str=None
    return o_str