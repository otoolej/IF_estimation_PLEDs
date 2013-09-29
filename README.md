# Methods to Estimate the Instantaneous Frequency (IF) of EEG Periodic Discharges

## Overview 

This collection of M-files (computer code) implements a method to estimate the IF for periodic discharges on EEG [1].  Requires Matlab or Octave (programming environments) and the [fast_TFDs]((https://github.com/otoolej/fast_TFDs/zipball/master) package.


## quick start
Add this direction to the Matlab path:
  
{% highlight matlab %}

>> load_paths_all;

{% endhighlight  %}

Generate a figure from [2]:

{% highlight matlab %}
  >> generate_Fig4;
{% endhighlight %}

Estimate the IF using the proposed method (log-lag time-frequency filtering):

{% highlight matlab %}
  >>   b=load('PLED_example_epoch.mat');
  >>   [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(b.x,b.Fs);

  >>   figure(1); clf; 
  >>   n=1:length(iflaw);
  >>   plot(n.*t_scale,iflaw.*f_scale);
  >>   xlim([10 30]); ylim([0 5]);
  >>   xlabel('time (seconds)'); 
  >>   ylabel('frequency (Hz)'); 

{% endhighlight %}

## Requirements

Either Matlab (R2012 or newer, http://www.mathworks.co.uk/products/matlab/) or Octave
(v3.6 or newer, http://www.gnu.org/software/octave/index.html, with the 'octave-signal'
add-on package).  In addition, the fast time-frequency distribution algorithms are
required [[3(#ref_fasttfds)]], and are available for free [download](https://github.com/otoolej/fast_TFDs/zipball/master).



## Test computer setup
- hardware:  Intel(R) Xeon(R) CPU E5-1603 0 @ 2.80GHz; 8GB memory.
- operating system: Ubuntu GNU/Linux x86_64 distribution (Raring,
  13.04), with Linux kernel 3.5.0-28-generic
- software: Octave 3.6.4 (using Gnuplot 4.6 patchlevel 1), with 
  'octave-signal' toolbox and Matlab (R2009b, R2012a, and R2013a)

---

## References

1. in preparation  

2. <a id="ref_spmag"></a> B. Boashash, A. Ghazem, J.M. O' Toole, Time-frequency processing of non-stationary signals to aid diagnosis: highlights from medical applications, IEEE Signal Processing Magazine, 2013, in press

3. <a id="ref_acha"></a> J.M. O' Toole and B. Boashash, "Fast and memory-efficient algorithms for computing quadratic time--frequency distributions", Applied and Computational Harmonic Analysis, Feb. 2013, doi:10.1016/j.acha.2013.01.003

---

## Contact

- John M. O' Toole,  
- Neonatal Brain Research Group, 
  Department of Paediatrics and Child Health,
  University College Dublin,
  Western Gateway Building, Room 2.17,
  Cork, Ireland
- Email j.otoole@ieee.org; jotoole@ucc.ie	

