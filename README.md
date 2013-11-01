# Methods to Estimate the Instantaneous Frequency (IF) of EEG Periodic Discharges


This collection of M-files (computer code) implements a method to estimate the IF for periodic discharges on EEG [[1]](#references).  Requires Matlab or Octave (programming environments) and the [fast_TFDs]((https://github.com/otoolej/fast_TFDs/zipball/master) package.


# Contents
- [quick start](#quick-start)
- [requirements](#requirements)
- [test computer setup](#test-computer-setup)
- [licence](#licence)
- [references](#references)
- [contact](#contact)


# quick start
The first and necessary step is to add paths so Matlab can find all the code; use the `load_paths_all` function:
```matlab
  >> load_paths_all;
```


### files
All Matlab files (.m files) have a description and an example in the header. To read this
header, type `help <filename.m>` in Matlab.  Directory structure is as follows: 
```
.
├── data                 # some signal examples in .mat files
├── methods              # IF estimation methods
├── synth_PLED_signal    # generate the synthetic (PLED-like) signals
│   └── generate_data
└── utils                # miscellaneous collection of utilities 
```


### example 1
An example, from [[1]](#references), compares the proposed IF estimation method with other methods:
```matlab
  >> demo_compare_3methods;
```

### example 2
To generate the PLED (periodic lateralized epilepiform discharges) example figure from [[1]](#references):
```matlab
  >> generate_Fig_IFest_exampl('a');
  >> generate_Fig_IFest_exampl('b');
  >> generate_Fig_IFest_exampl('c');
  >> generate_Fig_IFest_exampl('d');
  >> generate_Fig_IFest_exampl('e');
```

## synthetic signals

Call function `gen_synth_signals_plus_noise` (in the `synth_PLED_signals` directory) to
generate the synthetic signals:
```matlab
   snr=0; Fs=50; 
   N_iters=10; WAVEFORM_MONOPHASIC=0;
   [sigs,IFs]=gen_synth_signals_plus_noise(snr,N_iters,Fs,'CGN',WAVEFORM_MONOPHASIC);

   figure(1); clf; Fs=50; 
   n=(1:size(sigs,2))./Fs; 
   subplot(2,1,1); plot(n,sigs(1,:)); xlim([5 35]);
   xlabel('time (seconds)');  ylabel('amplitude');
   subplot(2,1,2); plot(n,IFs(1,:)); xlim([5 35]);
   xlabel('time (seconds)');  ylabel('frequency (Hz)');
```
 
These synthetic signals use the Duffing waveform from [[4]](#references), stored in data
file `data/duffing_waveforms.mat` directory.

### example 3
Estimate the IF on synthetic signal using the proposed method (log-lag time-frequency filtering):
```matlab
   b=load('synth_signal_example_0dB.mat');
   [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(b.x,b.Fs);
   
   % plot:
   figure(1); clf;  hold all;
   n=1:length(iflaw);
   plot(n.*t_scale,iflaw.*f_scale);
   plot( (1:length(b.true_IF))./b.Fs,b.true_IF);
   legend('proposed','true IF');
   xlim([10 30]); ylim([0 2]);
   xlabel('time (seconds)'); 
   ylabel('frequency (Hz)');
```



# requirements
Either Matlab (R2012 or newer, http://www.mathworks.co.uk/products/matlab/) or Octave (v3.6 or newer, http://www.gnu.org/software/octave/index.html, with the 'octave-signal' add-on package).  In addition, the fast time-frequency distribution algorithms are required [[3]](#references), and are available for free [download](https://github.com/otoolej/fast_TFDs/zipball/master).



# test computer setup
- hardware:  Intel(R) Xeon(R) CPU E5-1603 0 @ 2.80GHz; 8GB memory.
- operating system: Ubuntu GNU/Linux x86_64 distribution (Raring, 13.04), with Linux kernel 3.5.0-28-generic 
- software: Octave 3.6.4 (using Gnuplot 4.6 patchlevel 1), with 'octave-signal' toolbox and Matlab (R2009b, R2012a, and R2013a)

---

# licence

```
Copyright (c) 2011 and 2012, John O' Toole, University of Deusto
Copyright (c) 2013, John O' Toole, University College Cork
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

  Neither the name of the University of Deusto nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```


# references

1. in preparation  

2. B. Boashash, A. Ghazem, and J.M. O' Toole. "Time-Frequency Processing of Nonstationary Signals: Advanced TFD Design to Aid Diagnosis with Highlights from Medical Applications." IEEE Signal Processing Magazine, vol. 30, no. 6, pp. 108--119, November 2013.

3. J.M. O' Toole and B. Boashash, "Fast and memory-efficient algorithms for computing quadratic time--frequency distributions", Applied and Computational Harmonic Analysis, vol. 35, no. 2, pp. 350--358, 2013.

4. N.J. Stevenson, M. Mesbah, G.B. Boylan, P.B. Colditz, and B. Boashash, "A nonlinear model of newborn EEG with non-stationary inputs," Annals of Biomedical Engineering, vol. 38, no. 9, pp. 3010--3021, Sep. 2010.

---

# contact

John M. O' Toole

Neonatal Brain Research Group,  
Irish Centre for Fetal and Neonatal Translational Research,  
Department of Paediatrics and Child Health,  
University College Dublin,  
Western Gateway Building, Room 2.17,  
Cork, Ireland


