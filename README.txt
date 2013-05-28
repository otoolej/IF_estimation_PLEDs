Methods to Estimate the Instantaneous Frequency (IF) of EEG Periodic
Discharges


Table of Contents
─────────────────

1 Overview
.. 1.1 files:
2 start
3 Requirements
4 Copyright
5 Version
6 Test computer setup
7 References
8 Contact


1 Overview
══════════

  This collection of M-files (computer code) implements a method to
  estimate the IF for periodic discharges on EEG [1].  Matlab or Octave
  (programming environments) are required.


1.1 files:
──────────

  estIF_epoch_overlapadd.m: main function to estimate IF with overlap and
                            add
  : ….. TO FINISH …


2 start
═══════

  Add this direction to the Matlab path:

  >> load_paths_all;

  Generate a figure from [2]:

  >> generate_Fig4;

  Estimate the IF using the proposed method (log-lag time-frequency
  filtering):

  >>   b=load('PLED_example_epoch.mat');
  >>   [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(b.x,b.Fs);

  >>   figure(1); clf; n=1:length(iflaw);
  >>   plot(n.*t_scale,iflaw.*f_scale); xlim([10 30]); ylim([0 5]);
  >>   xlabel('time (seconds)'); ylabel('frequency (Hz)');


3 Requirements
══════════════

  Either Matlab (R2012 or newer,
  [http://www.mathworks.co.uk/products/matlab/]) or Octave (v3.6 or
  newer, [http://www.gnu.org/software/octave/index.html], with the
  'octave-signal' add-on package).  In addition, the fast time-frequency
  distribution algorithms are required [2], and are available for free
  download at [https://spideroak.com/browse/share/fast_DTFDs/fastTFDs].


4 Copyright
═══════════

  Copyright (c) 2013 John M. O' Toole, University College Cork All
  rights reserved.  Email: j.otoole@ieee.org

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  ‣ Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  ‣ Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  ‣ Neither the name of the University College Cork nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY JOHN M. O' TOOLE ''AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JOHN M. O' TOOLE BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
  IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


5 Version
═════════

  ⁃ Version: 0.12
  ⁃ Last update [2013-05-28 Tue]


6 Test computer setup
═════════════════════

  • hardware: Intel(R) Xeon(R) CPU E5-1603 0 @ 2.80GHz; 8GB memory.
  • operating system: Ubuntu GNU/Linux x86_64 distribution (Raring,
    13.04), with Linux kernel 3.5.0-28-generic
  • software: Octave 3.6.4 (using Gnuplot 4.6 patchlevel 1), with
    'octave-signal' toolbox and Matlab (R2009b, R2012a, and R2013a)


7 References
════════════

  [2] B. Boashash, A. Ghazem, J.M. O' Toole, Time-frequency processing
  of non-stationary signals to aid diagnosis: highlights from medical
  applications, IEEE Signal Processing Magazine, 2013, in press


8 Contact
═════════

  • John M. O' Toole,
  • Neonatal Brain Research Group, Department of Paediatrics and Child
    Health, University College Dublin, Western Gateway Building, Room
    2.17, Cork, Ireland
  • Email j.otoole@ieee.org; jotoole@ucc.ie


