# Brain-Heart-Interaction-indexes
The current repository contains an implementation of the Brain-Heart Interaction model proposed by Catrambone et al in 2019 ("Time-Resolved Directional Brain–Heart Interplay Measurement Through Synthetic Data Generation Models").

BHImodel2() is a function that quantify directional Brain-Heart Interplay (BHI) 
adopting the BHI model proposed by Catrambone et al.(2019) [1].

INPUT variables:
- TFR_EEG = Time course of EEG PSD, it must be a matrix: Channels X time points. It should be already filtered in the desired band of interest (psi)
- TFR_HRV = Time course of HRV PSD. It should be already filtered in the desired band of interest. (phi)
- FS      = Sampling Frequency of the two TFRs
- RR      = HRV time series (in seconds)
- win_RR  = windows length (in seconds) in which the RR model is reconstructed (default = 15 secs)
- window  = windows length (in seconds) in which the parameters are calculated (by default, it is assured that: window*FS >= 15 )

OUTPUT variables:
- HeartToBrain = Functional coupling index (c_rrTOeeg(T)) from the given Phi-band of HRV to the given Psi-band of the EEG
- BrainToHF, BrainToLF  = Functional coupling indexes from the given Psi-band of EEG to the LF or HF HRV bands
- HeartToBrain_sigma, HeartToBrain_mc = other parameters that can be extracted from the model in [1]

It is assumed that the signals given as input the BHImodel2() function are artifact free, from algotirhmic and physiological point of view (e.g. ectopic beats, eye blink, movement, ecc.)

---------------------------------------------------------------------------------------------

 This code is an implementation of the following scientific research, where an interested reader can find all the mathematical and teorethical details:
 [1] Catrambone Vincenzo, Alberto Greco, Nicola Vanello, Enzo Pasquale Scilingo, and Gaetano Valenza. "Time-Resolved Directional Brain/Heart Interplay Measurement Through Synthetic Data Generation Models." Annals of biomedical engineering 47, no. 6 (2019): 1479-1489.
 
---------------------------------------------------------------------------------------------

Copyright (C) 2019 Vincenzo Catrambone, Gaetano Valenza

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

If you use this program in support of published research, please include a citation of the reference above. If you use this code in a software package, please explicitly inform the end users of this copyright notice and ask them to cite the reference above in their published research.

---------------------------------------------------------------------------------------------

To use the software from Matlab or Octave, simply call the BHImodel function in the src/folder. Type 'help BHImodel' from the Matlab/Octave command window for help on the command's syntax and input/output arguments.

The software does not come with a GUI. Assuming 'x' and 'y' being the time course, sampled at 1 Hz, of Power spectral density of an EEG channel in the theta band (4-8 Hz) and HRV series respetively (integrated in LF+HF bands), and 'z' being the original HRV series, the following example performs the BHImodel analysis (with default parameters) and plots the results:

Fs = 1;

[HeartToTheta, ThetaToHF, ThetaToLF] = BHImodel(x,y,Fs,z);

figure, hold all

plot(1:length(HeartToTheta))'/Fs, HeartToTheta); 

plot(1:length(ThetaToHF))'/Fs, ThetaToHF);

plot(1:length(ThetaToHF))'/Fs, ThetaToLF);

---------------------------------------------------------------------------------------------
