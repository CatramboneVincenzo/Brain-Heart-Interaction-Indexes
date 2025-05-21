# Brain-Heart Interaction Indexes
% This function quantifies directional Brain-Heart Interplay (BHI) through the model proposed by Catrambone et al.(2019) [1].

% INPUT variables:
% TFR_EEG = Time course of EEG power spectral density (PSD). This must be a matrix (Dimension: Channels X time)
%           samples. Each series in each row should be filtered in the desired frequency band of interest (psi)
% TFR_HRV = Time course of HRV PSD (Dimension: 1 X time). This should be filtered in the desired frequency band of interest (phi: e.g., LF or HF band)
% FS      = Sampling Frequency of the two TFRs
% RR      = HRV series (expressed in seconds)
% win_RR  = windows length (expressed in seconds) in which the heartbeat generation model (IPFM) is reconstructed (default = 15s)
% window  = windows length (in seconds) in which the parameters are calculated (default: window*FS >= 15 )
% TV: time-varying flag. 0 for punctual estimate, 1 for time-resolved estimate.

% OUTPUT variables:
% - HeartToBrain = Functional coupling index (c_rrTOeeg(T)) from HRV Phi-band to EEG Psi-band
% - BrainToHF, BrainToLF  = Functional coupling indices from  EEG Psi-band to  HRV-LF or  HRV-HF bands
% - HeartToBrain_sigma, HeartToBrain_mc = model parameters to be used for fitting evaluation [1]
% 
% This software assumes that input series are all artifact free, e.g., heartbeat dynamics free of algotirhmic and/or physiological artifacts; e.g.
% EEG series free of artifacts from eye blink, movement, etc.
% ---------------------------------------------------------------------------------------------
%  This code implements the theoretical dissertation published in:
%  [1] Catrambone Vincenzo, Alberto Greco, Nicola Vanello, Enzo Pasquale Scilingo,
%  and Gaetano Valenza. "Time-Resolved Directional Brain–Heart Interplay Measurement 
%  Through Synthetic Data Generation Models." Annals of biomedical engineering 47, no. 6 (2019): 1479-1489.
% ---------------------------------------------------------------------------------------------
% Copyright (C) 2019 Vincenzo Catrambone, Gaetano Valenza
% 
% This program is a free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
% 
% If you use this program in support of published research, please include a
% citation of the reference above. If you use this code in a software package,
% please explicitly inform the end users of this copyright notice and ask them
% to cite the reference above in their published research.
% ---------------------------------------------------------------------------------------------
% To use the software from Matlab or Octave, simply call the BHImodel function in
% the src/folder. Type 'help BHImodel' from the Matlab/Octave command window for
% help on the command's syntax and input/output arguments.
% 
% The software does not come with a GUI. 
% Assuming the series sampled at 1 Hz, with 'x' as the time-varying PSD of a given EEG channel  
% integrated in the theta band (4-8 Hz) and 'y' as the time-varying PSD of an HRV series integrated in 0.04-0.4 Hz,
% and 'z' as the HRV series, the following example performs the BHImodel 
% analysis (with default parameters) and plots the function outcomes:

---------------------------------------------------------------------------------------------
