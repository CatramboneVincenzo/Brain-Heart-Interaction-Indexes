function [] = SDGM_LFHF(EEG,f_lims,RR,RRi,t_RRi,FS_rri,FS_bhi,TV,file_output)

%% Sinthetic Data Generation model implementation for Brain-Heart Interplay (BHI) time-resolved estimation
% Inputs:  
% EEG: EEGLAB-like structure with the preprocessed EEG signals
% f_lims: a matrix (Nb x 2) identifying the Nb frequency bands limits to be used
% in EEG analysis. Exaxmple: (to use theta, alpha, and beta bands) f_lims = [4 8; 8 12; 12 30];
% RR: HRV series
% RRi: HRV series interpolated
% t_RRi: time vectors corrisponding to RRi samples
% FS_rri: sampling frequency of HRV interpolation (usually 4 Hz)
% FS_bhi: sampling frequency of the output BHI series (usually < 10 Hz)
% TV: time-varying flag. 0 for punctual estimate, 1 for time-resolved estimate.
% file_output: the function automatically saves a 'BHI' variable into the path specified by 'file_output'

% The variable BHI has the following structure:
% BHI.BtH: brain-to-heart estimates, a matrix of (Nb x 2 x Nt_bhi) dimensions
% BHI.HtB: heart-to-brain estimates, a matrix of (2 x Nb x Nt_bhi) dimensions
% BHI.FS_bhi: sampling frequency of the BHI estimates
% BHI.time: time vector of Nt_bhi elements
% BHI.bands: frequency limits of the EEG bands being analyzed
% BHI.parameters.EEG_PSD_wind: time window through which the EEG-PSD STFT has been calculated
% BHI.parameters.EEG_PSD_step: time step through which the EEG-PSD STFT has been calculated
% BHI.channels: EEGLAB-like structure for EEG channels locations
% BHI.events: EEGLAB-like structure for EEG events
%
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
