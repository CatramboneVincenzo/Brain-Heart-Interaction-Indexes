function [HeartToBrain, BrainToLF, BrainToHF, time_bhi, window, HeartToBrain_sigma, HeartToBrain_mc]...
    = BHImodel_wT(TFR_EEG, TFR_HRV, FS, RR, win_RR, window, time_tfr, TV)

% This function quantifies directional Brain-Heart Interplay (BHI) 
% through the model proposed by Catrambone et al.(2019) [1].

% INPUT variables:
% TFR_EEG = Time course of EEG power spectral density (PSD). This must be a matrix (Dimension: Channels X time)
%           samples. Each series in each row should be filtered in the desired frequency band of
%           interest (psi)
% TFR_HRV = Time course of HRV PSD (Dimension: 1 X time). This should be filtered in the
% desired frequency band of interest (phi: e.g., LF or HF band)
% FS      = Sampling Frequency of the two TFRs
% RR      = HRV series (expressed in seconds)
% win_RR  = windows length (expressed in seconds) in which the heartbeat generation model (IPFM) is
% reconstructed (default = 15s)
% window  = windows length (in seconds) in which the parameters are calculated (default: window*FS >= 15 )
% TV: time-varying flag. 0 for punctual estimate, 1 for time-resolved estimate.

% OUTPUT variables:
% - HeartToBrain = Functional coupling index (c_rrTOeeg(T)) from 
% HRV Phi-band to EEG Psi-band
% - BrainToHF, BrainToLF  = Functional coupling indices from  
%  EEG Psi-band to  HRV-LF or  HRV-HF bands
% - HeartToBrain_sigma, HeartToBrain_mc = model parameters to be used for fitting evaluation [1]
% 
% This software assumes that input series 
% are all artifact free, e.g., heartbeat dynamics free of algotirhmic and/or physiological artifacts; e.g.
% EEG series free of artifacts from eye blink, movement, etc.
% ---------------------------------------------------------------------------------------------
%  This code implements the theoretical dissertation published in:
%  [1] Catrambone Vincenzo, Alberto Greco, Nicola Vanello, Enzo Pasquale Scilingo,
%  and Gaetano Valenza. "Time-Resolved Directional Brain–Heart Interplay Measurement 
%  Through Synthetic Data Generation Models." 
%  Annals of biomedical engineering 47, no. 6 (2019): 1479-1489.
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

%% checking input variables

if nargin<5
    disp('Three arguments are needed at least: the PSD time course of EEG and HRV signal; and the HRV time series');
    return
elseif nargin >= 5 && nargin < 9
    switch nargin
        case 5
            win_RR = 15; window = ceil(15/FS);
        case 6
            window = ceil(15/FS);
        case 7
            TV = 1;
    end 
elseif nargin > 8
    disp('Too many input arguments');
    return
end

[r,c] = size(TFR_HRV);
if (c==1)&&(r~=1)
    TFR_HRV = TFR_HRV';
elseif ((r~=1)&&(c~=1))||(length(size(TFR_HRV))>2)
    error('Time course of HRV PSD must be a row vector: 1 X time points')
end

if length(size(TFR_EEG))>2
    error('Time course of EEG PSD must be a 2D matrix: Channels X time points')
end
[Nch,Nt] = size(TFR_EEG);
if (Nt==1)&&(Nch~=1)
    TFR_EEG = TFR_EEG';
    [Nch,Nt] = size(TFR_EEG);
end

if (c~=Nt)
    error('The two PSDs, i.e. of EEG and HRV, must be homologously sampled and related to the same time vector, so they must have the same length.');
end

if (log10(abs(median(RR)))<-1)||(log10(abs(median(RR)))>0.7)
    error('HRV signal measure unit must be in seconds!')
end

wind = window*FS;
if TV
if wind < 15
    window = ceil(15/FS);
    wind = window*FS;
    disp(['The time window used for BHI estimation has been modified to ' num2str(window)...
        'secs, as minimum window allowing robust results with the chosen sampling rate']);
end
end

if (numel(TV)>1)||(TV>1)||(TV<0)
    error('TV input must be 0 (no time-varying, only punctual estimate), or 1 (time-resolved estimates required)')
end

%% RR model parameter estimation
omega_lf = 2*pi*0.1;                % LF mean frequency
omega_hf = 2*pi*0.25;               % HF mean Frequency
rr_cum = cumsum(RR);

index_old = 1;
index_new = find(rr_cum > 1 + win_RR,1)-1;
CS = zeros(fix(rr_cum(end)-win_RR),1);
CP = CS;

for i = 1:fix(rr_cum(end)-win_RR)
    HR = 1/mean(RR(index_old+1:index_new));                                % Time-varying Heart Rate 
    gamma = sin(omega_hf/(2*HR))-sin(omega_lf/(2*HR));                     % gamma parameter of the IFPM model
    MM = [sin(omega_hf/(2*HR))*omega_lf*HR/(sin(omega_lf/(2*HR))*4)   -sqrt(2)*omega_lf*HR/(8*sin(omega_lf/(2*HR)));
        -sin(omega_lf/(2*HR))*omega_hf*HR/(sin(omega_hf/(2*HR))*4)    sqrt(2)*omega_hf*HR/(8*sin(omega_hf/(2*HR)))];
    L = max(RR(index_old:index_new))-min(RR(index_old:index_new));         % estimation of the Poincarè plot indices
    W = sqrt(2)*max(abs(RR(index_old+1:index_new)-RR(index_old:index_new-1)));
    CC = 1/gamma*MM*[L; W];
    CS(i) = CC(1);   CP(i) = CC(2);
    index_old = find(rr_cum > i,1);
    index_new = find(rr_cum > i + win_RR,1)-1;
end

% normalization of the parameters for computational reasons and interpolation
% CSr = CS'/std(CS);  CPr = CP'/std(CP);
time_cs = median(1:win_RR)+(1:length(CS));
CSr = interp1(time_cs, CS, time_tfr, 'cubic');
CPr = interp1(time_cs, CP, time_tfr, 'cubic');

TFR_EEG = (TFR_EEG-min(TFR_EEG,[],2))./(max(TFR_EEG,[],2) - min(TFR_EEG,[],2));
TFR_HRV = (TFR_HRV-min(TFR_HRV,[],2))./(max(TFR_HRV,[],2) - min(TFR_HRV,[],2));
CSr = (CSr-min(CSr))./(max(CSr) - min(CSr));
CPr = (CPr-min(CPr))./(max(CPr) - min(CPr));

%% model running for each EEG channel

if Nch > 1
    parfor ch = 1:Nch
        [HeartToBrain(ch,:), BrainToLF(ch,:), BrainToHF(ch,:),HeartToBrain_sigma(ch,:),HeartToBrain_mc(ch,:),time_bhi(ch)] = ...
            BHI_InsideModel(TFR_EEG(ch,:), TFR_HRV, CPr, CSr, wind, time_tfr, TV);
    end
    time_bhi = time_bhi(1);
else
    [HeartToBrain, BrainToLF, BrainToHF, HeartToBrain_sigma,HeartToBrain_mc, time_bhi] = ...
        BHI_InsideModel(TFR_EEG, TFR_HRV, CPr, CSr, wind, time_tfr, TV);
end

end
