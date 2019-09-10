function [HeartToBrain, BrainToLF, BrainToHF, HeartToBrain_sigma, HeartToBrain_mc]...
    = BHImodel(TFR_EEG, TFR_HRV, FS, RR, win_RR, window)
% BHImodel2() is a function that quantify directional Brain-Heart Interplay (BHI) 
% adopting the BHI model proposed by Catrambone et al.(2019) [1].
% 
% INPUT variables:
% TFR_EEG = Time course of EEG PSD, it must be a matrix: Channels X time
%           points. It should be already filtered in the desired band of
%           interest (psi)
% TFR_HRV = Time course of HRV PSD. It should be already filtered in the
% desired band of interest. (phi)
% FS      = Sampling Frequency of the two TFRs
% RR      = HRV time series (in seconds)
% win_RR  = windows length (in seconds) in which the RR model is
% reconstructed (default = 15 secs)
% window  = windows length (in seconds) in which the parameters are
% calculated (by default, it is assured that: window*FS >= 15 )
% 
% OUTPUT variables:
% - HeartToBrain = Functional coupling index (c_rrTOeeg(T)) from the given
% Phi-band of HRV to the given Psi-band of the EEG
% - BrainToHF, BrainToLF  = Functional coupling indexes from the given 
% Psi-band of EEG to the LF or HF HRV bands
% - HeartToBrain_sigma, HeartToBrain_mc = other parameters that can be
%   extracted from the model in [1]
% 
% It is assumed that the signals given as input the BHImodel2() function
% are artifact free, from algotirhmic and physiological point of view (e.g.
% ectopic beats, eye blink, movement, ecc.)
% ---------------------------------------------------------------------------------------------
%  This code is an implementation of the following scientific research,
%  where an interested reader can find all the mathematical and teorethical details:
%  [1] Catrambone Vincenzo, Alberto Greco, Nicola Vanello, Enzo Pasquale Scilingo,
%  and Gaetano Valenza. "Time-Resolved Directional Brain–Heart Interplay Measurement 
%  Through Synthetic Data Generation Models." 
%  Annals of biomedical engineering 47, no. 6 (2019): 1479-1489.
% ---------------------------------------------------------------------------------------------
% Copyright (C) 2019-2020 Vincenzo Catrambone, Gaetano Valenza
% 
% This program is free software; you can redistribute it and/or modify it under
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
% The software does not come with a GUI. Assuming 'x' and 'y' being the
% time course, sampled at 1 Hz, of Power spectral density of an EEG channel 
% in the theta band (4-8 Hz) and HRV series respetively (integrated in LF+HF bands),
% and 'z' being the original HRV series, the following example performs the BHImodel 
% analysis (with default parameters) and plots the results:
% 
% Fs = 1;
% [HeartToTheta, ThetaToHF, ThetaToLF] = BHImodel(x,y,Fs,z);
% figure, hold all
% plot(1:length(HeartToTheta))'/Fs, HeartToTheta);
% plot(1:length(ThetaToHF))'/Fs, ThetaToHF);
% plot(1:length(ThetaToHF))'/Fs, ThetaToLF);
% ---------------------------------------------------------------------------------------------

%% checking input variables

if nargin<5
    disp('Three arguments are needed at least: the PSD time course of EEG and HRV signal; and the HRV time series');
    return
elseif nargin >= 5 && nargin < 6
    switch nargin
        case 5
            win_RR = 15; window = ceil(15/FS);
        case 6
            window = ceil(15/FS);
    end 
elseif nargin > 6
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
if wind < 15
    window = ceil(15/FS);
    wind = window*FS;
    disp(['The time window used for BHI estimation has been modified to ' num2str(window)...
        'secs, as minimum window allowing robust results with the chosen sampling rate']);
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
CSr = CS'/std(CS);  CPr = CP'/std(CP);
CSr = interp1(1:length(CSr), CSr, 1/FS:1/FS:Nt/FS, 'spline');
CPr = interp1(1:length(CPr), CPr, 1/FS:1/FS:Nt/FS, 'spline');
TFR_EEG = (sqrt(TFR_EEG));

%% model running for each EEG channel

if Nch > 1
    parfor ch = 1:Nch
        [HeartToBrain(ch,:), BrainToLF(ch,:), BrainToHF(ch,:),HeartToBrain_sigma(ch,:),HeartToBrain_mc(ch,:)] = ...
            BHI_InsideModel2(TFR_EEG(ch,:), TFR_HRV, CPr, CSr, wind);
    end
else
    [HeartToBrain, BrainToLF, BrainToHF, HeartToBrain_sigma,HeartToBrain_mc] = ...
        BHI_InsideModel2(TFR_EEG, TFR_HRV, CPr, CSr, wind);
end

end

function [HToB, BToLF, BToHF, HToB_sigma, HToB_mc] = BHI_InsideModel2(TFR_ch, TFR_rr, CPr, CSr, window)

Nt = length(TFR_ch);
Cs1 = 0.25; Cp1 = 0.24;

for i = 1:window
    
%% EEG parameter estimation, during the first period in which the RR model connot be calculated

arx_data = iddata(TFR_ch(i:i+window)', TFR_rr(i:i+window)',1);                      % iddata-format is necessary for the arx function
model_eegP = arx(arx_data,[1 1 1]);                                                 % here the model is estimated
HToB_sigma(i) = sqrt(model_eegP.NoiseVariance); HToB_mc(i) = -model_eegP.A(2);      % the parameters are extracted
HToB(i) = model_eegP.B(2);
medianTime_P_eeg(1,i) = median(TFR_ch(i:i+window));                                 % this operation is needed for the following RR-model
end
%% after a first period (equal to the window length), both the inverse models are running, so the parameters are parallelly calculated 
for i = window+1:min([length(CPr),Nt-window, length(TFR_rr)-window])
    
    %% ARX EEG parameter estimation
    arx_data = iddata(TFR_ch(i:i+window)', TFR_rr(i:i+window)',1); 
    model_eegP = arx(arx_data,[1 1 1]);
    HToB_sigma(i) = sqrt(model_eegP.NoiseVariance); HToB_mc(i) = -model_eegP.A(2);
    HToB(i) = model_eegP.B(2);
    medianTime_P_eeg(1,i) = median(TFR_ch(i:i+window));

    %% IPFM RR parameter estimation, the interaction are modelled separately with all the EEG bands
    if i-window <= length(CPr)-window-1
        % Brain to HF
        BToHF(i-window) = median((CPr(i-window:i)-Cp1)./medianTime_P_eeg(i-window:i));
        % Brain to LF
        BToLF(i-window) = median((CSr(i-window:i)-Cs1)./medianTime_P_eeg(i-window:i));
    else
        % Brain to HF
        BToHF(i-window) = BToHF(i-window-1);
        % Brain to LF
        BToLF(i-window) = BToLF(i-window-1);
    end
end

end