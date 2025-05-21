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

[R,C] = size(f_lims);
if (R==2)&&(C~=2)
    f_lims = f_lims';
elseif (C~=2)
    error('Variable f_lims must contain on each row lowerbound and upperbound of the frequency bands which are going to be analyzed');
end

Nbands = size(f_lims,1);
window = 2; % seconds
win_RR = 15; % RR samples 
Nch = EEG.nbchan;

w = EEG.srate*window; % window size (2 second) in #samples
step = round(EEG.srate/FS_bhi); % step size in #samples

%% PSD EEG
n = length(EEG.data(1,:)); % length of the signal
l = floor((n-w)/step);  % # number of segments
PSD_EEG = zeros(Nch,l,Nbands);
for bb = 1:Nbands
    for i = 1:l
        PSD_EEG(:,i,bb) = bandpower(EEG.data(:,((i-1)*step+1:step*(i-1)+w))',EEG.srate,f_lims(bb,:))';
    end
end
t_tfr_eeg = (0:l-1)/FS_bhi+1;

%% HRV
TFR_RR = wvd(hilbert(RRi-mean(RRi)),FS_rri,'smoothedPseudo');
f = linspace(0, FS_bhi/2, size(TFR_RR,1));
t_tfr_rr = linspace(t_RRi(1),t_RRi(end),size(TFR_RR,2));
% t_bhi = t_RRi(1):1/FS_bhi:t_RRi(end);

[~,LFidx] = find((0.04<=f)&(f<=0.15));
[~,HFidx] = find((0.15<=f)&(f<=0.4));

LF = trapz(TFR_RR(LFidx,:));
HF = trapz(TFR_RR(HFidx,:));

PSD_HRV(1,:) = interp1(t_tfr_rr,LF,t_tfr_eeg,'cubic');
PSD_HRV(2,:) = interp1(t_tfr_rr,HF,t_tfr_eeg,'cubic');

Ind_toKeep = (~isnan(PSD_HRV(1,:)))&(~isnan(PSD_HRV(2,:)));

%% BHI

time_tfr = t_tfr_eeg(Ind_toKeep);
PSD_EEG = PSD_EEG(:,Ind_toKeep,:);
PSD_HRV = PSD_HRV(:,Ind_toKeep);

%% SDG model 1: LF HF
for eeg_b = 1:Nbands
for rr_b = 1:2
    [HtB(eeg_b,rr_b,:,:), BtH(eeg_b,1,:,:), BtH(eeg_b,2,:,:), time_bhi, window, ~, ~] = ...
    BHImodel_wT(squeeze(PSD_EEG(:,:,eeg_b)), PSD_HRV(rr_b,:), FS_bhi, RR, win_RR, window, time_tfr, TV); 
end
end

BHI.BtH = BtH;
BHI.HtB = HtB;
BHI.FS_bhi = FS_bhi;
BHI.time = time_bhi;
BHI.bands = f_lims;
BHI.parameters.BHI_window = window;
BHI.parameters.EEG_PSD_wind = w;
BHI.parameters.EEG_PSD_step = step;
BHI.channels = EEG.chanlocs;

try
    BHI.events = EEG.event;
catch
    
end
    
%% save files
save(file_output, 'PSD_EEG', 'PSD_HRV', 'BHI');
