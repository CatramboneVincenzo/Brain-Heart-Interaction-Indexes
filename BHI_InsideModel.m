function [HToB, BToLF, BToHF, HToB_sigma, HToB_mc, time_bhi] = BHI_InsideModel(TFR_ch, TFR_rr, CPr, CSr, window, time_tfr, TV)

Nt = length(TFR_ch);
Cs1 = 0; % .25; 
Cp1 = 0; 

if TV

for i = 1:min([length(CPr),Nt-window, length(TFR_rr)-window])
    
%% EEG parameter estimation, during the first period in which the RR model connot be calculated
    arx_data = iddata(TFR_ch(i:i+window)', TFR_rr(i:i+window)',1);                      % iddata-format is necessary for the arx function
    model_eegP = arx(arx_data,[1 1 1]);                                                 % here the model is estimated
    HToB_sigma(i) = sqrt(model_eegP.NoiseVariance); HToB_mc(i) = -model_eegP.A(2);      % the parameters are extracted
    HToB(i) = model_eegP.B(2);
    time_HToB(i) = time_tfr(window/2+i);
    medianTime_P_eeg(1,i) = median(TFR_ch(i:i+window));     % this operation is needed for the following RR-model

    try        % Brain to HF
        BToHF(i) = (CPr(i)-Cp1)./medianTime_P_eeg(i);
        % Brain to LF
        BToLF(i) = (CSr(i)-Cs1)./medianTime_P_eeg(i);
    catch
        % Brain to HF
        BToHF(i) = BToHF(i-1);
        % Brain to LF
        BToLF(i) = BToLF(i-1);
    end
end
    time_BToH = time_HToB(~isnan(BToLF));
    BToLF = BToLF(~isnan(BToLF));
    BToHF = BToHF(~isnan(BToHF));
    
    HToB = conv(HToB,2/window*ones(1,window/2),'same');
    BToHF = conv(BToHF,2/window*ones(1,window/2),'same');
    BToLF = conv(BToLF,2/window*ones(1,window/2),'same');

else

    arx_data = iddata(TFR_ch', TFR_rr',1);                      % iddata-format is necessary for the arx function
    model_eegP = arx(arx_data,[1 1 1]);                                                 % here the model is estimated
    HToB_sigma = sqrt(model_eegP.NoiseVariance); 
    HToB_mc = -model_eegP.A(2);      % the parameters are extracted
    HToB = model_eegP.B(2);
    time_HToB = 1;
    medianTime_P_eeg = nanmedian(TFR_ch);     % this operation is needed for the following RR-model

    % Brain to HF
    BToHF = nanmedian((CPr-Cp1)./TFR_ch);
    % Brain to LF
    BToLF = nanmedian((CSr-Cs1)./TFR_ch);
    time_BToH = 1;
end

    time_bhi.time_HToB = time_HToB;
    time_bhi.time_BToH = time_BToH;

end