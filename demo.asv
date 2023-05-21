clear; clc;

addpath 'Datasets';
addpath 'functions';

dataset = 1;

if dataset == 1
    load 'dataset_1.mat';
elseif dataset == 2
    load 'dataset_2.mat';
elseif dataset == 3
    load 'mit_bih_arrhythmia_database.mat';
    load 'mitbihQRSnumber.mat';
end

i = 1;

ECG = signals(i,:);
t_axis=(0:length(ECG)-1)/Fs;

if dataset == 1
    [noiseECG, baseline, noise] = addNoise(ECG, Fs, dataset, 0.01, 0.2, 10, 0.2);
elseif dataset == 2
    [noiseECG, baseline, noise] = addNoise(ECG, Fs, dataset, 0.02, 0.8, 18, 0.1);
elseif dataset == 3
    noiseECG = ECG;
end

[mmfClean, mmfBaseline, mmfNoise] = MMF(noiseECG, Fs);

[mfClean, mfBaseline, mfNoise] = MF(noiseECG, Fs);

if dataset == 1 || dataset == 2

    %filtering the signal with WF (Wavelet Filter)
    [wfClean, wfBaseline] = WF(noiseECG);
    
    BCRwf = norm(wfBaseline) / norm(baseline);
    BCRmmf = norm(mmfBaseline) / norm(baseline);
    
    % Noise Suppression Ratio
    NSRmmf = norm(mmfNoise) / norm(noise);
    NSRmf = norm(mfNoise) / norm(noise);
    
    % Signal to Distortion Ratio
    SDRmmf = norm(ECG - mmfClean) / norm(mmfClean);
    SDRmf = norm(ECG - mfClean) / norm(mfClean);
    
    figure(1)
    subplot(2,1,1)
    plot(t_axis, baseline, t_axis, wfBaseline)
    title('WF detected baseline')
    subplot(2,1,2)
    plot(t_axis, baseline, t_axis, mmfBaseline)
    title('MMF/MF detected baseline')



elseif dataset == 3

    % performing QRS detection to evaluate ECGs from dataset 3
    [numQRS, QRS] = QRSdetection(noiseECG, Fs);
    [numQRSmmf, QRSmmf] = QRSdetection(mmfClean, Fs);
    [numQRSmf, QRSmf] = QRSdetection(mfClean, Fs);
    % computing CDR
    CDR = 100 * abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRS)) / mitbihQRSnumber(i);
    CDRmmf = 100 * abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRSmmf)) / mitbihQRSnumber(i);
    CDRmf = 100 * abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRSmf)) / mitbihQRSnumber(i);
    

end


%%
 figure(2)
 subplot(1,2,1)
 plot(t_axis, ECG, t_axis, mmfClean)
 title('MMF')
 subplot(1,2,2)
 plot(t_axis, ECG, t_axis, mfClean)
 title('MF')

%%
load 'mit_bih_arrhythmia_database.mat';
ECG2 = signals(2,:);
[mmfClean2, mmfBaseline, mmfNoise] = MMF(ECG2, Fs);

[mfClean2, mfBaseline, mfNoise] = MF(ECG2, Fs);


figure(2)

subplot(2,2,1)
plot(t_axis, noiseECG, t_axis, mmfClean)
title('MMF (dataset 1, ECG 1)')

subplot(2,2,2)
plot(t_axis, noiseECG, t_axis, mfClean)
title('MF (dataset 1, ECG 1)')

t_axis=(0:length(ECG2)-1)/Fs;

subplot(2,2,3)
plot(t_axis, ECG2, t_axis, mmfClean2)
title('MMF (dataset 3, ECG 1)')

subplot(2,2,4)
plot(t_axis, ECG2, t_axis, mfClean2)
title('MF (dataset 3, ECG 1)')
