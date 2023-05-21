%% filterECG.m is the main script of the project. It lets choose the dataset
% and one ECG from the dataset to be conditioned. Statistics are displaied
% on the console, while original and conditioned ECGs are plotted. Details
% are given through the code. Generally only the variable "dataset" and "i"
% are the ones which need to be modified in order to select a different
% subject and ECG.
clear;

addpath 'Datasets';
addpath 'functions';

% Choose the dataset to use with the variable dataset:
% dataset = 1, dataset_1 with clean ECG;
% dataset = 2, dataset_2 with clean ECG with arrhythmia;
% dataset = 3, mit_bih_arrhytmia_database;
dataset = 1;

if dataset == 1
    load 'dataset_1.mat';
elseif dataset == 2
    load 'dataset_2.mat';
elseif dataset == 3
    load 'mit_bih_arrhythmia_database.mat';
    load 'mitbihQRSnumber.mat';
end

% Select the ECG to filter with the variable i;
% dataset_1, 24 signals;
% dataset_2, 24 signals;
% mit_bih_arrhytmia_database, 48 signals;
i = 1;

ECG = signals(i,:);
t_axis=(0:length(ECG)-1)/Fs;
disp("---------- Dataset " + dataset + ", Observation " + i + " ----------")

% Adding noise to signal from dataset_1 or dataset_2:
if dataset == 1
    % Apart from ECG, Fs and dataset, the arguments passed to the addNoise
    % function are used to initialize the Gaussians, the slanted line and
    % the cosine to produce the noise and the baseline drift. The values
    % used are the one provided by the paper for dataset 1 and 2 respectively
    % in paragraphs 4.1.1 and 4.1.2.
    % A further explanation of the parameters is provided in addNoise.m
    [noiseECG, baseline, noise] = addNoise(ECG, Fs, dataset, 0.01, 0.2, 10, 0.2);
elseif dataset == 2
    [noiseECG, baseline, noise] = addNoise(ECG, Fs, dataset, 0.02, 0.8, 18, 0.1);
elseif dataset == 3
    noiseECG = ECG;
end

% filtering the signal with MMF (Modified Morphological Filter, Paper filter)
[mmfClean, mmfBaseline, mmfNoise] = MMF(noiseECG, Fs);

%filtering the signal with standard MF (Morphological Filter)
[mfClean, mfBaseline, mfNoise] = MF(noiseECG, Fs);

% As dataset 3 can be evaluated only on the based of the Correct Detection
% Rate of the QRS peaks and the WF performs only baseline correction, ECGs
% from dataset 3 are not conditioned with WF.
if dataset == 1 || dataset == 2

    %filtering the signal with WF (Wavelet Filter)
    [wfClean, wfBaseline] = WF(noiseECG);
    
    % Baseline Correction Ratio
    % As baseline correction is performed the same way for both MMF and MF,
    % here the Baseline Correction Ratio is computed only for the first one.
    BCRwf = norm(wfBaseline) / norm(baseline);
    BCRmmf = norm(mmfBaseline) / norm(baseline);
    disp("Baseline Correction Ratio for MMF conditioning: " + BCRmmf);
    disp("Baseline Correction Ratio for WF conditioning: " + BCRwf);
    
    % Noise Suppression Ratio
    NSRmmf = norm(mmfNoise) / norm(noise);
    NSRmf = norm(mfNoise) / norm(noise);
    disp("Noise Suppression Ratio for MMF conditioning: " + NSRmmf);
    disp("Noise Suppression Ratio for MF conditioning: " + NSRmf);
    
    % Signal to Distortion Ratio
    SDRmmf = norm(ECG - mmfClean) / norm(mmfClean);
    SDRmf = norm(ECG - mfClean) / norm(mfClean);
    disp("Signal-to-Noise Ratio for MMF conditioning: " + SDRmmf);
    disp("Signal-to-Noise for MF conditioning: " + SDRmf);
    
    % Plotting results

    figure(1)
    % Here is plotte the original signal corrupted by baseline drift and
    % noise and the same filtered with MMF, MF and WF
    subplot(3, 1, 1)
    plot(t_axis, noiseECG, t_axis, mmfClean)
    title('MMF Conditioning')
    xlabel('s')
    ylabel('mV')
    legend('non conditioned', 'conditioned')
    subplot(3, 1, 2)
    plot(t_axis, noiseECG, t_axis, mfClean)
    title('MF Conditioning')
    xlabel('s')
    ylabel('mV')
    subplot(3, 1, 3)
    plot(t_axis, noiseECG, t_axis, wfClean)
    title('WF Conditioning')
    xlabel('s')
    ylabel('mV')

    figure(2)
    % Here the original baseline drift and noise added to the clean signal
    % are plotted together with the ones detected by the filters
    subplot(3, 2, 1)
    plot(t_axis, baseline, t_axis, mmfBaseline)
    title('Baseline detected by MMF / MF')
    legend('added', 'detected')
    subplot(3, 2, 2)
    plot(t_axis, baseline, t_axis, wfBaseline)
    title('Baseline detected by WF')
    subplot(3, 2, 3)
    plot(t_axis, noise, t_axis, mmfNoise)
    title('Noise detected by MMF')
    subplot(3, 2, 4)
    plot(t_axis, noise, t_axis, mfNoise)
    title('Noise detected by MF')
    subplot(3, 2, 5)
    % Here is plotted the original clean signal with the ones cleaned by
    % MMf and MF
    plot(t_axis, ECG, t_axis, mmfClean)
    legend('original non corrupted', 'conditioned')
    title('MMF conditioning')
    subplot(3, 2, 6)
    plot(t_axis, ECG, t_axis, mfClean)
    title('MF conditioning')

elseif dataset == 3

    % performing QRS detection to evaluate ECGs from dataset 3
    [numQRS, QRS] = QRSdetection(noiseECG, Fs);
    [numQRSmmf, QRSmmf] = QRSdetection(mmfClean, Fs);
    [numQRSmf, QRSmf] = QRSdetection(mfClean, Fs);
    % computing CDR
    CDR = 100 * abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRS)) / mitbihQRSnumber(i);
    CDRmmf = 100 * abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRSmmf)) / mitbihQRSnumber(i);
    CDRmf = 100 * abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRSmf)) / mitbihQRSnumber(i);
    
    disp("Number of QRS complexes in the original signal: " + mitbihQRSnumber(i))
    disp("Correct Detection Rate for QRS peaks in the original signal: " + CDR + "% (" + numQRS + ")");
    disp("Correct Detection Rate for QRS peaks after MMF conditioning: " + CDRmmf + "% (" + numQRSmmf + ")");
    disp("Correct Detection Rate for QRS peaks after MF conditioning: " + CDRmf + "% (" + numQRSmf + ")");

    % Plotting results

    figure(1)
    % Comparison between orginal signal and filtered ones
    subplot(2, 1, 1)
    plot(t_axis, ECG, t_axis, mmfClean)
    legend('non conditioned', 'conditioned')
    title('MMF Conditioning')
    xlabel('s')
    ylabel('mV')
    subplot(2, 1, 2)
    plot(t_axis, ECG, t_axis, mfClean)
    title('MF Conditioning')
    xlabel('s')
    ylabel('mV')

    figure(2)
    subplot(3, 1, 1)
    plot(t_axis, ECG, t_axis, QRS, 'r^','markerfacecolor',[1 0 0])
    title('QRS peaks detected in the original signal')
    subplot(3, 1, 2)
    plot(t_axis, mmfClean, t_axis, QRSmmf, 'r^','markerfacecolor',[1 0 0])
    title('QRS peaks detected after MMF conditioning')
    subplot(3, 1, 3)
    plot(t_axis, mfClean, t_axis, QRSmf, 'r^','markerfacecolor',[1 0 0])
    title('QRS peaks detected after MF conditioning')

end