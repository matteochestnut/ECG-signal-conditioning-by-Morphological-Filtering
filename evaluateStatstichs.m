%% Evalueting statistics for datasets 1 and 2:
% 1) fixing theta = zero,  A varies between 0.2 and 1, average BCR across
%   all the observations is evaluated for WF and MMF;
% 2) fixing A = 0, theta varies between 15 and 75 degrees, average BCR
%   across all the observations is evaluated for WF and MMF;
% 3) fixing K = 10 (the ratio between v1 and v2, the std. dev, of the
%   two gaussians making the noise, with v1 = 0.1), e varies between 0.1
%   and 0.5, average NSR and SDR across all the observations is evaluated
%   for MF and MMF;
% 4) fixing e = 0.2, K varies between 5 an 25 (v1 = 0.1), average NSR and
%   SDR across all the observations is evaluated for MF and MMF.

% Evaluating statistics for mit bih arrhythmia dataset:
% Computing average correct detection rate (CDR) for QRS peaks of
% the original signals, signals filtered with WF and signals filtered
% with MMF. Given a signal from the dataset, its number of QRS peaks noted
% in its metadata; the number of QRS peaks of all the signals is save in a
% vector in mitbohQRSnumber.mat. CDR evaluation takes into account only the
% difference between the number of QRS peaks detected and the actual number
% of QRS peaks of the signal, true positive and false positive detections
% are not defined.

clear;

addpath 'Datasets';
addpath 'functions';

% choose the dataset for which evaluate the statistics:
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

%% Choose statistics to evaluate, if datset 1 or 2 is chose then stat goes
% from 1 to 4 accordgly with what written on top, if mit bih dataset is
% chose, stat is not used.
stat = 1;

if dataset == 1 || dataset == 2
    if stat == 1
        for i = 1 : 1 : size(signals, 1)
            count = 1;
            ECG = signals(i,:);
            for A = 0.2 : 0.1 : 1
                [noiseECG, baseline, ~] = addNoise(ECG, Fs, dataset, 0, A, 10, 0.2);
                [mmfClean, mmfBaseline, ~] = MMF(noiseECG, Fs);
                [wfClean, wfBaseline] = WF(noiseECG);
                BCRmmf(i,count) = norm(mmfBaseline) / norm(baseline);
                BCRwf(i,count) = norm(wfBaseline) / norm(baseline);
                count = count + 1;
            end
        end
        BCRmmf = mean(BCRmmf, 1);
        BCRwf = mean(BCRwf, 1);
        
        disp("---------- Dataset " + dataset + " ----------")
        disp("For theta = 0 and A varying between 0.2 and 1")
        disp("Average Baseline Correction Rate for MMF conditioning: " + mean(BCRmmf));
        disp("Average Baseline Correction Rate for WF conditioning: " + mean(BCRwf));
        figure(1)
        plot(0.2 : 0.1 : 1, BCRmmf, 0.2 : 0.1 : 1, BCRwf)
        legend('BCR for MMF', 'BCR for WF')
        xlabel('A')
        ylabel('BCR')
    elseif stat == 2
        for i = 1 : 1 : size(signals, 1)
            count = 1;
            ECG = signals(i,:);
            for theta = 15 : 1 : 75
                [noiseECG, baseline, ~] = addNoise(ECG, Fs, dataset, theta, 0, 10, 0.2);
                [mmfClean, mmfBaseline, ~] = MMF(noiseECG, Fs);
                [wfClean, wfBaseline] = WF(noiseECG);
                BCRmmf(i,count) = norm(mmfBaseline) / norm(baseline);
                BCRwf(i,count) = norm(wfBaseline) / norm(baseline);
                count = count + 1;
            end
        end
        BCRmmf = mean(BCRmmf, 1);
        BCRwf = mean(BCRwf, 1);

        disp("---------- Dataset " + dataset + " ----------")
        disp("For A = 0 and theta varying between 15 and 75")
        disp("Average Baseline Correction Rate for MMF conditioning: " + mean(BCRmmf));
        disp("Average Baseline Correction Rate for WF conditioning: " + mean(BCRwf));
        figure(1)
        plot(15 : 1 : 75, BCRmmf, 15 : 1 : 75, BCRwf)
        legend('BCR for MMF', 'BCR for WF')
        xlabel('theta')
        ylabel('BCR')
    elseif stat == 3
        for i = 1 : 1 : size(signals, 1)
            count = 1;
            ECG = signals(i,:);
            for e = 0.1 : 0.1 : 0.5
                [noiseECG, ~, noise] = addNoise(ECG, Fs, dataset, 0.01, 0.2, 10, e);
                [mmfClean, ~, mmfNoise] = MMF(noiseECG, Fs);
                [mfClean, ~, mfNoise] = MF(noiseECG, Fs);
                NSRmmf(i,count) = norm(mmfNoise) / norm(noise);
                NSRmf(i,count) = norm(mfNoise) / norm(noise);
                SDRmmf(i,count) = norm(ECG - mmfClean) / norm(mmfClean);
                SDRmf(i, count) = norm(ECG - mfClean) / norm(mfClean);
                count = count + 1;
            end
        end
        NSRmmf = mean(NSRmmf, 1);
        NSRmf = mean(NSRmf, 1);
        SDRmmf = mean(SDRmmf, 1);
        SDRmf = mean(SDRmf, 1);

        disp("---------- Dataset " + dataset + " ----------")
        disp("For K = 10 and e varying between 0.1 and 0.5")
        disp("Average Noise Suppression Rate for MMF conditioning: " + mean(NSRmmf));
        disp("Average Noise Suppression Rate for MF conditioning: " + mean(NSRmf));
        disp("Average Signal Distortion Rateo for MMF conditioning: " + mean(SDRmmf));
        disp("Average Signal Distortion Rateo for MF conditioning: " + mean(SDRmf));
        figure(1)
        plot(0.1 : 0.1 : 0.5, NSRmmf, 0.1 : 0.1 : 0.5, NSRmf)
        legend('NSR for MMF', 'NSR for MF')
        xlabel('e')
        ylabel('NSR')
        figure(2)
        plot(0.1 : 0.1 : 0.5, SDRmmf, 0.1 : 0.1 : 0.5, SDRmf)
        legend('SDR for MMF', 'SDR for MF')
        xlabel('e')
        ylabel('SDR')
    elseif stat == 4
        for i = 1 : 1 : size(signals, 1)
            count = 1;
            ECG = signals(i,:);
            for K = 5 : 1 : 25
                [noiseECG, ~, noise] = addNoise(ECG, Fs, dataset, 0.01, 0.2, K, 0.2);
                [mmfClean, ~, mmfNoise] = MMF(noiseECG, Fs);
                [mfClean, ~, mfNoise] = MF(noiseECG, Fs);
                NSRmmf(i,count) = norm(mmfNoise) / norm(noise);
                NSRmf(i,count) = norm(mfNoise) / norm(noise);
                SDRmmf(i,count) = norm(ECG - mmfClean) / norm(mmfClean);
                SDRmf(i, count) = norm(ECG - mfClean) / norm(mfClean);
                count = count + 1;
            end
        end
        NSRmmf = mean(NSRmmf, 1);
        NSRmf = mean(NSRmf, 1);
        SDRmmf = mean(SDRmmf, 1);
        SDRmf = mean(SDRmf, 1);

        disp("---------- Dataset " + dataset + " ----------")
        disp("For e = 0.2 and K varying between 5 and 25")
        disp("Average Noise Suppression Rate for MMF conditioning: " + mean(NSRmmf));
        disp("Average Noise Suppression Rate for MF conditioning: " + mean(NSRmf));
        disp("Average Signal Distortion Rateo for MMF conditioning: " + mean(SDRmmf));
        disp("Average Signal Distortion Rateo for MF conditioning: " + mean(SDRmf));
        figure(1)
        plot(5 : 1 : 25, NSRmmf, 5 : 1 : 25, NSRmf)
        legend('NSR for MMF', 'NSR for MF')
        xlabel('K')
        ylabel('NSR')
        figure(2)
        plot(5 : 1 : 25, SDRmmf, 5 : 1 : 25, SDRmf)
        legend('SDR for MMF', 'SDR for MF')
        xlabel('K')
        ylabel('SDR')
    end
else
    for i = 1 : 1 : size(signals, 1)
        noiseECG = signals(i,:);
        [mmfClean(i,:), ~, ~] = MMF(noiseECG, Fs);
        [mfClean(i,:), ~, ~] = MF(noiseECG, Fs);
        [numQRS(i), ~] = QRSdetection(noiseECG, Fs);
        [numQRSmmf(i), ~] = QRSdetection(mmfClean(i,:), Fs);
        [numQRSmf(i), ~] = QRSdetection(mfClean(i,:), Fs);
        CDR(i) = 100 .* abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRS(i))) / mitbihQRSnumber(i);
        CDRmmf(i) = 100 .* abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRSmmf(i))) / mitbihQRSnumber(i);
        CDRmf(i) = 100 .* abs(mitbihQRSnumber(i) - abs(mitbihQRSnumber(i) - numQRSmf(i))) / mitbihQRSnumber(i);
    end
    CDR = mean(CDR);
    CDRmmf = mean(CDRmmf);
    CDRmf = mean(CDRmf);
    disp( "---------- Dataset " + dataset + " ----------" + newline + ...
        "Average CDR for originals signal in mit bih arrhythmia dataset (no conditioning): " + ...
        CDR + newline + "Average CDR for signals filtered with MMF: " + ...
        CDRmmf + newline + "Average CDR for signal filtered with MF: " + CDRmf);
end