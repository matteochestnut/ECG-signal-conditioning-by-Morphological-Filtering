function [numQRS, QRS] = QRSdetection(ECG, Fs)

% [1] ECG signal conditioning by morphological Filtering, Y. Sun
%   K. L. Chan, S. M. Krishnan, 2002
% [2] An Approach to QRS Complex Detection Using Mathematical Morphology,
%   P. E. Trahanias, 1993
% [3] QRS Detection Based on Improved Adaptive Threshold, Y. Yu, X. Lu,
%   M. Pan, 2018
% [4] Impulsive Noise Suppression and Background Normalization of
%   Electrocardiogram Signals Using Morphological Operators, CH. H. Chu,
%   E. J. Delp, 1989

% In the paper from which I needed to implement the algorithm, [1], for the
% evaluation of the morphological filtering for the mit bih arrhythmia
% dataset, the QRS complexes of the original ECG and the filterd one are
% detected and then the Correct Detection Rate is computed. The QRS
% detection algorithm used in that paper is the one described by [2] which
% use a morphological filter itself to preprocess the signal. After the
% preprocessing an adaptive threshold technique is used to identify the QRS
% complexes. the technique proposed by [2] to detect QRS complexes is
% described in the paper "Bottom-up approach to the ECG pattern-recognition 
% problem", by the same P. Trahanias et al. Instead I decided to implement
% a much more modern technique described in [3].

    % Preprocessing of the ECG with standard morphological filtering is
    % carried out as in [2], which perform noise reduction and baseline
    % drift removal similarly to [1]. Both the papers implemens a modified
    % version of the standard morphological filter described in [4].
    % Structuring element for noise reduction, the same used in [4] and [1]
    B1 = offsetstrel([0, 1, 1, 1, 0]);
    f1 = imopen(ECG, B1);
    f1 = imclose(f1, B1);
    f2 = imclose(ECG,B1);
    f2 = imopen(f2,B1);
    f = 0.5 * (f1 + f2); % averaging to remove noise
    % Structuring element for baseline correction, in [2] the structuring
    % element used is long accordingly with a temporal length of 22 ms,
    % here it has been built the same way

    offset = zeros([1,int64(0.4*Fs)]);
    B2 = offsetstrel(offset);
    f3 = imopen(f, B2);
    f4 = imclose(f3, B2);
    % Baseline drift removal
    signal = f - f4;
    
    % Now the rest of the preprocessing and the QRS detection stages are
    % implemented accordingly to [3]
    % Differentiate and squaring the signal to amplify and make steeper the
    % QRS complexes
    signalDifferentiated = gradient(signal);
    signalSquared = signalDifferentiated.^2;
    % The integral is perform with a moving windo summation,
    % the window is 66.7 ms long, so 66.7/100 * Fs samples
    signalIntegrated = movsum(signalSquared, (33/1000)*Fs);
    signalIntegrated(signalIntegrated<0.01) = 0;

    % Peak detection phase
    % The QRS peaks are traversed, if a sample has a major amplitude than
    % the previous one then that amplitude value is memorized, when the
    % amplitude drops by half of the previous maximum recorded a QRS peak
    % is recorded in the position of the latter maximum
    peaks = [];
    maximum = signalIntegrated(1);
    maximum_index = 1;
    previous_maximum = 1;
    for i = 2 : 1 : length(signalIntegrated)
        if signalIntegrated(i) > (signalIntegrated(i-1))
            maximum = signalIntegrated(i);
            maximum_index = i;
        end
        if (signalIntegrated(i) < maximum*0.75) && (maximum_index - previous_maximum > 222*Fs/1000)
            peaks(maximum_index) = maximum;
            previous_maximum = maximum_index;
        end
    end
    
    % The peaks registered are traversed with an adaptive threshold
    maximum_val = max(peaks(1:300));
    signal_peak = 0.13 * maximum_val;
    noise_peak = 0.1 * signal_peak;
    threshold = 0.25 * signal_peak + 0.75 * noise_peak;

    for j = 1 : 1 : length(peaks)
        if peaks(j) > threshold
            signal_peak = 0.13 * peaks(j);
            QRS(j) = 1;
        else
            noise_peak = 0.1 * signal_peak;
        end
        % Threshold update
        threshold = 0.25 * signal_peak + 0.75 * noise_peak;
    end

    % QRS complexes detected are counted
    numQRS = sum(QRS);

    % filling QRS vector
    QRS = [QRS zeros([1, length(ECG)-length(QRS)])];
    % substituting the peaks with their actual value in the ECG for
    % plotting
    for k = 1 : 1 : length(QRS)
        if QRS(k) == 1
            QRS(k) = abs(ECG(k));
        end
        if QRS(k) == 0
            QRS(k) = nan;
        end
    end

end