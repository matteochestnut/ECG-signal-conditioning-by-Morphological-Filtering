function [cleanECG, detectedBaseline] = WF(ECG)

% [1] ECG signal conditioning by morphological Filtering, Y. Sun
%   K. L. Chan, S. M. Krishnan, 2002
% [2] Processing and Analysis of ECG Signal Using Nonorthogonal Wavelet
% Transform, WW. Dai, Z. Yang, S. L. Lim, O. Mikhailova, J. Chee

% In the paper from which I needed to implement the algorithm [1], a comparison
% between the MMF filter and a wavelet filter is proposed with repsect to
% the baseline removal stage. The wavelet filter is implemented following
% the algorithm described in [2].

    % [2] propose a 12 level wavelet transform, while testing the
    % performance of the filter, I found that using 20 level gives more
    % precision in finding the frequency band where the baseline drift is
    % dominant
    levels = 20;
    % The paper use a Nonorthogonal wavelet, unfotunately Matlab allows
    % only built-in orthogonal and biorthogonal wavelets, so I used the
    % Daubechies 3 wavelet
    [C, L] = wavedec(ECG, levels, 'db3');
    for i = 7 : 1 : levels
        % the seventh layer is the one from which frequency band the baseline
        % drift usually occurs, the for cicle finds the local maxima in
        % each level coefficients (detail coefficients)
        cd = abs(detcoef(C,L,i));
        localPeaks(i) = length(findpeaks(cd));
    end
    
    % a ratio between adjacents number of local maxima per level is
    % computed
    A = zeros([1,levels-1]);
    for i = 1 : 1 : length(A)
        A(i) = localPeaks(i) / localPeaks(i+1); 
    end
    [~, index] = max(A);
    % from level index + 1 is where the baseline drift is dominant so the
    % coefficients of the levels from index + 1 to 20 are put to 0, then
    % the reconstruction is performed
    Cnew = C;
    tmp = cumsum(L);
    Cnew( tmp(end - (levels + 1)) + 1 : tmp(end - (index + 1)) ) = 0;   
    cleanECG=waverec(Cnew,L,'db3');
    detectedBaseline = ECG - cleanECG;
end