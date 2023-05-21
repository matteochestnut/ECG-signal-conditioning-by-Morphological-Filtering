function [cleanECG, detectedBaseline, detectedNoise] = MMF(ECG, Fs)

% [1] ECG signal conditioning by morphological Filtering, Y. Sun
%   K. L. Chan, S. M. Krishnan, 2002
% [2] Impulsive Noise Suppression and Background Normalization of
%   Electrocardiogram Signals Using Morphological Operators, CH. H. Chu,
%   E. J. Delp, 1989

% The modified morphological filter proposed in [1] is implemented using
% the functions available in the Image Processing Toolbox. The structuring
% elements have been implemented with the offsetstrel function, as the
% the function strel with which struturing elements are normaly made,
% creates vectors of only 0s and 1s, while the procedure use also the
% vector [0, 1, 5, 1, 0] as structuring element. This structuring element
% has been modified to [0, 1, 2, 1, 0], according to having ECGs which
% amplitude does not reach 5 millivolts as the ones used in [2], where this
% structuring element was proposed. The same thing has been done for the
% standard MF filter.

    %baseline removal
    Bo = offsetstrel(zeros([1,0.2*Fs]));
    Bc = offsetstrel(zeros([1,0.2*Fs*1.5]));
    open = imopen(ECG, Bo);
    detectedBaseline = imclose(open, Bc);
    ECGb = ECG - detectedBaseline;
    
    % noise reduction
    B1 = offsetstrel([0, 1, 2, 1, 0]);
    B2 = offsetstrel(zeros([1,5]));
    f1 = imdilate(ECGb, B1);
    f1 = imerode(f1, B2);
    f2 = imerode(ECGb,B1);
    f2 = imdilate(f2,B2);
    cleanECG = 0.5 * (f1 + f2);
    detectedNoise = ECGb - cleanECG;

end