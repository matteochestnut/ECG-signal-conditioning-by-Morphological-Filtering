function [cleanECG, detectedBaseline, detectedNoise] = MF(ECG, Fs)

% [1] Impulsive Noise Suppression and Background Normalization of
%   Electrocardiogram Signals Using Morphological Operators, CH. H. Chu,
%   E. J. Delp, 1989
% [2] ECG signal conditioning by morphological Filtering, Y. Sun
%   K. L. Chan, S. M. Krishnan, 2002

% A standard morphological filter to remove baseline drift and perform
% noise reduction is applied as described in [1], some variations have been
% apllied: in [1] the noise reduction stage is performed first, but
% to produce comparable result it has been switched with the base drift
% removal stage as in [2]; also the baseline drift removal in [1] is
% performed by averaging the result of two operation: 1) opening - closing
% on the signal, 2) closing opening on the signal; in the paper itself it
% is noted that only one branch of the procedure can be performed, in [2]
% only the first operation (opening - closing) is performed, so it is in
% this code. The last change is using the vector [0, 1, 1, 1, 0] as
% structuring element for the noise reduction stage instead of [0, 1, 5, 1,
% 0], this change is due to a performance improvement noticed while testing
% the MF filter. In the original paper from which the triangle structuring
% element is propose, [1], ECGs reach an amplitude of 5 millivolts. Despite
% [2] still using [0, 1, 5, 1, 0] as structuring element, I noticed that
% the filter was completely unable to clean the ECG from any noise; from
% this the changing in the structuring element.

    %baseline removal
    Bo = offsetstrel(zeros([1,0.2*Fs]));
    Bc = offsetstrel(zeros([1,0.2*Fs*1.5]));
    open = imopen(ECG, Bo);
    detectedBaseline = imclose(open, Bc);
    ECGb = ECG - detectedBaseline;
    
    % noise reduction
    B = offsetstrel([0, 1, 1, 1, 0]);
    f1 = imopen(ECGb, B);
    f1 = imclose(f1, B);
    f2 = imclose(ECGb,B);
    f2 = imopen(f2,B);
    cleanECG = 0.5 * (f1 + f2);
    detectedNoise = ECGb - cleanECG;

end