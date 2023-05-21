function [signal_new, baseline, noise] = addNoise(signal, Fs, dataset, theta, A, K, e)

% [1] ECG signal conditioning by morphological Filtering, Y. Sun
%   K. L. Chan, S. M. Krishnan, 2002

% theta, slope of baseline drift
% A amplitude of the cosine to be added to the slanted line
% e, scale factor for the Gaussians
% K, ratio between the standard deviations of the Gaussians, v1 and v2,
%   v1 is always fixed to 0.1

    N = length(signal); %ECG length
    t_axis=(0:N-1)/Fs;
    
    v1 = 0.1;
    v2 = v1 * K;
    G1 = normrnd(0, v1, [1,N]); % simulate background noise
    G2 = normrnd(0, v2, [1,N]); % simulate impulsive noises
    % In [1] the scaling factor 0.3 for the noise is not present, I decide to adopt
    % it in order to produce noise which amplitude was as similar as possible
    % to the one plotted in the paper's figures which is much smaller from
    % the one actually produced by the equation provided in [1]
    noise = 0.3 * ((1 - e) * G1 + e * G2); % noise
    
    % coefficient of the slanted line which simulate the baseline drift
    m = tan(theta);
    % B, bias of the baseline drift
    if dataset == 1
        B = -0.6;
    else
        B = 0;
    end
    slanted_line = t_axis .* m;
    % The cosine frequency has been chose arbitrarily, the paper does not
    % provide any detail on the frequency used, so a popular frequency for
    % an actual baseline drift has been used.
    sinusoidal_sig = cos(2*pi*t_axis*0.5 + theta);
    baseline = B + slanted_line + A .* sinusoidal_sig;

    signal_new = signal + noise + baseline;

end