function [w, eTrain] = spatialLMSRefineSteeringSamples(Xs, dSamples, p, wInit)

% Spatial LMS applied at sample rate.
%
% Xs: Nsamples x M matrix
% dSamples: sample-rate reference signal
% wInit: M x 1 initial steering weights

[N, M] = size(Xs);

trainLen = min(length(dSamples), N);

w = wInit(:);

if length(w) ~= M
    error('wInit must have one coefficient per hydrophone.');
end

eTrain = zeros(trainLen, 1);

for n = 1:trainLen

    x = Xs(n,:).';

    y_n = w' * x;

    eTrain(n) = dSamples(n) - y_n;

    if isfield(p, 'normalizeLMS') && p.normalizeLMS
        px = sum(abs(x).^2) + 1e-8;
        step = p.mu / px;
    else
        step = p.mu;
    end

    w = w + step * x * conj(eTrain(n));
end

end