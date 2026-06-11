function yNorm = phaseGainNormalize(y, preSym)

preLen = length(preSym);

if length(y) < preLen
    yNorm = y;
    return;
end

yNorm = y;

yPre = yNorm(1:preLen);

%removes residual phase shift
phaseEst = angle(sum(yPre .* conj(preSym)));
yNorm = yNorm * exp(-1j * phaseEst);

yPre = yNorm(1:preLen);

%ajust the power of the output to be compatible with the demodulation
gainEst = sqrt(mean(abs(preSym).^2)) / ...
    (sqrt(mean(abs(yPre).^2)) + 1e-12);
yNorm = yNorm * gainEst;

end