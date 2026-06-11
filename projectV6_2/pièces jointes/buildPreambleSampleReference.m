function dSamples = buildPreambleSampleReference(preSym, txFilter, rxFilter, p)

preTx = upfirdn(preSym, txFilter, p.sps, 1);

preMf = upfirdn(preTx, rxFilter, 1, 1);

totalDelay = length(rxFilter) - 1;

dSamples = preMf(totalDelay+1:end);

targetLen = length(preSym) * p.sps;

if length(dSamples) >= targetLen
    dSamples = dSamples(1:targetLen);
end

dSamples = dSamples(:);

end