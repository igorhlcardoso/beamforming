function rx = rxFrontEndSIMO(rxWave, preambleBits, p, rxFilter)

[~, M] = size(rxWave);

preSym = qpskBitsToSymbols(preambleBits, p);

% For identical TX/RX FIR filters, total group delay is:
% groupDelayTx + groupDelayRx = length(filter) - 1
totalDelay = length(rxFilter) - 1;

rx.frameSymbolsCell = cell(1,M);
rx.corrVecCell = cell(1,M);
rx.frameStartSym = zeros(1,M);
rx.rxMfNoTransientCell = cell(1,M);
rx.rxSymNoCorrCell = cell(1,M);
rx.bestPhase = zeros(1,M);

frameSymCell = cell(1,M);
minFrameLen = inf;

if ~isfield(p, 'useBestSamplingPhase')
    p.useBestSamplingPhase = true;
end

for m = 1:M
%matched filter
    rxMf = upfirdn(rxWave(:,m), rxFilter, 1, 1);%up and downsampling does not happen here, they are 1 and 1

    if length(rxMf) <= totalDelay

        rxMfNoTransient = complex([]);
        rxSymNoCorr = complex([]);
        corrVec = [];
        frameStartSym = 1;
        frameSymbols = complex([]);
        bestPhase = 1;

    else
%remove the filter delay
        rxMfNoTransient = rxMf(totalDelay+1:end);

        if p.useBestSamplingPhase

            bestMetric = -inf;
            bestCorrVec = [];
            bestFrameStartSym = 1;
            bestRxSymNoCorr = complex([]);
            bestPhase = 1;
% sampling phases that are possible
            for ph = 1:p.sps

                rxSymCandidate = rxMfNoTransient(ph:p.sps:end);%downsampling

                if length(rxSymCandidate) < length(preSym)
                    continue;
                end
%calculate the correlation with the preamble
                corrCandidate = abs(conv(rxSymCandidate, ...
                    flipud(conj(preSym)), 'valid'));
%choose the phase with the highest peak
                [metric, idx] = max(corrCandidate);

                if metric > bestMetric
                    bestMetric = metric;
                    bestCorrVec = corrCandidate;
                    bestFrameStartSym = idx;
                    bestRxSymNoCorr = rxSymCandidate;
                    bestPhase = ph;
                end
            end

            rxSymNoCorr = bestRxSymNoCorr;
            corrVec = bestCorrVec;
            frameStartSym = bestFrameStartSym; 

        else

            bestPhase = 1;

            rxSymNoCorr = rxMfNoTransient(bestPhase:p.sps:end);

            if length(rxSymNoCorr) >= length(preSym)
                corrVec = abs(conv(rxSymNoCorr, ...
                    flipud(conj(preSym)), 'valid'));

                [~, frameStartSym] = max(corrVec);
            else
                corrVec = [];
                frameStartSym = 1;
            end
        end

        numFrameSym = length(preSym) + p.payloadLengthBits/p.k;
        lastIdx = frameStartSym + numFrameSym - 1;

        if lastIdx <= length(rxSymNoCorr)
            frameSymbols = rxSymNoCorr(frameStartSym:lastIdx);
        else
            frameSymbols = rxSymNoCorr(frameStartSym:end);
        end
    end

    rx.frameSymbolsCell{m} = frameSymbols; %the beggining of the frame
    rx.corrVecCell{m} = corrVec;
    rx.frameStartSym(m) = frameStartSym;
    rx.rxMfNoTransientCell{m} = rxMfNoTransient;
    rx.rxSymNoCorrCell{m} = rxSymNoCorr;
    rx.bestPhase(m) = bestPhase; %store the best phase

    frameSymCell{m} = frameSymbols;
    minFrameLen = min(minFrameLen, length(frameSymbols));
end

if isinf(minFrameLen) || minFrameLen == 0
    rx.frameSymbolsMatrix = complex([]);
else
    X = zeros(minFrameLen, M); % X matrix  entry of the beamforming
    for m = 1:M
        X(:,m) = frameSymCell{m}(1:minFrameLen);
    end
    rx.frameSymbolsMatrix = X;
end

end