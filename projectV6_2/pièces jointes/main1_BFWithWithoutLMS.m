clc; clear; close all;
rng(2026);

allResults = struct();

spacingCases = 1:2;   % 1 = lambda/2, 2 = 25 cm

% Channel mode
channelMode = 'direct';      % 'direct' or 'multipath'

% Plot controls
plotMainBerFigures = true;
plotSignalFigures  = false;
plotWeightFigures  = false;


% BER x SNR Doppler configuration

% Controls only the normal BER x Eb/N0 simulation.
% The Doppler residual sweep remains independent.
useDopplerInEbNoSweep = false;
compensateDopplerInEbNoSweep = false;

% Doppler velocity used in the normal BER x Eb/N0 simulation
% One value per hydrophone.
dopplerVelEbNo = [0.005 0.004 0.003 0.006];

% Steering vector configuration

% If true, the beamformer starts from an imperfect steering vector.
% This is useful to show that the LMS refines/converges from an initial error.
useImperfectSteering = true;

steeringErr = struct();
steeringErr.thetaOffsetDeg = 2;       % AoA error in degrees
steeringErr.phaseErrorStd_rad = 0.20; % random phase error per hydrophone
steeringErr.gainErrorStd_dB = 0.50;   % random gain error per hydrophone
steeringErr.positionErrorStd_m = 0;   % optional position error
steeringErr.seed = 1234;


% Doppler residual sweep configuration
doDopplerSweep = false;

% For d = 25 cm, only this branch keeps residual Doppler.
% The other branches are perfectly compensated.
dopplerBranch = 1;

% Compare these SNRs in the Doppler residual plot
EbNoDopplerList = [4 10];

% Residual velocity after Doppler compensation
dopplerResidualList = [0 ...
    0.001 0.002 0.003 0.004 0.005 ...
    0.006 0.007 0.008 0.009 0.010 ...
    0.011 0.012 0.013 0.014 0.015 ...
    0.0155 0.016 0.017 0.018 0.019 0.020];

% True physical Doppler before compensation
vTrueDoppler = max(dopplerResidualList);

numFramesDopplerSweep = 500;

% Loop over array spacings
for spacingCase = spacingCases

    %  Physical parameters
    env.c = 1500;              % sound speed [m/s]
    env.fc = 26e3;             % carrier frequency [Hz]
    env.fsADC = 156250;        % ADC sampling frequency [Hz]
    env.decim = 8;
    env.fs = env.fsADC / env.decim;
    env.lambda = env.c / env.fc;
 
    %  Array parameters 
    array.M = 4;
    %test
    %array.thetaDeg = 0;
    array.thetaDeg = 20;
    array.thetaRad = deg2rad(array.thetaDeg);

    switch spacingCase
        case 1
            array.d = env.lambda/2;
            spacingLabel = 'd_lambda_over_2';
        case 2
            array.d = 0.25;
            spacingLabel = 'd_25cm';
        otherwise
            error('Unknown spacing case.');
    end

    array.positions = (0:array.M-1) * array.d;
 
    %  Modulation, frame and receiver parameters
    p.M = 4;
    p.k = log2(p.M);
    p.sps = 8;
    p.Rs = env.fs / p.sps;
    p.phaseOffset = pi/4;

    p.trainingLengthSymbols = 1024;
    p.preambleLengthBits = p.trainingLengthSymbols * p.k;
    p.payloadLengthBits = 4000;

    p.EbNo_dB = 0:7;
    p.numFrames = 20;
    p.exampleSnrIndex = 4;

    p.useBestSamplingPhase = true ;

    p.useFractionalDelayCorrection = true;

    % Interpolation method used after parabolic TOA estimation.
    % Options: 'linear', 'pchip', 'spline'
    p.fracDelayInterpMethod = 'pchip';

    % Search window around the integer timing estimate.
    % This avoids locking on a wrong correlation peak far away.
    %p.fracDelaySearchHalfWindow = 2*p.sps;
    p.fracDelaySearchHalfWindow = 1;

    p.mu = 0.005;
    p.normalizeLMS = true;

   % Doppler control for the normal BER x Eb/N0 simulation
p.useDoppler = useDopplerInEbNoSweep;
p.compensateDoppler = compensateDopplerInEbNoSweep;
p.dopplerVel = dopplerVelEbNo;

if length(p.dopplerVel) ~= array.M
    error('p.dopplerVel must have one value per hydrophone.');
end
    p.channelMode = channelMode;

     
    %  Pulse shaping filters

    % Rectangular pulse shaping filter
    txFilter = ones(p.sps,1) / sqrt(p.sps);
    rxFilter = txFilter;

    % RRC filter :
    % rolloff = 0.25;
    % span = 8;
    % txFilter = rcosdesign(rolloff, span, p.sps, 'sqrt').';
    % rxFilter = txFilter;

    %  Channel definition 
    ch = buildChannel(channelMode, array, env, p);
 
    %  Preamble and training reference
    preambleBits = randi([0 1], p.preambleLengthBits, 1);
    preSym = qpskBitsToSymbols(preambleBits, p);
    dTrainSamples = buildPreambleSampleReference(preSym, txFilter, rxFilter, p);

    % Steering vector
aTrue = steeringVector(array.thetaDeg, array, env);
wTrue = aTrue / array.M;

if useImperfectSteering
    aInit = steeringVectorImperfect(array.thetaDeg, array, env, steeringErr);
    w0 = aInit / array.M;

    fprintf('Initial steering vector is imperfect.\n');
    fprintf('AoA offset = %.2f deg | phase std = %.2f rad | gain std = %.2f dB\n', ...
        steeringErr.thetaOffsetDeg, ...
        steeringErr.phaseErrorStd_rad, ...
        steeringErr.gainErrorStd_dB);
else
    w0 = wTrue;
    fprintf('Initial steering vector is perfect.\n');
end


    %  Theoretical curves 
    berTheorySISO = berawgn(p.EbNo_dB, 'psk', p.M, 'nondiff');
    berTheoryBF = berawgn(p.EbNo_dB + 10*log10(array.M), ...
        'psk', p.M, 'nondiff');
 
    fprintf('Spacing case %d | d = %.4f m = %.2f lambda | channel = %s\n', ...
        spacingCase, array.d, array.d/env.lambda, channelMode);
    fprintf('AoA = %.1f deg | sps = %d | frames = %d | mu = %.4f\n', ...
        array.thetaDeg, p.sps, p.numFrames, p.mu);

    %  BER vs Eb/N0
    [berBranches, berBFNoLMS, berBFLMS, example] = runEbNoSweep( ...
    env, array, ch, p, ...
    preambleBits, preSym, dTrainSamples, ...
    txFilter, rxFilter, w0, wTrue);


    %  Store normal BER results
    T = table(p.EbNo_dB(:), berTheorySISO(:), berTheoryBF(:), ...
        berBranches(:,1), berBranches(:,2), berBranches(:,3), berBranches(:,4), ...
        berBFNoLMS(:), berBFLMS(:), ...
        'VariableNames', {'EbNo_dB','Theory_SISO','Theory_BF', ...
        'Branch1','Branch2','Branch3','Branch4','BF_No_LMS','BF_LMS'});

    results = struct();
    results.env = env;
    results.array = array;
    results.p = p;
    results.ch = ch;
    results.berTheorySISO = berTheorySISO;
    results.berTheoryBF = berTheoryBF;
    results.berBranches = berBranches;
    results.berBFNoLMS = berBFNoLMS;
    results.berBFLMS = berBFLMS;
    results.table = T;
    results.example = example;

    allResults.(spacingLabel) = results;

    %  Plot normal BER results
    if plotMainBerFigures
        plotMainBerAndConvergence(p, array, env, channelMode, ...
            berTheorySISO, berTheoryBF, ...
            berBranches, berBFNoLMS, berBFLMS, example);
    end

    if plotSignalFigures
        plotSignalDiagnostics(example, p);
    end

    if plotWeightFigures
        plotWeightDiagnostics(example, array);
    end

    % Doppler residual sweep
if doDopplerSweep

    pSweep = p;
    pSweep.numFrames = numFramesDopplerSweep;

    switch spacingCase

        case 1
            % d = lambda/2:
            % normal array case.
            % All branches have the same residual Doppler after compensation.
            fprintf('\nRunning Doppler sweep for normal array | d = lambda/2\n');
            fprintf('True Doppler velocity before compensation = %.4f m/s\n', vTrueDoppler);
            fprintf('Residual Doppler applied to all branches\n');

            dopplerResults = runDopplerResidualSweep( ...
                env, array, ch, pSweep, ...
                preambleBits, preSym, dTrainSamples, ...
                txFilter, rxFilter, w0, ...
                EbNoDopplerList, dopplerResidualList, ...
                vTrueDoppler, dopplerBranch, ...
                'allBranches');

            dopplerResults.caseLabel = 'lambda_over_2_all_branches';
            allResults.dopplerLambdaOver2 = dopplerResults;

        case 2
            % d = 25 cm:
            fprintf('\nRunning Doppler sweep for d = 25 cm\n');
            fprintf('True Doppler velocity before compensation = %.4f m/s\n', vTrueDoppler);
  

            dopplerResults = runDopplerResidualSweep( ...
                env, array, ch, pSweep, ...
                preambleBits, preSym, dTrainSamples, ...
                txFilter, rxFilter, w0, ...
                EbNoDopplerList, dopplerResidualList, ...
                vTrueDoppler, dopplerBranch, ...
                'allBranches');

            dopplerResults.caseLabel = 'd_25cm';
            allResults.doppler25cm = dopplerResults;
    end
end
end

% Doppler Residual Plot

if doDopplerSweep && ...
        isfield(allResults, 'dopplerLambdaOver2') && ...
        isfield(allResults, 'doppler25cm')

    plotDopplerResidualResults( ...
        allResults.dopplerLambdaOver2, ...
        allResults.doppler25cm);
end

if doDopplerSweep && isfield(allResults, 'dopplerLambdaOver2')
    plotDopplerResidualMap(allResults.dopplerLambdaOver2, ...
        'Doppler residual map | d = \lambda/2');
end

if doDopplerSweep && isfield(allResults, 'doppler25cm')
    plotDopplerResidualMap(allResults.doppler25cm, ...
        'Doppler residual map | d = 25 cm');
end


save('allResults_BF_clean_FormB.mat', 'allResults');


% local functions

function ch = buildChannel(channelMode, array, env, p)

    switch channelMode

        case 'direct'
            ch.Npaths = 1;
            ch.delaysSec = 0;
            ch.gains = 1;
            ch.phases = 0;
            ch.thetaDeg = array.thetaDeg;
            ch.thetaRad = array.thetaRad;

        case 'multipath'
            geo.Zs = 5;
            geo.Zr = 20;
            geo.HauteurEau = 30;
            geo.DistanceHoriz = 200;
            geo.Celerite = env.c;

            v_toa = computeGeometricTOA(geo);

            ch.Npaths = 4;
            ch.toaAbs = v_toa(1:ch.Npaths);
            ch.delaysSec = ch.toaAbs - ch.toaAbs(1);

            %ch.gains = [1.0 0.6 0.25 0.1];
            %ch.gains = [1.0 1.0 1.0 1.0];

            ch.gains = [1.0 0.2 0.25 0.1];
            ch.gains = ch.gains / sqrt(sum(abs(ch.gains).^2));

            ch.phases = [0.00 0.60 -0.80 1.10];
            %angles that correspond to a better array response
             %ch.thetaDeg = [20 -31 59 -10];
             %intermediate
            %ch.thetaDeg = [20 -31 55 -10];
            %same
            %ch.thetaDeg = [20 20 20 20];
            %angles that correspond to a worst array response
             ch.thetaDeg = [20 -35 55 -10];
            ch.thetaRad = deg2rad(ch.thetaDeg);

        otherwise
            error('Unknown channel mode.');
    end

    if length(ch.gains) ~= ch.Npaths || ...
       length(ch.phases) ~= ch.Npaths || ...
       length(ch.thetaDeg) ~= ch.Npaths || ...
       length(ch.delaysSec) ~= ch.Npaths
        error('Channel vectors must have length equal to ch.Npaths.');
    end

    % Useful diagnostic stored in channel structure
    ch.delaySymbols = ch.delaysSec * p.Rs;
end

function [berBranches, berBFNoLMS, berBFLMS, example] = runEbNoSweep( ...
    env, array, ch, p, ...
    preambleBits, preSym, dTrainSamples, ...
    txFilter, rxFilter, w0, wTrue)

    nSNR = length(p.EbNo_dB);

    berBranches = zeros(nSNR, array.M);
    berBFNoLMS = zeros(nSNR, 1);
    berBFLMS = zeros(nSNR, 1);

    example = struct();

    for iSNR = 1:nSNR

        EbNo = p.EbNo_dB(iSNR);

        nErrBranch = zeros(1,array.M);
        nBitsBranch = zeros(1,array.M);
        nErrBFNoLMS = 0;
        nBitsBFNoLMS = 0;
        nErrBFLMS = 0;
        nBitsBFLMS = 0;

        for iFrame = 1:p.numFrames

           [frameResult, frameOK] = simulateOneFrameEbNo( ...
    env, array, ch, p, ...
    preambleBits, preSym, dTrainSamples, ...
    txFilter, rxFilter, w0, wTrue, EbNo);


            if ~frameOK
                continue;
            end

            nErrBranch = nErrBranch + frameResult.nErrBranch;
            nBitsBranch = nBitsBranch + frameResult.nBitsBranch;

            nErrBFNoLMS = nErrBFNoLMS + frameResult.nErrBFNoLMS;
            nBitsBFNoLMS = nBitsBFNoLMS + frameResult.nBitsBFNoLMS;

            nErrBFLMS = nErrBFLMS + frameResult.nErrBFLMS;
            nBitsBFLMS = nBitsBFLMS + frameResult.nBitsBFLMS;

            if iSNR == p.exampleSnrIndex && iFrame == 1
                example = frameResult.example;
                example.EbNo = EbNo;
                example.env = env;
                example.array = array;
                example.p = p;
                example.ch = ch;
            end
        end

        berBranches(iSNR,:) = nErrBranch ./ max(nBitsBranch,1);
        berBFNoLMS(iSNR) = nErrBFNoLMS / max(nBitsBFNoLMS,1);
        berBFLMS(iSNR) = nErrBFLMS / max(nBitsBFLMS,1);

        fprintf('Eb/N0 = %2d dB | BF no LMS = %.3e | BF LMS = %.3e\n', ...
            EbNo, berBFNoLMS(iSNR), berBFLMS(iSNR));
    end
end

function [frameResult, frameOK] = simulateOneFrameEbNo( ...
    env, array, ch, p, ...
    preambleBits, preSym, dTrainSamples, ...
    txFilter, rxFilter, w0, wTrue, EbNo_dB)

    dataBits = randi([0 1], p.payloadLengthBits, 1);

    [~, ~, txWave] = txChain(dataBits, preambleBits, p, txFilter);

    chanOut = applyChannelSIMO(txWave, ch, array, env);

    % Optional Doppler generation for BER x Eb/N0 simulation
    if p.useDoppler
        for m = 1:array.M
            [dopOut, ~] = f_dopplergeneration( ...
                chanOut(:,m), ...
                p.dopplerVel(m), ...
                env.fs, env.c, env.fc);

            chanOut(:,m) = matchLength(dopOut(:), size(chanOut,1));
        end
    end

    snrSample_dB = EbNo_dB + 10*log10(p.k) - 10*log10(p.sps);
    rxWave = addAwgnPerBranch(chanOut, snrSample_dB, array.M);

    % Optional Doppler compensation for BER x Eb/N0 simulation
    if p.useDoppler && p.compensateDoppler
        for m = 1:array.M
            compOut = f_dopplercompensation( ...
                rxWave(:,m), ...
                p.dopplerVel(m), ...
                env.fs, env.c, env.fc);

            rxWave(:,m) = matchLength(compOut(:), size(rxWave,1));
        end
    end

    [frameResult, frameOK] = receiverAndBerFromRxWave( ...
        rxWave, dataBits, array, p, ...
        preambleBits, preSym, dTrainSamples, ...
        rxFilter, w0);

    if frameOK
        frameResult.example = buildExampleStruct(frameResult, w0, wTrue, dTrainSamples);
    else
        frameResult.example = struct();
    end
end


function dopplerResults = runDopplerResidualSweep( ...
    env, array, ch, p, ...
    preambleBits, preSym, dTrainSamples, ...
    txFilter, rxFilter, w0, ...
    EbNoList, residualList, vTrueDoppler, dopplerBranch, dopplerMode)

    nVel = length(residualList);
    nEb = length(EbNoList);

    berArray = NaN(nVel, nEb);
    nBitsArray = NaN(nVel, nEb);

    for iEb = 1:nEb

        EbNoThis = EbNoList(iEb);
        fprintf('  Eb/N0 = %.1f dB\n', EbNoThis);

        for iVel = 1:nVel

            vResidual = residualList(iVel);
            rng(200000 + 1000*iEb + iVel);

            [berArray(iVel,iEb), nBitsArray(iVel,iEb)] = ...
                runOneBERPointDoppler( ...
                    env, array, ch, p, ...
                    preambleBits, preSym, dTrainSamples, ...
                    txFilter, rxFilter, w0, ...
                    EbNoThis, vTrueDoppler, vResidual, ...
                    dopplerBranch, dopplerMode);

            fprintf('    residual = %.4f m/s | BF+LMS = %.3e\n', ...
                vResidual, berArray(iVel,iEb));
        end
    end

    dopplerResults.velocityList = residualList;
    dopplerResults.dopplerFrequencyList = residualList * env.fc / env.c;
    dopplerResults.EbNoList = EbNoList;
    dopplerResults.berArray = berArray;
    dopplerResults.nBitsArray = nBitsArray;
    dopplerResults.dopplerBranch = dopplerBranch;
    dopplerResults.dopplerMode = dopplerMode;
    dopplerResults.vTrueDoppler = vTrueDoppler;
    dopplerResults.arraySpacing_m = array.d;
    dopplerResults.arraySpacing_lambda = array.d/env.lambda;
end

function [berArray, nBitsArray] = ...
    runOneBERPointDoppler( ...
        env, array, ch, p, ...
        preambleBits, preSym, dTrainSamples, ...
        txFilter, rxFilter, w0, ...
        EbNo_dB, vTrueDoppler, vResidual, ...
        dopplerBranch, dopplerMode)

    nErrArray = 0;
    nBitsArray = 0;

    snrSample_dB = EbNo_dB + 10*log10(p.k) - 10*log10(p.sps);

    for iFrame = 1:p.numFrames

        dataBits = randi([0 1], p.payloadLengthBits, 1);

        [~, ~, txWave] = txChain(dataBits, preambleBits, p, txFilter);

        chanOut = applyChannelSIMO(txWave, ch, array, env);

        chanOut = applyDoppler( ...
            chanOut, env, array.M, ...
            vTrueDoppler, vResidual, ...
            dopplerBranch, dopplerMode);

        rxWave = addAwgnPerBranch(chanOut, snrSample_dB, array.M);

        [frameResult, frameOK] = receiverAndBerFromRxWave( ...
            rxWave, dataBits, array, p, ...
            preambleBits, preSym, dTrainSamples, ...
            rxFilter, w0);

        if ~frameOK
            continue;
        end

        nErrArray = nErrArray + frameResult.nErrBFLMS;
        nBitsArray = nBitsArray + frameResult.nBitsBFLMS;
    end

    berArray = nErrArray / max(nBitsArray,1);
end


function chanOut = applyDoppler( ...
    chanOut, env, M, vTrueDoppler, vResidual, ...
    dopplerBranch, dopplerMode)

    N = size(chanOut,1);

    for m = 1:M

        % True physical Doppler generated on every branch
        [dopOut, ~] = f_dopplergeneration( ...
            chanOut(:,m), vTrueDoppler, env.fs, env.c, env.fc);

        dopOut = matchLength(dopOut(:), N);

        switch dopplerMode

            case 'allBranches'
                % Normal array case:
                % every branch has the same residual Doppler after compensation.
                vEstimated = vTrueDoppler - vResidual;

            case 'oneBranch'
                % One-branch residual case:
                % only dopplerBranch keeps residual Doppler.
                % All other branches are perfectly compensated.
                if m == dopplerBranch
                    vEstimated = vTrueDoppler - vResidual;
                else
                    vEstimated = vTrueDoppler;
                end

            otherwise
                error('Unknown dopplerMode.');
        end

        % Doppler compensation
        compOut = f_dopplercompensation( ...
            dopOut, vEstimated, env.fs, env.c, env.fc);

        chanOut(:,m) = matchLength(compOut(:), N);
    end
end

function [frameResult, frameOK] = receiverAndBerFromRxWave( ...
    rxWave, dataBits, array, p, ...
    preambleBits, preSym, dTrainSamples, ...
    rxFilter, w0)

    preLenSym = length(preSym);
    trainSampleLen = length(dTrainSamples);

    numFrameSym = preLenSym + p.payloadLengthBits/p.k;
    numFrameSamples = numFrameSym * p.sps;

    rxBranch = rxFrontEndSIMO(rxWave, preambleBits, p, rxFilter);

    [Xsamples, frameOK,timingInfo] = buildAlignedSampleMatrix(rxBranch, numFrameSamples, array, p, dTrainSamples);

    frameResult = struct();

    if ~frameOK
        return;
    end

    % Branch BER
    nErrBranch = zeros(1,array.M);
    nBitsBranch = zeros(1,array.M);

    for m = 1:array.M
        frameSym = Xsamples(1:p.sps:end,m);
        frameSym = frameSym(1:numFrameSym);
        frameSym = phaseGainNormalize(frameSym, preSym);

        payloadSym = frameSym(preLenSym+1:end);
        rxBitsBranch = qpskSymbolsToBits(payloadSym, p);

        minLen = min(length(dataBits), length(rxBitsBranch));

        nErrBranch(m) = sum(dataBits(1:minLen) ~= rxBitsBranch(1:minLen));
        nBitsBranch(m) = minLen;
    end

    %% BF without LMS
    yBFNoLMS_samples = Xsamples * conj(w0);
    yBFNoLMS = yBFNoLMS_samples(1:p.sps:end);
    yBFNoLMS = yBFNoLMS(1:numFrameSym);
    yBFNoLMS = phaseGainNormalize(yBFNoLMS, preSym);

    payloadNoLMS = yBFNoLMS(preLenSym+1:end);
    rxBitsNoLMS = qpskSymbolsToBits(payloadNoLMS, p);

    minLenNoLMS = min(length(dataBits), length(rxBitsNoLMS));
    nErrBFNoLMS = sum(dataBits(1:minLenNoLMS) ~= rxBitsNoLMS(1:minLenNoLMS));
    nBitsBFNoLMS = minLenNoLMS;

    % BF with sample-domain LMS
    XtrainSamples = Xsamples(1:trainSampleLen,:);

    [wLMS, eTrain] = spatialLMSRefineSteeringSamples( ...
        XtrainSamples, dTrainSamples, p, w0);

    yBFLMS_samples = Xsamples * conj(wLMS);
    yBFLMS = yBFLMS_samples(1:p.sps:end);
    yBFLMS = yBFLMS(1:numFrameSym);
    yBFLMS = phaseGainNormalize(yBFLMS, preSym);

    payloadLMS = yBFLMS(preLenSym+1:end);
    rxBitsLMS = qpskSymbolsToBits(payloadLMS, p);

    minLenLMS = min(length(dataBits), length(rxBitsLMS));
    nErrBFLMS = sum(dataBits(1:minLenLMS) ~= rxBitsLMS(1:minLenLMS));
    nBitsBFLMS = minLenLMS;

    % Output
    frameResult.nErrBranch = nErrBranch;
    frameResult.nBitsBranch = nBitsBranch;

    frameResult.nErrBFNoLMS = nErrBFNoLMS;
    frameResult.nBitsBFNoLMS = nBitsBFNoLMS;

    frameResult.nErrBFLMS = nErrBFLMS;
    frameResult.nBitsBFLMS = nBitsBFLMS;

    frameResult.rxBranch = rxBranch;
    frameResult.Xsamples = Xsamples;
    frameResult.yBFNoLMS = yBFNoLMS;
    frameResult.yBFLMS = yBFLMS;
    frameResult.yBFNoLMS_samples = yBFNoLMS_samples;
    frameResult.yBFLMS_samples = yBFLMS_samples;
    frameResult.wLMS = wLMS;
    frameResult.eTrain = eTrain;
    frameResult.preLenSym = preLenSym;
    frameResult.numFrameSym = numFrameSym;
    frameResult.timingInfo = timingInfo;
end

function [Xsamples, frameOK, timingInfo] = buildAlignedSampleMatrix( ...
    rxBranch, numFrameSamples, array, p, dTrainSamples)

    Xsamples = zeros(numFrameSamples, array.M);
    frameOK = true;

    timingInfo = struct();
    timingInfo.integerStart = NaN(1,array.M);
    timingInfo.corrPeakIndex = NaN(1,array.M);
    timingInfo.fracOffset = NaN(1,array.M);
    timingInfo.fractionalStart = NaN(1,array.M);
    timingInfo.peakAmplitude = NaN(1,array.M);

    useFrac = isfield(p, 'useFractionalDelayCorrection') && ...
              p.useFractionalDelayCorrection;

    if isfield(p, 'fracDelayInterpMethod')
        interpMethod = p.fracDelayInterpMethod;
    else
        interpMethod = 'pchip';
    end

    if isfield(p, 'fracDelaySearchHalfWindow')
        searchHalfWindow = p.fracDelaySearchHalfWindow;
    else
        searchHalfWindow = 2*p.sps;
    end

    ref = dTrainSamples(:);
    refLen = length(ref);

    for m = 1:array.M

        rxMfNoTransient = rxBranch.rxMfNoTransientCell{m};
        rxMfNoTransient = rxMfNoTransient(:);

        if isempty(rxMfNoTransient)
            frameOK = false;
            return;
        end

        % Original integer timing estimate from rxFrontEndSIMO
        integerStart = rxBranch.bestPhase(m) + ...
            (rxBranch.frameStartSym(m)-1) * p.sps;
        timingInfo.integerStart(m) = integerStart;
        
        % If the frame starts too close to the beginning, parabolic interpolation
% is unreliable because the correlation peak is at the boundary.
if integerStart <= 2
    sampleStart = integerStart;
    lastSample = sampleStart + numFrameSamples - 1;

    if sampleStart < 1 || lastSample > length(rxMfNoTransient)
        frameOK = false;
        return;
    end

    Xsamples(:,m) = rxMfNoTransient(sampleStart:lastSample);

    timingInfo.corrPeakIndex(m) = sampleStart;
    timingInfo.fracOffset(m) = 0;
    timingInfo.fractionalStart(m) = sampleStart;
    timingInfo.peakAmplitude(m) = NaN;

    continue;
end


        if ~useFrac

            % Original behavior: integer extraction only
            sampleStart = integerStart;
            lastSample = sampleStart + numFrameSamples - 1;

            if sampleStart < 1 || lastSample > length(rxMfNoTransient)
                frameOK = false;
                return;
            end

            Xsamples(:,m) = rxMfNoTransient(sampleStart:lastSample);

            timingInfo.corrPeakIndex(m) = sampleStart;
            timingInfo.fracOffset(m) = 0;
            timingInfo.fractionalStart(m) = sampleStart;
            timingInfo.peakAmplitude(m) = NaN;

        else

            % Correlation between received matched-filtered signal
            % and known sample-domain preamble reference.
            if length(rxMfNoTransient) < refLen
                frameOK = false;
                return;
            end

            Rxx = conv(rxMfNoTransient, flipud(conj(ref)), 'valid');

            % Search around the previous integer timing estimate.
            % This avoids locking on a wrong correlation peak.
            idxMin = max(2, round(integerStart) - searchHalfWindow);
            idxMax = min(length(Rxx)-1, round(integerStart) + searchHalfWindow);

            if idxMin >= idxMax
                frameOK = false;
                return;
            end

            % [~, localIdx] = max(abs(Rxx(idxMin:idxMax)));
            % peakIdx = idxMin + localIdx - 1;
            % 
            % % Use your parabolic interpolation function.
            % % Fs = 1 here because we want the TOA in sample-index units,
            % % not in seconds.
            % [TOA_samples, AMP, fracOffset, ~] = P_InterpParaboTiming(Rxx, peakIdx);
            % 
            % if ~isfinite(TOA_samples) || ~isfinite(fracOffset) || abs(fracOffset) > 0.45
            % TOA_samples = peakIdx;
            %     fracOffset = 0;
            % AMP = abs(Rxx(peakIdx));
            % end
            % 
            % fractionalStart = TOA_samples;

            % Conservative fractional correction:
% use the integer timing from rxFrontEndSIMO as the reference.
% The parabolic interpolation must only refine this point locally.

peakIdx = round(integerStart);

% If peakIdx is too close to the correlation boundary, do not interpolate.
if peakIdx <= 1 || peakIdx >= length(Rxx)
    fracOffset = 0;
    AMP = abs(Rxx(max(1,min(length(Rxx),peakIdx))));
    fractionalStart = integerStart;
else
    [~, AMP, fracOffset, ~] = P_InterpParaboTiming(Rxx, peakIdx);

    % Safety: reject large or invalid corrections.
    if ~isfinite(fracOffset) || abs(fracOffset) > 0.35
        fracOffset = 0;
        AMP = abs(Rxx(peakIdx));
    end

    fractionalStart = integerStart + fracOffset;
end

            timingInfo.corrPeakIndex(m) = peakIdx;
            timingInfo.fracOffset(m) = fracOffset;
            timingInfo.fractionalStart(m) = fractionalStart;
            timingInfo.peakAmplitude(m) = AMP;

            % Extract the frame at fractional sample positions.
            sampleGrid = fractionalStart + (0:numFrameSamples-1).';
            n = (1:length(rxMfNoTransient)).';

            if sampleGrid(1) < 1 || sampleGrid(end) > length(rxMfNoTransient)
                frameOK = false;
                return;
            end

            Xsamples(:,m) = interp1( ...
                n, rxMfNoTransient, sampleGrid, interpMethod, 0);
        end
    end
end


% function rxWave = addAwgnPerBranch(chanOut, snrSample_dB, M)
% 
%     rxWave = zeros(size(chanOut));
% 
%     for m = 1:M
%         rxWave(:,m) = awgn(chanOut(:,m), snrSample_dB, 'measured');
%     end
% end

function rxWave = addAwgnPerBranch(chanOut, snrSample_dB, M)

    % Common AWGN applied equally to all hydrophones.
    % This is only for testing spatially correlated/common noise.

    rxWave = zeros(size(chanOut));

    % Reference signal power measured over all branches
    sigPower = mean(abs(chanOut(:)).^2);

    % Convert SNR from dB to linear
    snrLinear = 10^(snrSample_dB/10);

    % Required noise power
    noisePower = sigPower / snrLinear;

    % Generate one common complex noise vector
    commonNoise = sqrt(noisePower/2) * ...
        (randn(size(chanOut,1),1) + 1j*randn(size(chanOut,1),1));

    % Add the same noise to all branches
    for m = 1:M
        rxWave(:,m) = chanOut(:,m) + commonNoise;
    end
end

function example = buildExampleStruct(frameResult, w0, wTrue, dTrainSamples)

    example = struct();
    example.rxBranch = frameResult.rxBranch;
    example.yBFNoLMS = frameResult.yBFNoLMS;
    example.yBFLMS = frameResult.yBFLMS;
    example.yBFNoLMS_samples = frameResult.yBFNoLMS_samples;
    example.yBFLMS_samples = frameResult.yBFLMS_samples;

    example.w0 = w0;
    example.wTrue = wTrue;
    example.wLMS = frameResult.wLMS;

    example.eTrain = frameResult.eTrain;
    example.dTrainSamples = dTrainSamples;

    example.preLenSym = frameResult.preLenSym;
    example.numFrameSym = frameResult.numFrameSym;
    example.timingInfo = frameResult.timingInfo;

    example.Xsamples = frameResult.Xsamples;
    example.rxBranch = frameResult.rxBranch;
end


function plotMainBerAndConvergence(p, array, env, channelMode, ...
    berTheorySISO, berTheoryBF, ...
    berBranches, berBFNoLMS, berBFLMS, example)

    figure('Name', sprintf('BER and LMS convergence | d = %.2f lambda | %s', ...
        array.d/env.lambda, channelMode));

    subplot(1,2,1);
    semilogy(p.EbNo_dB, berTheorySISO, 'k-', 'LineWidth', 1.6); hold on;
    semilogy(p.EbNo_dB, berTheoryBF, 'g-', 'LineWidth', 1.6);

    for m = 1:array.M
        semilogy(p.EbNo_dB, berBranches(:,m), '--', 'LineWidth', 1.1);
    end

    semilogy(p.EbNo_dB, berBFNoLMS, 'bo-', 'LineWidth', 1.8);
    semilogy(p.EbNo_dB, berBFLMS, 'ro-', 'LineWidth', 1.8);

    grid on;
    xlabel('E_b/N_0 (dB)');
    ylabel('BER');
    title(sprintf('BER | d = %.2f \\lambda | %s', ...
        array.d/env.lambda, channelMode));

    legend('Theory SISO', 'Theory BF + array gain', ...
        'Branch 1', 'Branch 2', 'Branch 3', 'Branch 4', ...
        'BF without LMS', 'BF with sample-LMS', ...
        'Location', 'southwest');

    subplot(1,2,2);
    plot(abs(example.eTrain), 'LineWidth', 1.2); hold on;
    plot(movmean(abs(example.eTrain), 64), 'LineWidth', 1.8);
    grid on;
    xlabel('Training sample index');
    ylabel('|e[n]|');
    title(sprintf('Sample-time spatial LMS convergence | d = %.2f \\lambda', ...
        array.d/env.lambda));
    legend('|e[n]|', 'Moving average', 'Location', 'northeast');
end

function plotDopplerResidualResults(dopplerLambdaOver2, doppler25cm)

    vList = dopplerLambdaOver2.velocityList;
    EbList = dopplerLambdaOver2.EbNoList;

    berLambda = dopplerLambdaOver2.berArray;
    ber25cm = doppler25cm.berArray;

    nBitsLambda = dopplerLambdaOver2.nBitsArray;
    nBits25cm = doppler25cm.nBitsArray;

    % Visual floor for zero-BER points in semilogy
    berLambdaPlot = berLambda;
    ber25cmPlot = ber25cm;

    for iEb = 1:length(EbList)
        floorLambda = 0.5 / max(nBitsLambda(1,iEb),1);
        floor25cm = 0.5 / max(nBits25cm(1,iEb),1);

        berLambdaPlot(:,iEb) = max(berLambda(:,iEb), floorLambda);
        ber25cmPlot(:,iEb) = max(ber25cm(:,iEb), floor25cm);
    end

    figure('Name','BER versus residual Doppler');

    semilogy(vList, berLambdaPlot(:,1), 'o-', 'LineWidth', 1.8); hold on;
    semilogy(vList, berLambdaPlot(:,2), 'o--', 'LineWidth', 1.8);

    semilogy(vList, ber25cmPlot(:,1), 's-', 'LineWidth', 1.8);
    semilogy(vList, ber25cmPlot(:,2), 's--', 'LineWidth', 1.8);

    grid on;
    xlabel('Residual Doppler velocity after compensation (m/s)');
    ylabel('BER');

    title(sprintf(['BER vs residual Doppler after compensation\n' ...
        '\\lambda/2: residual on full array']));

    legend( ...
        sprintf('d = \\lambda/2 | array normal | E_b/N_0 = %.0f dB', EbList(1)), ...
        sprintf('d = \\lambda/2 | array normal | E_b/N_0 = %.0f dB', EbList(2)), ...
        sprintf('d = 25 cm |  E_b/N_0 = %.0f dB', ...
            EbList(1)), ...
        sprintf('d = 25 cm | E_b/N_0 = %.0f dB', ...
             EbList(2)), ...
        'Location', 'best');
end

function y = matchLength(x, N)

    x = x(:);

    if length(x) >= N
        y = x(1:N);
    else
        y = [x; zeros(N-length(x),1)];
    end
end

function plotSignalDiagnostics(example, p)

    figure('Name','Signal diagnostics');

    payloadNoLMS = example.yBFNoLMS(example.preLenSym+1:end);
    payloadLMS = example.yBFLMS(example.preLenSym+1:end);

    subplot(2,2,1);
    Nplot = min(3000, length(payloadNoLMS));
    scatter(real(payloadNoLMS(1:Nplot)), imag(payloadNoLMS(1:Nplot)), '.');
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('Constellation - BF without LMS | Eb/N0 = %d dB', example.EbNo));

    subplot(2,2,2);
    Nplot = min(3000, length(payloadLMS));
    scatter(real(payloadLMS(1:Nplot)), imag(payloadLMS(1:Nplot)), '.');
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('Constellation - BF with sample-LMS | Eb/N0 = %d dB', example.EbNo));

    subplot(2,2,3);
    plotEyeOnAxis(gca, real(example.yBFNoLMS), p.sps, 'Eye - BF without LMS');

    subplot(2,2,4);
    plotEyeOnAxis(gca, real(example.yBFLMS), p.sps, 'Eye - BF with sample-LMS');
end

function plotWeightDiagnostics(example, array)

    figure('Name','Beamforming coefficients');

    subplot(2,2,1);
    stem(1:array.M, abs(example.wTrue), 'k', 'filled', 'LineWidth', 1.4); hold on;
    stem(1:array.M, abs(example.w0), 'b', 'filled', 'LineWidth', 1.4);
    stem(1:array.M, abs(example.wLMS), 'r', 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('|w_m|');
    title('Beamforming coefficient magnitude');
    legend('True weights', 'Initial weights', 'LMS refined weights', ...
        'Location', 'best');

    subplot(2,2,2);
    stem(1:array.M, unwrap(angle(example.wTrue)), 'k', 'filled', 'LineWidth', 1.4); hold on;
    stem(1:array.M, unwrap(angle(example.w0)), 'b', 'filled', 'LineWidth', 1.4);
    stem(1:array.M, unwrap(angle(example.wLMS)), 'r', 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('Phase of w_m (rad)');
    title('Beamforming coefficient phase');
    legend('True weights', 'Initial weights', 'LMS refined weights', ...
        'Location', 'best');

    subplot(2,2,3);
    stem(1:array.M, abs(example.w0 - example.wTrue), 'b', 'filled', 'LineWidth', 1.4); hold on;
    stem(1:array.M, abs(example.wLMS - example.wTrue), 'r', 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('Error relative to true weights');
    title('Initial error vs LMS residual error');
    legend('|w_0 - w_{true}|', '|w_{LMS} - w_{true}|', ...
        'Location', 'best');

    subplot(2,2,4);
    plot(abs(example.eTrain), 'LineWidth', 1.2); hold on;
    plot(movmean(abs(example.eTrain), 64), 'LineWidth', 1.8);
    grid on;
    xlabel('Training sample index');
    ylabel('|e[n]|');
    title('LMS convergence during preamble');
    legend('|e[n]|', 'Moving average', 'Location', 'northeast');
end

function plotEyeOnAxis(ax, sigSym, sps, ttl)

    sigSym = sigSym(:);
    sigEye = interp(sigSym, sps);

    nTraces = 80;
    L = 2*sps;

    cla(ax);
    hold(ax, 'on');
    grid(ax, 'on');

    idx = 1;
    count = 0;

    while idx + L - 1 <= length(sigEye) && count < nTraces
        seg = sigEye(idx:idx+L-1);
        plot(ax, linspace(-0.5,0.5,L), seg, 'LineWidth', 0.8);
        idx = idx + sps;
        count = count + 1;
    end

    xlabel(ax, 'Time');
    ylabel(ax, 'Amplitude');
    title(ax, ttl);
end

function plotDopplerResidualMap(dopplerResults, figTitle)

    vList = dopplerResults.velocityList(:).';   % X axis
    EbList = dopplerResults.EbNoList(:);        % Y axis

    % berArray is stored as:
    % rows = velocity points
    % cols = Eb/N0 points
    %
    % For surf/imagesc, it is more convenient to use:
    % rows = Eb/N0
    % cols = velocity
    berMat = dopplerResults.berArray.';         % size = nEb x nVel

    % Visual floor for zero-BER points
    berPlot = berMat;

    for iEb = 1:length(EbList)
        floorBer = 0.5 / max(dopplerResults.nBitsArray(1,iEb), 1);
        berPlot(iEb,:) = max(berPlot(iEb,:), floorBer);
    end

    % Because BER spans many decades, plot log10(BER)
    Z = log10(berPlot);

    [X, Y] = meshgrid(vList, EbList);

    %% 3D surface plot
    figure('Name', [figTitle ' | surf']);
    surf(X, Y, Z);
    shading interp;
    grid on;
    colorbar;
    xlabel('Residual Doppler speed (m/s)');
    ylabel('E_b/N_0 (dB)');
    zlabel('log_{10}(BER)');
    title([figTitle ' | 3D surface']);

end