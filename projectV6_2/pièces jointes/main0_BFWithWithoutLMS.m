clc; clear; close all;
rng(2026);

allResults = struct();

spacingCases = 1:2;

dopplerVel = [0.005 0.004 0.003 0.006];

%% Doppler sweep configuration - one residual Doppler branch
doDopplerSweep = true;

% Only this branch has residual Doppler.
% Other branches are assumed perfectly Doppler compensated.
dopplerBranch = 1;

% Fixed Eb/N0 values for comparison
EbNoDopplerList = [4 10];

% Reference Doppler velocity applied only to dopplerBranch
dopplerNominalList = [0 ...
    0.001 0.002 0.003 0.004 0.005 ...
    0.006 0.007 0.008 0.009 0.010 ...
    0.011 0.012 0.013 0.014 0.015 ...
    0.0155 0.016 0.017 0.018 0.019 0.020];

% Increase if the high-SNR curve is still noisy
numFramesDopplerSweep = 100;

% Automatic knee detection parameters
kneeFactor = 3;      % curve is considered to rise when BER > 3 x baseline
minErrForKnee = 5;   % avoid detecting knee from only 1 random error

for spacingCase = spacingCases

    %% Physical parameters
    env.c = 1500;
    env.fc = 26e3;

    env.fsADC = 156250;
    env.decim = 8;
    env.fs = env.fsADC / env.decim;

    env.lambda = env.c / env.fc;

    %% Array parameters
    array.M = 4;

        %test
         %array.thetaDeg = 0;
    array.thetaDeg = 20;

    array.thetaRad = deg2rad(array.thetaDeg);

    switch spacingCase
        case 1
            array.d = env.lambda/2;
            spacingLabel = 'd_lambda_over_2';
            spacingLegend = '\lambda/2';

        case 2
            array.d = 0.25;
            spacingLabel = 'd_25cm';
            spacingLegend = '25 cm';
    end

    %array.positions = zeros(1,array.M);
    array.positions = (0:array.M-1) * array.d;

    %% Modulation parameters
    p.M = 4;
    p.k = log2(p.M);
    p.sps = 8;
    p.Rs = env.fs / p.sps;
    p.phaseOffset = pi/4;

    %test
    % rolloff = 0.25;
    % span = 8;
    %txFilter = rcosdesign(rolloff, span, p.sps,'sqrt');
    txFilter = ones(p.sps,1) / sqrt(p.sps);
    rxFilter = txFilter;

    %% Frame parameters
    p.trainingLengthSymbols = 1024;
    p.preambleLengthBits = p.trainingLengthSymbols * p.k;
    p.payloadLengthBits = 4000;

    %% Monte Carlo parameters
    p.EbNo_dB = 0:7;
    p.numFrames = 200;
    p.exampleSnrIndex = 4;

    %% Receiver and LMS parameters
    p.useBestSamplingPhase = true;
    p.mu = 0.05;
    p.normalizeLMS = true;

    %% Channel mode
    %p.channelMode = 'multipath';
    p.channelMode = 'direct';


    %% Doppler parameters
    p.useDoppler = false;
    p.compensateDoppler = false;

    if length(dopplerVel) ~= array.M
        error('dopplerVel must have one value per hydrophone.');
    end

    %% Channel definition
    switch p.channelMode

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
            ch.gains = [1.0 0.20 0.25 0.1];


            ch.gains = ch.gains / sqrt(sum(abs(ch.gains).^2));

            ch.phases = [0.00 0.60 -0.80 1.10];

            ch.thetaDeg = [20 -35 55 -10];
            %angles that correspond to a weak gain in the array response
             %ch.thetaDeg = [20 -31 59 -10];
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

    delaySym = ch.delaysSec * p.Rs;

    disp(table((1:ch.Npaths).', ch.delaysSec(:), delaySym(:), ...
        ch.gains(:), ch.phases(:), ch.thetaDeg(:), ...
        'VariableNames', {'Path','Delay_s','Delay_symbols','Gain','Phase_rad','AoA_deg'}));

    %% Preamble
    preambleBits = randi([0 1], p.preambleLengthBits, 1);
    preSym = qpskBitsToSymbols(preambleBits, p);
    preLenSym = length(preSym);

    dTrainSamples = buildPreambleSampleReference(preSym, txFilter, rxFilter, p);
    trainSampleLen = length(dTrainSamples);

    %% Theoretical BER curves
    berTheorySISO = berawgn(p.EbNo_dB, 'psk', p.M, 'nondiff');

    berTheoryBF = berawgn(p.EbNo_dB + 10*log10(array.M), ...
        'psk', p.M, 'nondiff');

    %% Storage
    berBranches = zeros(length(p.EbNo_dB), array.M);
    berBFNoLMS = zeros(size(p.EbNo_dB));
    berBFLMS = zeros(size(p.EbNo_dB));

    example = struct();

    fprintf('\nBFWithWithoutLMS | spacing case %d | channel = %s\n', ...
        spacingCase, p.channelMode);
    fprintf('d = %.4f m = %.2f lambda\n', array.d, array.d/env.lambda);
    fprintf('AoA = %.1f deg\n', array.thetaDeg);
    fprintf('Training = %d symbols\n', p.trainingLengthSymbols);
    fprintf('LMS mu = %.4f\n', p.mu);
    fprintf('Beamforming LMS applied at sample time\n');
    fprintf('Doppler enabled = %d | known-velocity compensation enabled = %d\n', ...
        p.useDoppler, p.compensateDoppler);


     %% Steering vector selection
useImperfectSteering = false;

    if useImperfectSteering
       aTrue = steeringVector(array.thetaDeg, array, env);

    err = struct();
    err.thetaOffsetDeg = 2;
    err.phaseErrorStd_rad = 0.2;
    err.gainErrorStd_dB = 0;
    err.positionErrorStd_m = 0;
    err.seed = 1234;

    aInit = steeringVectorImperfect(array.thetaDeg, array, env, err);

    wTrue = aTrue / array.M;
    w0 = aInit / array.M;
    else
    a0 = steeringVector(array.thetaDeg, array, env);
    w0 = a0 / array.M;
    end



    %% BER versus Eb/N0
    for iSNR = 1:length(p.EbNo_dB)

        EbNo = p.EbNo_dB(iSNR);

        nErrBranch = zeros(1,array.M);
        nBitsBranch = zeros(1,array.M);

        nErrBFNoLMS = 0;
        nBitsBFNoLMS = 0;

        nErrBFLMS = 0;
        nBitsBFLMS = 0;

        for iFrm = 1:p.numFrames
%transmission
            dataBits = randi([0 1], p.payloadLengthBits, 1);

            [~, ~, txWave] = txChain(dataBits, preambleBits, p, txFilter);


%apply the channel
            chanOut = applyChannelSIMO(txWave, ch, array, env);

       %test: verify if the hydrphones receive the same signal
            %chanOut = repmat(txWave, 1, array.M);
        %test
%         aTest = steeringVector(array.thetaDeg, array, env);
% 
% chanOut = zeros(length(txWave), array.M);
% 
% for m = 1:array.M
%     chanOut(:,m) = aTest(m) * txWave;
% end
%only Doppler
            if p.useDoppler
                for m = 1:array.M
                    [dopOut, ~] = f_dopplergeneration( ...
                        chanOut(:,m), dopplerVel(m), ...
                        env.fs, env.c, env.fc);

                    chanOut(:,m) = matchLength(dopOut(:), size(chanOut,1));
                end
            end

            snrSample_dB = EbNo + 10*log10(p.k) - 10*log10(p.sps);

            rxWave = zeros(size(chanOut));

            for m = 1:array.M
                rxWave(:,m) = awgn(chanOut(:,m), snrSample_dB, 'measured');
            end
%Doppler + compensation
            if p.useDoppler && p.compensateDoppler
                for m = 1:array.M
                    compOut = f_dopplercompensation( ...
                        rxWave(:,m), dopplerVel(m), ...
                        env.fs, env.c, env.fc);

                    rxWave(:,m) = matchLength(compOut(:), size(rxWave,1));
                end
            end
%reception
            rxBranch = rxFrontEndSIMO(rxWave, preambleBits, p, rxFilter);

            numFrameSym = preLenSym + p.payloadLengthBits/p.k;
            numFrameSamples = numFrameSym * p.sps;

            Xsamples = zeros(numFrameSamples, array.M);
            validSamples = true;

            for m = 1:array.M

                rxMfNoTransient = rxBranch.rxMfNoTransientCell{m};

                sampleStart = rxBranch.bestPhase(m) + ...
                    (rxBranch.frameStartSym(m)-1) * p.sps;

                lastSample = sampleStart + numFrameSamples - 1;

                if isempty(rxMfNoTransient) || lastSample > length(rxMfNoTransient)
                    validSamples = false;
                    break;
                end

                Xsamples(:,m) = rxMfNoTransient(sampleStart:lastSample);
            end

            if ~validSamples
                continue;
            end

            for m = 1:array.M

                frameSym = Xsamples(1:p.sps:end,m);
                frameSym = frameSym(1:numFrameSym);

                frameSym = phaseGainNormalize(frameSym, preSym);

                payloadSym = frameSym(preLenSym+1:end);
                rxBitsBranch = qpskSymbolsToBits(payloadSym, p);

                minLen = min(length(dataBits), length(rxBitsBranch));

                nErrBranch(m) = nErrBranch(m) + ...
                    sum(dataBits(1:minLen) ~= rxBitsBranch(1:minLen));

                nBitsBranch(m) = nBitsBranch(m) + minLen;
            end

            yBFNoLMS_samples = Xsamples * conj(w0);

            yBFNoLMS = yBFNoLMS_samples(1:p.sps:end);
            yBFNoLMS = yBFNoLMS(1:numFrameSym);
            yBFNoLMS = phaseGainNormalize(yBFNoLMS, preSym);

            payloadNoLMS = yBFNoLMS(preLenSym+1:end);
            rxBitsNoLMS = qpskSymbolsToBits(payloadNoLMS, p);

            minLen = min(length(dataBits), length(rxBitsNoLMS));

            nErrBFNoLMS = nErrBFNoLMS + ...
                sum(dataBits(1:minLen) ~= rxBitsNoLMS(1:minLen));

            nBitsBFNoLMS = nBitsBFNoLMS + minLen;

            XtrainSamples = Xsamples(1:trainSampleLen,:);

            [wLMS, eTrain] = spatialLMSRefineSteeringSamples( ...
                XtrainSamples, dTrainSamples, p, w0);

            yBFLMS_samples = Xsamples * conj(wLMS);

            yBFLMS = yBFLMS_samples(1:p.sps:end);
            yBFLMS = yBFLMS(1:numFrameSym);
            yBFLMS = phaseGainNormalize(yBFLMS, preSym);

            payloadLMS = yBFLMS(preLenSym+1:end);
            rxBitsLMS = qpskSymbolsToBits(payloadLMS, p);

            minLen = min(length(dataBits), length(rxBitsLMS));

            nErrBFLMS = nErrBFLMS + ...
                sum(dataBits(1:minLen) ~= rxBitsLMS(1:minLen));

            nBitsBFLMS = nBitsBFLMS + minLen;

            if iSNR == p.exampleSnrIndex && iFrm == 1

                example.EbNo = EbNo;
                example.env = env;
                example.array = array;
                example.p = p;
                example.ch = ch;
                example.preLenSym = preLenSym;
                example.rxBranch = rxBranch;

                example.yBFNoLMS = yBFNoLMS;
                example.yBFLMS = yBFLMS;

                example.yBFNoLMS_samples = yBFNoLMS_samples;
                example.yBFLMS_samples = yBFLMS_samples;

                example.w0 = w0;
                example.wLMS = wLMS;
                example.eTrain = eTrain;

                example.dTrainSamples = dTrainSamples;
                example.dopplerVel = dopplerVel;
            end
        end

        for m = 1:array.M
            berBranches(iSNR,m) = nErrBranch(m) / max(nBitsBranch(m),1);
        end

        berBFNoLMS(iSNR) = nErrBFNoLMS / max(nBitsBFNoLMS,1);
        berBFLMS(iSNR) = nErrBFLMS / max(nBitsBFLMS,1);

        fprintf('Eb/No = %2d dB | BF no LMS = %.4e | BF LMS = %.4e\n', ...
            EbNo, berBFNoLMS(iSNR), berBFLMS(iSNR));
    end


    %% Doppler velocity sweep - one branch with residual Doppler
if doDopplerSweep && spacingCase == 1

    fprintf('\nRunning one-branch residual Doppler sweep for ideal array d = lambda/2\n');

    pSweep = p;
    pSweep.numFrames = numFramesDopplerSweep;
    pSweep.useDoppler = true;
    pSweep.compensateDoppler = false;

    nVel = length(dopplerNominalList);
    nEb  = length(EbNoDopplerList);

    berDopplerArray = NaN(nVel, nEb);   % BF+LMS using 4 hydrophones
    berDopplerBranch = NaN(nVel, nEb);  % affected branch only

    nBitsArray = NaN(nVel, nEb);
    nBitsBranch = NaN(nVel, nEb);

    for iEb = 1:nEb

        EbNoThis = EbNoDopplerList(iEb);

        fprintf('  Eb/N0 = %.1f dB\n', EbNoThis);

        for iVel = 1:nVel

            vResidual = dopplerNominalList(iVel);

            % Common random seed for fair comparison between array and branch
            seedThisPoint = 200000 + iVel;

            rng(seedThisPoint);

            [berArrayTmp, berBranchTmp, bitsArrayTmp, bitsBranchTmp] = ...
                runOneBERPoint_OneBranchResidualDoppler( ...
                    env, array, ch, pSweep, ...
                    preambleBits, preSym, dTrainSamples, ...
                    txFilter, rxFilter, ...
                    w0, EbNoThis, ...
                    vResidual, dopplerBranch);

            berDopplerArray(iVel,iEb) = berArrayTmp;
            berDopplerBranch(iVel,iEb) = berBranchTmp;

            nBitsArray(iVel,iEb) = bitsArrayTmp;
            nBitsBranch(iVel,iEb) = bitsBranchTmp;

            fprintf('    v = %.4f m/s | BF+LMS BER = %.3e | Branch %d BER = %.3e\n', ...
                vResidual, berArrayTmp, dopplerBranch, berBranchTmp);
        end
    end

    allResults.dopplerOneBranch.velocityList = dopplerNominalList;
    allResults.dopplerOneBranch.EbNoList = EbNoDopplerList;
    allResults.dopplerOneBranch.berArray = berDopplerArray;
    allResults.dopplerOneBranch.berBranch = berDopplerBranch;
    allResults.dopplerOneBranch.nBitsArray = nBitsArray;
    allResults.dopplerOneBranch.nBitsBranch = nBitsBranch;
    allResults.dopplerOneBranch.dopplerBranch = dopplerBranch;
    allResults.dopplerOneBranch.arraySpacing_m = array.d;
    allResults.dopplerOneBranch.arraySpacing_lambda = array.d/env.lambda;
end




    %% Result table
    T = table(p.EbNo_dB(:), berTheorySISO(:), berTheoryBF(:), ...
        berBranches(:,1), berBranches(:,2), berBranches(:,3), berBranches(:,4), ...
        berBFNoLMS(:), berBFLMS(:), ...
        'VariableNames', {'EbNo_dB','Theory_SISO','Theory_BF', ...
        'Branch1','Branch2','Branch3','Branch4','BF_No_LMS','BF_LMS'});

    disp(T);

    %% BER and convergence plots
    figure('Name', sprintf('BER and LMS convergence | %s | %s', spacingLabel, p.channelMode));

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
        array.d/env.lambda, p.channelMode));

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
        example.array.d/example.env.lambda));
    legend('|e[n]|', 'Moving average', 'Location', 'northeast');

    %% Signal plots
    figure('Name', sprintf('Signal plots | %s | %s', spacingLabel, p.channelMode));

    payloadNoLMS = example.yBFNoLMS(example.preLenSym+1:end);
    payloadLMS = example.yBFLMS(example.preLenSym+1:end);

    subplot(2,2,1);
    Nplot = min(3000, length(payloadNoLMS));
    scatter(real(payloadNoLMS(1:Nplot)), imag(payloadNoLMS(1:Nplot)), '.');
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('Constellation - BF without LMS | Eb/No = %d dB', example.EbNo));

    subplot(2,2,2);
    Nplot = min(3000, length(payloadLMS));
    scatter(real(payloadLMS(1:Nplot)), imag(payloadLMS(1:Nplot)), '.');
    grid on; axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('Constellation - BF with sample-LMS | Eb/No = %d dB', example.EbNo));

    subplot(2,2,3);
    plotEyeOnAxis(gca, real(example.yBFNoLMS), p.sps, ...
        sprintf('Eye - BF without LMS | d = %.2f \\lambda', ...
        example.array.d/example.env.lambda));

    subplot(2,2,4);
    plotEyeOnAxis(gca, real(example.yBFLMS), p.sps, ...
        sprintf('Eye - BF with sample-LMS | d = %.2f \\lambda', ...
        example.array.d/example.env.lambda));

    %% Beamforming coefficients
    figure('Name', sprintf('Beamforming coefficients | %s | %s', spacingLabel, p.channelMode));

    subplot(2,2,1);
    stem(1:array.M, abs(example.w0), 'filled', 'LineWidth', 1.4); hold on;
    stem(1:array.M, abs(example.wLMS), 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('|w_m|');
    title('Beamforming coefficient magnitude');
    legend('Initial steering weights', 'Sample-LMS refined weights', 'Location', 'best');

    subplot(2,2,2);
    stem(1:array.M, angle(example.w0), 'filled', 'LineWidth', 1.4); hold on;
    stem(1:array.M, angle(example.wLMS), 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('Phase of w_m (rad)');
    title('Beamforming coefficient phase');
    legend('Initial steering weights', 'Sample-LMS refined weights', 'Location', 'best');

    subplot(2,2,3);
    stem(1:array.M, abs(example.wLMS - example.w0), ...
        'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('|w_{LMS} - w_0|');
    title('Magnitude of LMS correction');

    subplot(2,2,4);
    stem(1:array.M, example.rxBranch.bestPhase, ...
        'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('Best sampling phase index');
    title('Best sampling phase per hydrophone');


%% Save results
    results = struct();
    results.env = env;
    results.array = array;
    results.p = p;
    results.ch = ch;
    results.dopplerVel = dopplerVel;
    results.berTheorySISO = berTheorySISO;
    results.berTheoryBF = berTheoryBF;
    results.berBranches = berBranches;
    results.berBFNoLMS = berBFNoLMS;
    results.berBFLMS = berBFLMS;
    results.table = T;
    results.example = example;

    allResults.(spacingLabel) = results;
end

%% Plot - BER versus Doppler velocity with one residual Doppler branch
if doDopplerSweep && isfield(allResults, 'dopplerOneBranch')

    vList = allResults.dopplerOneBranch.velocityList;
    EbList = allResults.dopplerOneBranch.EbNoList;

    berArray = allResults.dopplerOneBranch.berArray;
    berBranch = allResults.dopplerOneBranch.berBranch;

    nBitsArray = allResults.dopplerOneBranch.nBitsArray;
    nBitsBranch = allResults.dopplerOneBranch.nBitsBranch;

    figure('Name','BER versus one-branch residual Doppler');

    semilogy(vList, berArray(:,1), 'o-', 'LineWidth', 1.8); hold on;
    semilogy(vList, berArray(:,2), 'o--', 'LineWidth', 1.8);

    semilogy(vList, berBranch(:,1), 's-', 'LineWidth', 1.8);
    semilogy(vList, berBranch(:,2), 's--', 'LineWidth', 1.8);

    grid on;
    xlabel(sprintf('Residual Doppler velocity on branch %d (m/s)', dopplerBranch));
    ylabel('BER');

    title(sprintf('BER vs residual Doppler on one branch | d = \\lambda/2 | Branch %d affected', ...
        dopplerBranch));

    legend( ...
        sprintf('BF+LMS array | E_b/N_0 = %.0f dB', EbList(1)), ...
        sprintf('BF+LMS array | E_b/N_0 = %.0f dB', EbList(2)), ...
        sprintf('Branch %d only | E_b/N_0 = %.0f dB', dopplerBranch, EbList(1)), ...
        sprintf('Branch %d only | E_b/N_0 = %.0f dB', dopplerBranch, EbList(2)), ...
        'Location', 'best');

    %%  knee detection
    % The knee is defined as the first point where BER rises significantly
    % above the low-Doppler baseline.
    for iEb = 1:length(EbList)

        minObservableArray = 1 / max(nBitsArray(1,iEb),1);
        minObservableBranch = 1 / max(nBitsBranch(1,iEb),1);

        [vKneeArray, berKneeArray] = detectDopplerKnee( ...
            vList, berArray(:,iEb), minObservableArray, kneeFactor, minErrForKnee);

        [vKneeBranch, berKneeBranch] = detectDopplerKnee( ...
            vList, berBranch(:,iEb), minObservableBranch, kneeFactor, minErrForKnee);

        if ~isnan(vKneeArray)
            semilogy(vKneeArray, berKneeArray, 'kp', ...
                'MarkerFaceColor', 'k', 'MarkerSize', 10);
            xline(vKneeArray, 'k:', ...
                sprintf('Array knee %.0f dB: %.4f m/s', EbList(iEb), vKneeArray), ...
                'LabelVerticalAlignment', 'bottom');
        end

        if ~isnan(vKneeBranch)
            semilogy(vKneeBranch, berKneeBranch, 'kd', ...
                'MarkerFaceColor', 'k', 'MarkerSize', 8);
        end

        allResults.dopplerOneBranch.knee.array.v(iEb) = vKneeArray;
        allResults.dopplerOneBranch.knee.array.ber(iEb) = berKneeArray;
        allResults.dopplerOneBranch.knee.branch.v(iEb) = vKneeBranch;
        allResults.dopplerOneBranch.knee.branch.ber(iEb) = berKneeBranch;
    end
end

save('allResults_BFWithWithoutLMS_refactoredDoppler.mat', 'allResults');



function y = matchLength(x, N)

    x = x(:);

    if length(x) >= N
        y = x(1:N);
    else
        y = [x; zeros(N-length(x),1)];
    end
end


function [berArray, berBranch, nBitsArray, nBitsBranch] = ...
    runOneBERPoint_OneBranchResidualDoppler( ...
        env, array, ch, p, ...
        preambleBits, preSym, dTrainSamples, ...
        txFilter, rxFilter, ...
        w0, EbNo_dB, ...
        vResidual, dopplerBranch)

    preLenSym = length(preSym);
    payloadLenSym = p.payloadLengthBits / p.k;

    numFrameSym = preLenSym + payloadLenSym;
    numFrameSamples = numFrameSym * p.sps;

    trainSampleLen = length(dTrainSamples);

    nErrArray = 0;
    nBitsArray = 0;

    nErrBranch = 0;
    nBitsBranch = 0;

    snrSample_dB = EbNo_dB + 10*log10(p.k) - 10*log10(p.sps);

    for iFrame = 1:p.numFrames

        %% Transmitter
        dataBits = randi([0 1], p.payloadLengthBits, 1);

        [~, ~, txWave] = txChain(dataBits, preambleBits, p, txFilter);

        %% SIMO channel
        chanOut = applyChannelSIMO(txWave, ch, array, env);

        %% Residual Doppler on only one branch
        % Other branches are assumed perfectly Doppler compensated.
        if abs(vResidual) > 0
            [dopOut, ~] = f_dopplergeneration( ...
                chanOut(:,dopplerBranch), ...
                vResidual, env.fs, env.c, env.fc);

            chanOut(:,dopplerBranch) = matchLength(dopOut(:), size(chanOut,1));
        end

        %% AWGN
        rxWave = zeros(size(chanOut));

        for m = 1:array.M
            rxWave(:,m) = awgn(chanOut(:,m), snrSample_dB, 'measured');
        end

        %% Receiver front-end
        rxBranch = rxFrontEndSIMO(rxWave, preambleBits, p, rxFilter);

        %% Build aligned sample matrix Xsamples
        Xsamples = zeros(numFrameSamples, array.M);
        frameOK = true;

        for m = 1:array.M

            rxMfNoTransient = rxBranch.rxMfNoTransientCell{m};

            sampleStart = rxBranch.bestPhase(m) + ...
                (rxBranch.frameStartSym(m)-1)*p.sps;

            lastSample = sampleStart + numFrameSamples - 1;

            if sampleStart < 1 || lastSample > length(rxMfNoTransient)
                frameOK = false;
                break;
            end

            Xsamples(:,m) = rxMfNoTransient(sampleStart:lastSample);
        end

        if ~frameOK
            continue;
        end

        %% Single affected branch BER
        frameSymBranch = rxBranch.frameSymbolsCell{dopplerBranch};

        if length(frameSymBranch) >= numFrameSym

            frameSymBranch = frameSymBranch(1:numFrameSym);
            frameSymBranch = phaseGainNormalize(frameSymBranch, preSym);

            payloadBranch = frameSymBranch(preLenSym+1:end);
            rxBitsBranch = qpskSymbolsToBits(payloadBranch, p);

            minLenBranch = min(length(dataBits), length(rxBitsBranch));

            nErrBranch = nErrBranch + ...
                sum(dataBits(1:minLenBranch) ~= rxBitsBranch(1:minLenBranch));

            nBitsBranch = nBitsBranch + minLenBranch;
        end

        %% BF + sample-domain LMS
        XtrainSamples = Xsamples(1:min(trainSampleLen,size(Xsamples,1)),:);
        dTrain = dTrainSamples(1:size(XtrainSamples,1));

        [wLMS, ~] = spatialLMSRefineSteeringSamples( ...
            XtrainSamples, dTrain, p, w0);

        yBFLMS_samples = Xsamples * conj(wLMS);

        yBFLMS = yBFLMS_samples(1:p.sps:end);
        yBFLMS = yBFLMS(1:numFrameSym);

        yBFLMS = phaseGainNormalize(yBFLMS, preSym);

        payloadLMS = yBFLMS(preLenSym+1:end);
        rxBitsLMS = qpskSymbolsToBits(payloadLMS, p);

        minLenArray = min(length(dataBits), length(rxBitsLMS));

        nErrArray = nErrArray + ...
            sum(dataBits(1:minLenArray) ~= rxBitsLMS(1:minLenArray));

        nBitsArray = nBitsArray + minLenArray;
    end

    berArray = nErrArray / max(nBitsArray,1);
    berBranch = nErrBranch / max(nBitsBranch,1);
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