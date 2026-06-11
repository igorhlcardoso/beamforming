clc; clear; close all;
rng(2025);

allResults = struct();

for spacingCase = 1:2

    % Physical parameters
    env.c = 1500;
    env.fc = 26e3;
    env.fsADC = 156250;
    env.decim = 8;
    env.fs = env.fsADC / env.decim;
    env.lambda = env.c / env.fc;

    % Array parameters
    array.M = 4;
    array.thetaDeg = 20;
    array.thetaRad = deg2rad(array.thetaDeg);

    switch spacingCase
        case 1
            array.d = env.lambda/2;
            spacingLabel = 'd_lambda_over_2';

        case 2
            array.d = 0.25;
            spacingLabel = 'd_25cm';
    end

    array.positions = (0:array.M-1) * array.d;

    % Modulation parameters
    p.M = 4;
    p.k = log2(p.M);
    p.sps = 10;
    p.Rs = env.fs / p.sps;
    p.phaseOffset = pi/4;

    % Rectangular pulse
    txFilter = ones(p.sps,1) / sqrt(p.sps);
    rxFilter = txFilter;

    % Frame parameters
    p.trainingLengthSymbols = 1024;
    p.preambleLengthBits = p.trainingLengthSymbols * p.k;
    p.payloadLengthBits = 5000;

    % Monte Carlo parameters
    p.EbNo_dB = 0:7;
    p.numFrames = 20;
    p.exampleSnrIndex = 4;

    % Receiver and LMS parameters
    p.useBestSamplingPhase = true;
    p.mu = 0.05;
    p.normalizeLMS = true;

    % Channel mode
    p.channelMode = 'multipath';
    %p.channelMode = 'direct';

    % Doppler parameters
    p.useDoppler = false;
    p.compensateDoppler = false;

    dopplerVel = [0.005 0.004 0.003 0.006];

    if length(dopplerVel) ~= array.M
        error('dopplerVel must have one value per hydrophone.');
    end

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

            ch.gains = [1.0 0.6 0.25 0.1];
            ch.gains = ch.gains / sqrt(sum(abs(ch.gains).^2));

            ch.phases = [0.00 0.60 -0.80 1.10];

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

    delaySym = ch.delaysSec * p.Rs;

    disp(table((1:ch.Npaths).', ch.delaysSec(:), delaySym(:), ...
        ch.gains(:), ch.phases(:), ch.thetaDeg(:), ...
        'VariableNames', {'Path','Delay_s','Delay_symbols','Gain','Phase_rad','AoA_deg'}));

    % Preamble
    preambleBits = randi([0 1], p.preambleLengthBits, 1);
    preSym = qpskBitsToSymbols(preambleBits, p);
    preLenSym = length(preSym);

    % Sample-rate preamble reference
    dTrainSamples = buildPreambleSampleReference(preSym, txFilter, rxFilter, p);
    trainSampleLen = length(dTrainSamples);

    % Doppler compensation windows
    p.dopplerWinLen = 512;
    p.dopplerIdx1 = round(0.25 * length(dTrainSamples));
    p.dopplerIdx2 = round(0.75 * length(dTrainSamples));

    % Theoretical BER curves
    berTheorySISO = berawgn(p.EbNo_dB, 'psk', p.M, 'nondiff');
    berTheoryBF = berawgn(p.EbNo_dB + 10*log10(array.M), ...
        'psk', p.M, 'nondiff');

    % Storage
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
    fprintf('Doppler enabled = %d | Doppler compensation enabled = %d\n', ...
        p.useDoppler, p.compensateDoppler);

    % Deterministic steering vector
    a0 = steeringVector(array.thetaDeg, array, env);
    w0 = a0 / array.M;

    for iSNR = 1:length(p.EbNo_dB)

        EbNo = p.EbNo_dB(iSNR);

        nErrBranch = zeros(1,array.M);
        nBitsBranch = zeros(1,array.M);

        nErrBFNoLMS = 0;
        nBitsBFNoLMS = 0;

        nErrBFLMS = 0;
        nBitsBFLMS = 0;

        for iFrm = 1:p.numFrames

            % Transmitter
            dataBits = randi([0 1], p.payloadLengthBits, 1);
            [~, ~, txWave] = txChain(dataBits, preambleBits, p, txFilter);

            % SIMO channel
            chanOut = applyChannelSIMO(txWave, ch, array, env);

            % Doppler per branch
            if p.useDoppler
                for m = 1:array.M
                    [dopOut, ~] = f_dopplergeneration( ...
                        chanOut(:,m), ...
                        dopplerVel(m), ...
                        env.fs, ...
                        env.c, ...
                        env.fc);

                    if p.compensateDoppler
                        compOut = f_dopplercompensation(dopOut(:), ...
                            dopplerVel(m), env.fs, env.c, env.fc);
                        chanOut(:,m) = compOut(:);
                    else 
                        chanOut(:,m) = dopOut(:);
                    end

                end
            end

            % AWGN
            snrSample_dB = EbNo + 10*log10(p.k) - 10*log10(p.sps);

            rxWave = zeros(size(chanOut));

            for m = 1:array.M
                rxWave(:,m) = awgn(chanOut(:,m), snrSample_dB, 'measured');
            end

            % Receiver front-end
            rxBranch = rxFrontEndSIMO(rxWave, preambleBits, p, rxFilter);

            % Build sample-time matrix aligned by branch synchronization
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
             
            % Individual branch BERs after sample-time Doppler compensation
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

            % Beamforming without LMS at sample time
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

            % Beamforming with spatial LMS trained at sample time
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

            % Save one frame for plots
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
                %example.dopplerEstHz = dopplerEstHz;
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

    T = table(p.EbNo_dB(:), berTheorySISO(:), berTheoryBF(:), ...
        berBranches(:,1), berBranches(:,2), berBranches(:,3), berBranches(:,4), ...
        berBFNoLMS(:), berBFLMS(:), ...
        'VariableNames', {'EbNo_dB','Theory_SISO','Theory_BF', ...
        'Branch1','Branch2','Branch3','Branch4','BF_No_LMS','BF_LMS'});

    disp(T);

    % BER and LMS convergence
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

    % Doppler correction plot
    figure('Name', sprintf('Doppler correction | %s | %s', spacingLabel, p.channelMode));

    % for m = 1:array.M
    %     subplot(2,2,m);
    %     plotDopplerCorrectionOnAxis(gca, ...
    %         example.XsamplesBeforeDopplerComp(:,m), ...
    %         example.XsamplesAfterDopplerComp(:,m), ...
    %         example.dTrainSamples, ...
    %         example.env.fs, ...
    %         sprintf('Branch %d | v = %.4f m/s | fEst = %.4f Hz', ...
    %         m, example.dopplerVel(m), example.dopplerEstHz(m)));
    % end

    % Signal plots
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

    % Coefficients and synchronization
    figure('Name', sprintf('Beamforming coefficients | %s | %s', spacingLabel, p.channelMode));

    subplot(2,2,1);
    stem(1:example.array.M, abs(example.w0), 'filled', 'LineWidth', 1.4); hold on;
    stem(1:example.array.M, abs(example.wLMS), 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('|w_m|');
    title('Beamforming coefficient magnitude');
    legend('Initial steering weights', 'Sample-LMS refined weights', 'Location', 'best');

    subplot(2,2,2);
    stem(1:example.array.M, angle(example.w0), 'filled', 'LineWidth', 1.4); hold on;
    stem(1:example.array.M, angle(example.wLMS), 'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('Phase of w_m (rad)');
    title('Beamforming coefficient phase');
    legend('Initial steering weights', 'Sample-LMS refined weights', 'Location', 'best');

    subplot(2,2,3);
    stem(1:example.array.M, abs(example.wLMS - example.w0), ...
        'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('|w_{LMS} - w_0|');
    title('Magnitude of LMS correction');

    subplot(2,2,4);
    stem(1:example.array.M, example.rxBranch.bestPhase, ...
        'filled', 'LineWidth', 1.4);
    grid on;
    xlabel('Hydrophone index');
    ylabel('Best sampling phase index');
    title('Best sampling phase per hydrophone');

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

save('allResults_BFWithWithoutLMS.mat', 'allResults');

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

% function plotDopplerCorrectionOnAxis(ax, xBefore, xAfter, dRef, Fs, ttl)
% 
%     xBefore = xBefore(:);
%     xAfter = xAfter(:);
%     dRef = dRef(:);
% 
%     L = min([length(xBefore), length(xAfter), length(dRef)]);
% 
%     xBefore = xBefore(1:L);
%     xAfter = xAfter(1:L);
%     dRef = dRef(1:L);
% 
%     t = (0:L-1).' / Fs;
% 
%     phaseBefore = unwrap(angle(xBefore .* conj(dRef)));
%     phaseAfter = unwrap(angle(xAfter .* conj(dRef)));
% 
%     phaseBefore = phaseBefore - mean(phaseBefore(1:min(100,length(phaseBefore))));
%     phaseAfter = phaseAfter - mean(phaseAfter(1:min(100,length(phaseAfter))));
% 
%     phaseBeforeSmooth = movmean(phaseBefore, 128);
%     phaseAfterSmooth = movmean(phaseAfter, 128);
% 
%     cla(ax);
%     plot(ax, t, phaseBeforeSmooth, 'LineWidth', 1.2); hold(ax, 'on');
%     plot(ax, t, phaseAfterSmooth, 'LineWidth', 1.2);
%     grid(ax, 'on');
% 
%     xlabel(ax, 'Time in training sequence (s)');
%     ylabel(ax, 'Residual phase (rad)');
%     title(ax, ttl);
%     legend(ax, 'Before compensation', 'After compensation', 'Location', 'best');
% end