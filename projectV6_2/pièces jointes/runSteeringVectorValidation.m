clc; clear; close all;
rng(2026);

resultsDir = 'results_BFWithWithoutLMS';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Physical parameters
env.c = 1500;
env.fc = 26e3;

env.fsADC = 156250;
env.decim = 8;
env.fs = env.fsADC / env.decim;

env.lambda = env.c / env.fc;

%% Array
array.M = 4;
array.thetaDeg = 20;
array.thetaRad = deg2rad(array.thetaDeg);

spacingCases = [1 2];

%% Modulation
p.M = 4;
p.k = log2(p.M);
p.sps = 8;
p.Rs = env.fs / p.sps;
p.phaseOffset = pi/4;

p.trainingLengthSymbols = 1024;
p.preambleLengthBits = p.trainingLengthSymbols * p.k;
p.payloadLengthBits = 4000;

p.useBestSamplingPhase = true;

pulse = ones(p.sps,1) / sqrt(p.sps);

%% High SNR
EbNo = 50;

for iCase = 1:length(spacingCases)

    spacingCase = spacingCases(iCase);

    switch spacingCase
        case 1
            array.d = env.lambda/2;
            spacingLabel = 'lambda_over_2';

        case 2
            array.d = 0.25;
            spacingLabel = '25cm';
    end

    array.positions = (0:array.M-1) * array.d;

    %% Single direct path
    ch.Npaths = 1;
    ch.delaysSec = 0;
    ch.gains = 1;
    ch.phases = 0;
    ch.thetaDeg = array.thetaDeg;
    ch.thetaRad = array.thetaRad;

    %% Preamble
    preambleBits = randi([0 1], p.preambleLengthBits, 1);
    preSym = qpskBitsToSymbols(preambleBits, p);

    %% TX
    dataBits = randi([0 1], p.payloadLengthBits, 1);
    [~, ~, txWave] = txChain(dataBits, preambleBits, p, pulse);

    %% Channel
    chanOut = applyChannelSIMO(txWave, ch, array, env);

    %% AWGN
    snrSample_dB = EbNo + 10*log10(p.k) - 10*log10(p.sps);

    rxWave = zeros(size(chanOut));

    for m = 1:array.M
        rxWave(:,m) = awgn(chanOut(:,m), snrSample_dB, 'measured');
    end

    %% RX
    rxBranch = rxFrontEndSIMO(rxWave, preambleBits, p, pulse);

    %% Validate steering vector
    st = validateSteeringVector(rxBranch, preSym, array, env);

    %% Plot azimuth cut
    plotAzimuthCut(array, env, array.thetaDeg);

    %% Save
    save(fullfile(resultsDir, ...
        sprintf('steering_validation_%s.mat', spacingLabel)), ...
        'st', 'array', 'env', 'p', 'ch', 'rxBranch', 'preSym');
end
