function Y = applyChannelSIMO(x, ch, array, env)

x = x(:);
N = length(x);
M = array.M;

n = (0:N-1).';
t = n / env.fs;

maxSpatialDelay = 0;

for pth = 1:ch.Npaths
    maxSpatialDelay = max(maxSpatialDelay, ...
        max(abs(array.positions * sin(ch.thetaRad(pth)) / env.c)));
end

maxDelaySec = max(ch.delaysSec) + maxSpatialDelay;
maxDelaySamples = ceil(maxDelaySec * env.fs) + 4;

Y = zeros(N + maxDelaySamples, M);

for m = 1:M

    sensorPos = array.positions(m);

    y_m = zeros(N + maxDelaySamples, 1);

    for pth = 1:ch.Npaths

        thetaPath = ch.thetaRad(pth);
%delay
        tauSpatial = sensorPos * sin(thetaPath) / env.c;

            %test force samples to entire numbers
        tauTotal = ch.delaysSec(pth) + tauSpatial;
        %tauTotal = round(tauTotal * env.fs) / env.fs;
%phase
        phiSpatial = -2*pi*env.fc*tauSpatial;
        phiTotal = ch.phases(pth) + phiSpatial;

%temporal delay using interpolation
        tQuery = t - tauTotal;


        delayed = interp1(t, x, tQuery, 'linear', 0);
        %test
        %delayed = interp1(t, x, tQuery, 'spline', 0);



        pathSig = ch.gains(pth) * exp(1j*phiTotal) * delayed;

        y_m(1:length(pathSig)) = y_m(1:length(pathSig)) + pathSig;
    end

    Y(:,m) = y_m;
end

end