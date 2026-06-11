function st = validateSteeringVector(rxBranch, preSym, array, env)

M = array.M;
preLenSym = length(preSym);

%% Geometric delay vector
tauGeom = array.positions(:) * sin(array.thetaRad) / env.c;

%% Theoretical steering vector from geometric delay
aGeom = exp(-1j * 2*pi * env.fc * tauGeom);
aGeom = aGeom / aGeom(1);

%% Estimated steering vector from received preamble
aEst = zeros(M,1);

for m = 1:M

    frameSym = rxBranch.frameSymbolsCell{m};

    if length(frameSym) < preLenSym
        error('Frame too short to estimate steering vector.');
    end

    rxPre = frameSym(1:preLenSym);

    alpha = sum(rxPre .* conj(preSym));

    aEst(m) = alpha;
end

aEst = aEst / aEst(1);

%% Phase comparison
phaseGeom = unwrap(angle(aGeom));
phaseEst = unwrap(angle(aEst));

phaseError = angle(aEst .* conj(aGeom));

%% Table
st = table((1:M).', array.positions(:), tauGeom(:), tauGeom(:)*1e6, ...
    abs(aGeom(:)), angle(aGeom(:)), abs(aEst(:)), angle(aEst(:)), ...
    phaseError(:), ...
    'VariableNames', {'Sensor','Position_m','TauGeom_s','TauGeom_us', ...
    'Abs_Geom','Phase_Geom_rad','Abs_Est','Phase_Est_rad', ...
    'Phase_Error_rad'});

disp(st);

%% Plot phase
figure;
plot(1:M, phaseGeom, 'o-', 'LineWidth', 1.5); hold on;
plot(1:M, phaseEst, 's--', 'LineWidth', 1.5);
grid on;
xlabel('Hydrophone index');
ylabel('Phase (rad)');
title('Steering Vector Phase: Geometric vs Estimated');
legend('Geometric/theoretical', 'Estimated from preamble', ...
    'Location', 'best');

%% Plot magnitude
figure;
plot(1:M, abs(aGeom), 'o-', 'LineWidth', 1.5); hold on;
plot(1:M, abs(aEst), 's--', 'LineWidth', 1.5);
grid on;
xlabel('Hydrophone index');
ylabel('Magnitude');
title('Steering Vector Magnitude: Geometric vs Estimated');
legend('Geometric/theoretical', 'Estimated from preamble', ...
    'Location', 'best');

%% Plot phase error
figure;
stem(1:M, phaseError, 'filled', 'LineWidth', 1.5);
grid on;
xlabel('Hydrophone index');
ylabel('Phase error (rad)');
title('Estimated Steering Vector Phase Error');

end