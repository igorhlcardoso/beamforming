function a = steeringVectorImperfect(thetaDeg, array, env, err)
% steeringVectorImperfect
% Generates an imperfect steering vector for testing LMS convergence.
%
% Inputs:
%   thetaDeg : nominal steering angle in degrees
%   array    : array structure
%   env      : environment structure
%   err      : structure with optional imperfections:
%              err.thetaOffsetDeg
%              err.positionErrorStd_m
%              err.phaseErrorStd_rad
%              err.gainErrorStd_dB
%              err.seed
%
% Output:
%   a        : imperfect steering vector, normalized by a(1)

    if nargin < 4
        err = struct();
    end

    if isfield(err, 'seed') && ~isempty(err.seed)
        rng(err.seed);
    end

    % Default values
    thetaOffsetDeg = getFieldDefault(err, 'thetaOffsetDeg', 0);
    positionErrorStd_m = getFieldDefault(err, 'positionErrorStd_m', 0);
    phaseErrorStd_rad = getFieldDefault(err, 'phaseErrorStd_rad', 0);
    gainErrorStd_dB = getFieldDefault(err, 'gainErrorStd_dB', 0);

    % Imperfect steering angle
    thetaUsedDeg = thetaDeg + thetaOffsetDeg;
    thetaRad = deg2rad(thetaUsedDeg);

    % Imperfect hydrophone positions
    pos = array.positions(:);

    if positionErrorStd_m > 0
        pos = pos + positionErrorStd_m * randn(size(pos));
        pos(1) = array.positions(1); % keep H1 as reference
    end

    % Geometric delay
    tau = pos * sin(thetaRad) / env.c;

    % Ideal phase from imperfect geometry/angle
    a = exp(-1j * 2*pi * env.fc * tau);

    % Additional random phase mismatch per hydrophone
    if phaseErrorStd_rad > 0
        phaseErr = phaseErrorStd_rad * randn(size(a));
        phaseErr(1) = 0; % keep H1 reference unchanged
        a = a .* exp(1j * phaseErr);
    end

    % Additional random gain mismatch per hydrophone
    if gainErrorStd_dB > 0
        gainErr_dB = gainErrorStd_dB * randn(size(a));
        gainErr_dB(1) = 0; % keep H1 reference unchanged
        gainLin = 10.^(gainErr_dB/20);
        a = gainLin .* a;
    end

    % Normalize by first hydrophone
    a = a / a(1);
end

function value = getFieldDefault(s, fieldName, defaultValue)
    if isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end