function [TOA_samples, AMP, delta, ABC] = P_InterpParaboTiming(Rxx, IDEch)

% Parabolic interpolation around a correlation peak.
%
% Inputs:
%   Rxx   : correlation vector, real or complex
%   IDEch : integer peak index
%
% Outputs:
%   TOA_samples : interpolated peak position in sample-index units
%   AMP         : interpolated peak amplitude
%   delta       : fractional offset relative to IDEch
%   ABC         : parabola coefficients

    TOA_samples = NaN;
    AMP = NaN;
    delta = NaN;
    ABC = [NaN; NaN; NaN];

    Rxx = Rxx(:);
    N = length(Rxx);

    if isempty(Rxx) || ~isfinite(IDEch)
        return;
    end

    IDEch = round(IDEch);

    if IDEch <= 1 || IDEch >= N
        return;
    end

    % Use local coordinates x = -1, 0, +1
    MAT = [ ...
        1  -1  1; ...
        0   0  1; ...
        1   1  1];

    Yind = abs([Rxx(IDEch-1); Rxx(IDEch); Rxx(IDEch+1)]);

    ABC = MAT \ Yind;

    A = ABC(1);
    B = ABC(2);
    C = ABC(3);

    % If the parabola is flat or opens upward, the interpolation is unreliable.
    % For a true local maximum, A should be negative.
    if ~isfinite(A) || ~isfinite(B) || abs(A) < eps || A >= 0
        delta = 0;
        TOA_samples = IDEch;
        AMP = abs(Rxx(IDEch));
        return;
    end

    delta = -B/(2*A);

    % Safety: the vertex must remain inside the local interval.
    if ~isfinite(delta) || abs(delta) > 1
        delta = 0;
        TOA_samples = IDEch;
        AMP = abs(Rxx(IDEch));
        return;
    end

    TOA_samples = IDEch + delta;
    AMP = A*delta^2 + B*delta + C;
end
