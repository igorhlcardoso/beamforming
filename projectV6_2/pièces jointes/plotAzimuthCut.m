function a0 = plotAzimuthCut(array, env, theta0Deg)

thetaScan = -90:0.1:90;

a0 = steeringVector(theta0Deg, array, env);

response = zeros(size(thetaScan));

for i = 1:length(thetaScan)

    aTheta = steeringVector(thetaScan(i), array, env);

    % Matched spatial response:
    % response(theta) = |a(theta0)^H a(theta)| / M
    response(i) = abs(a0' * aTheta) / array.M;
end

response = response / max(response);
response_dB = 20*log10(response + 1e-12);

figure;
plot(thetaScan, response_dB, 'LineWidth', 1.5);
grid on;
xlabel('Azimuth angle (deg)');
ylabel('Normalized array response (dB)');
title(sprintf('Azimuth cut - steering to %.1f° | d = %.2f \\lambda', ...
    theta0Deg, array.d/env.lambda));
ylim([-40 0]);

end