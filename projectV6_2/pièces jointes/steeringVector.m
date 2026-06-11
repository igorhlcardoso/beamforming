function a = steeringVector(thetaDeg, array, env)

thetaRad = deg2rad(thetaDeg);
%geometric delay
tau = array.positions(:) * sin(thetaRad) / env.c;
%converts it to complex phase
a = exp(-1j * 2*pi * env.fc * tau);

a = a / a(1);

end