function v_toa = computeGeometricTOA(st_param)

v_toa = [];
d_ReflecteurZ = 0;

for NN = 1:3

    Z(1) = 2*(NN-1)*(st_param.HauteurEau - d_ReflecteurZ) + ...
        st_param.Zs - st_param.Zr;

    Z(2) = 2*(NN-1)*(st_param.HauteurEau - d_ReflecteurZ) + ...
        st_param.Zs + st_param.Zr;

    Z(3) = 2*NN*(st_param.HauteurEau - d_ReflecteurZ) - ...
        st_param.Zs - st_param.Zr;

    Z(4) = 2*NN*(st_param.HauteurEau - d_ReflecteurZ) - ...
        st_param.Zs + st_param.Zr;

    Ranges = sqrt(Z.^2 + st_param.DistanceHoriz^2);

    v_toa = [v_toa, Ranges / st_param.Celerite];
end

v_toa = sort(v_toa(:)).';

end