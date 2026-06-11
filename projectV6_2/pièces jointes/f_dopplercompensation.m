function v_out = f_dopplercompensation(v_in,d_dop,d_Fech,d_cel,d_Fporteuse)
% compensation d'une vitesse doppler constante par re-echantillonnage uniforme sur un signal
% 
% input :   v_in  - signal d'interet (vecteur) 
%           d_dop - vitesse doppler (en m/s)
%           d_Fech - frequence d'echantillonnage de v_in
%           d_cel - celerite du son dans l'eau
%           d_Fporteuse - frequence centrale de v_in [Hz]
% 
% output :  v_out - v_in affecte d'un doppler d_dop
%           v_outcorr - v_in affecte d'un doppler d_dop utilisable pour
%           correlation sous forme de filtre adapte
% 
% ex : [v_out, v_outcorr] = f_dopplercompensation(v_in,2.5,19531.25,1.5E3,26E3)
% 
% version 1.0 - 2022/06/23
%
% Contributeur : N. GROLLIER
%
% tous droits reserves (iXblue)

fo_dop          = (1-d_dop/d_cel);
v_tdop          = (1:numel(v_in))*1/(d_Fech*fo_dop);
v_tf0           = (1:numel(v_in))*1/(d_Fech);
v_tf0           = v_tf0(v_tf0<=v_tdop(end)); % eviter extrapolation si DOP negatif
v_in_interp     = interp1(v_tdop,v_in,v_tf0,'spline');
v_out           = v_in_interp.*exp(-2*1i*pi*d_Fporteuse/d_Fech*(1:numel(v_in_interp))*(d_dop/d_cel));

% figure; hold all;
% plot(v_tf0,real(v_out))
% plot(v_tdop,real(v_in))

end