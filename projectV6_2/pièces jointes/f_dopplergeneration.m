function [v_out, v_outcorr] = f_dopplergeneration(v_in,d_dop,d_Fech,d_cel,d_Fporteuse)
% application d'un effet doppler constante par re-echantillonnage uniforme sur un signal
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
% ex : [v_out, v_outcorr] = f_dopplergeneration(v_in,2.5,192E3,1.5E3,26E3)
% 
% version 1.0 - 2022/06/23
%
% Contributeur : N. GROLLIER
%
% tous droits reserves (iXblue)

fo_dop      = (1-d_dop/d_cel);
v_tdop      = (0:numel(v_in)-1)*1/(d_Fech*fo_dop);
v_tf0       = (0:numel(v_in)-1)*1/(d_Fech);
% v_tf0       = v_tf0(v_tf0<=v_tdop(end)); % eviter extrapolation si DOP negatif
v_ininterp  = interp1(v_tf0,v_in,v_tdop,'spline',0);

v_out       = v_ininterp.*exp(-2*1i*pi*d_Fporteuse/d_Fech*(0:numel(v_ininterp)-1)*(-d_dop/d_cel));
v_outcorr   = conj(v_ininterp(end:-1:1).* exp(-2*1i*pi*d_Fporteuse/d_Fech*(0:numel(v_ininterp)-1)*(d_dop/d_cel)));