function [TOA,POA,AMP,ABC]=P_InterpParabo(Rxx,IDEch,Fs)

% Inititialisation des variables de sortie
TOA                 = 0.0;
POA                 = 0.0;
AMP                 = 0.0;

% Controle de coherence - Domaine [1-2;1023-1024]
if ((IDEch>1) && (IDEch<length(Rxx)))
    
    % Estimation du TOA
   %%% IDEch=399;
    % MAT                 = [ (IDEch-1)^2         (IDEch-1)^1         (IDEch-1)^0;...
    %     (IDEch-0)^2         (IDEch-0)^1         (IDEch-0)^0;...
    %     (IDEch+1)^2         (IDEch+1)^1         (IDEch+1)^0];
     MAT                 = [ (-1)^2         (-1)^1         (-1)^0;...
        (0)^2         (0)^1         (0)^0;...
        (1)^2         (1)^1         (1)^0];
    Yind                = [ Rxx(IDEch-1)        Rxx(IDEch)          Rxx(IDEch+1)]';
    ABC                 = MAT\abs(Yind);
    
    TOA                 = (-ABC(2)/2/ABC(1))./Fs+IDEch/Fs;
    TOA_local           = (-ABC(2)/2/ABC(1))./Fs;
    AMP                 = ABC(1)*( TOA_local.*Fs).^2 + ABC(2)*(TOA_local.*Fs) + ABC(3);

    %% Affichage Debug
    % M=abs(Yind(1));
    % N=abs(Yind(2));
    % O=abs(Yind(3));
    %
    %
    % A=(O-2*N+M)/2;
    % B=(O-M)/2;
    % C=N;
    % Sommet=-(B/(2*A));
    %
    %  TOA_interp=IDEch/Fs+Sommet/Fs;
    %  AMP_interp=A*Sommet^2+B*Sommet+C; % -(B^2-4*A*C)/(4*A)
    %
    % figure;
    % Sample=IDEch-20:IDEch+20;
    % plot(Sample,abs(Rxx(Sample)),'b+:')
    % hold on;
    % plot([IDEch-1:IDEch+1],abs([M N O]),'ro')
    %
    % XX=-1:0.1:1;
    %
    % P_XX=A.*XX.^2+B.*XX+C;
    % hold on;
    % plot(XX+IDEch,P_XX,'g')
    %
    % hold on;
    %
    % plot(TOA_interp*Fs,AMP_interp,'ks','MarkerFaceColor','g')

     
    % Estimation du POA
    Pxx                 = Rxx([floor(TOA.*Fs) ceil(TOA.*Fs)]);
    Pxx                 = atan2(imag(Pxx),real(Pxx));
    
%     figure
%      plot(Sample,atan2(imag(Rxx(Sample)),real(Rxx(Sample))),'b+:')
%      hold on;
%      plot([floor(TOA.*Fs) ceil(TOA.*Fs)],Pxx,'ro')
     
    if abs(Pxx(1)-Pxx(2))>pi
        % Ajustement de Pxx(2) au domaine de Pxx(1)
        if Pxx(1)>Pxx(2)
            Pxx(2)          = Pxx(2)+2*pi;
        else
            Pxx(2)          = Pxx(2)-2*pi;            
        end
    end

    % Interpolation lineaire
     POA                 = mod(interp1([0 1],Pxx,TOA.*Fs-floor(TOA.*Fs),'linear'),2*pi);
     
else
    fprintf('Erreur de domaine\n');
end