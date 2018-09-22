% Coverage and Rate Analysis for Millimeter-Wave
% Cellular Networks
% Tianyang Bai, Student Member, IEEE, and Robert W. Heath, Jr., Fellow, IEEE

fc = 28 ;                                     % GHz
W = 100 ;                                     % MHz
Q = 10;                                       % Noise figure
n_pdB = -174 +10*10*log10(W) + Q ;              % noise power
n_p = 10^(n_pdB/10);
alphaL = 2;                                   % LOS path loss exponent
alphaN =4 ;
betaL = 61.4;
betaN = 72;
CL = 10^(-betaL/10);
CN = 10^(-betaN/10);
N_L = 3;
N_N = 2 ;
beta = 1/141.4 ;
rc = 100;
%----- Base Station
Mr_dB = 10 ; 
Mr = 10^(Mr_dB/10);
mr_dB = -10;
mr = 10^(mr_dB/10);
theta_r = 90;
%-------User--------------
Mt_dB = 10 ;
%Mt_dB = 20 ;
Mt = 10^(Mt_dB/10);
mt_dB = -10;
mt = 10^(mt_dB/10);
% theta_t = 30;
theta_t = 45;

cr = theta_r/360;
ct = theta_t/360;
a = [Mr*Mt, Mr*mt, mr*Mt, mr*mt];
b = [cr*ct, cr*(1-ct), (1-cr)*ct, (1-cr)*(1-ct)];
R = -141.4*log(b);
hsq0 = gamrnd(N_L,1/N_L);
%L_R0 = CL*abs(R(1)).^(-alphaL);              
L_R0dB = 10*alphaL*log10(R(1));
L_R0 = 10.^(L_R0dB./10);
som = 0;
%L_LOS = CL*R.^(-alphaL) ;                   
L_LOSdB = 10*alphaL*log10(R);
L_LOS = 10.^(L_LOSdB./10);
L_NLOS = CN*abs(R).^(-alphaN) ;
thresh_dB = 5:1:40;
thresh = 10.^(thresh_dB./10);

hwait = waitbar(0,'Please wait');

for i = 1:length(thresh)
    d=0;
    for j=1:10000
        som =0;
        
        for l=2:4
            hsq=gamrnd(N_L,1/N_L);
            som = som + hsq*a(l)*L_LOS(l);
        end
%         test(j)=som;
%         %som_dB = 10*log10(som);
% 
         SINR(i) = (hsq0*Mr*Mt*L_R0)./(n_p+som);
%         
         if SINR(i) > thresh(i)
             d = d+1;
         end
    end
    Pc(i)=d/j;
    waitbar(i/length(thresh),hwait);
end
close(hwait);
%plot(thresh_dB,Pc,'b--o')
%plot(thresh_dB,Pc,'r-d')
plot(thresh_dB,Pc,'g->')
xlabel('SINR Threshold in dB')
ylabel('SINR Coverage Probability')
grid on
hold on
%legend('(Mt,mt,?t)=(10dB,-10dB,30º)','(Mt,mt,?t)=(20dB,-10dB,30º)','(Mt,mt,?t)=(10dB,-10dB,45º)')