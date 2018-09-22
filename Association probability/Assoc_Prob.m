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

%----- Base Station and user parameters---------------------
Mr_dB = 10 ; 
Mr = 10^(Mr_dB/10);
mr_dB = -10;
mr = 10^(mr_dB/10);
theta_r = 30;
%-------User--------------
Mt_dB = 10 ; 
Mt = 10^(Mt_dB/10);
mt_dB = -10;
mt = 10^(mt_dB/10);
theta_t = 90;

cr = theta_r/360;
ct = theta_t/360;
a = [Mr*Mt, Mr*mt, mr*Mt, mr*mt];
b = [cr*ct, cr*(1-ct), (1-cr)*ct, (1-cr)*(1-ct)];
beta = 1/141.4;

syms x
p(x) = exp(-beta.*x);


j=1;
%% SINR Coverage probability
rc=50:10:300;
hwait = waitbar(0,'Please wait ...');
for i=50:10:300;
lambda = 1/(pi*i^2);
syms r
p(r) = exp(-beta.*r);
BL = 1 - exp(-2*pi*lambda*int(r.*p(r),0,inf));
fL(x) = (2*pi*lambda.*x.*p(x).*exp(-2*pi*lambda*int(r.*p(r),0,x)))./BL;

BN = 1 - exp(-2*pi*lambda*int(r.*(1-p(r)),0,inf));
fN(x) = (2*pi*lambda.*x.*(1-p(x)).*exp(-2*pi*lambda*int(r.*(1-p(r)),0,x)))./BN;

psiL(x) = (CN/CL).^(1/alphaN).*(x.^(alphaL/alphaN));
psiN(x) = (CL/CN).^(1/alphaL).*(x.^(alphaN/alphaL));
%
%------ Association probability-----------------------
AL(j) = vpa(BL.*int(exp(-2*pi*lambda*int(r.*(1-p(r)),0,psiL(x))).*fL(x),0,inf),6);
AN(j)= 1 - AL(j);
j=j+1;
waitbar(i/300,hwait);
end

close(hwait);
plot(rc,AL,'b--o')
grid on 
hold on
plot(rc,AN,'r-o')
xlabel('Avg.Cell radius in meters')
ylabel('Association Probability')
legend('LOS Association Prob.','NLOS Association Prob.')