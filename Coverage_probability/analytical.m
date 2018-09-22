% Coverage and Rate Analysis for Millimeter-Wave
% Cellular Networks
% Tianyang Bai, Student Member, IEEE, and Robert W. Heath, Jr., Fellow, IEEE

fc = 28 ;                                     % GHz
W = 100*1d6 ;                                     % MHz
Q = 10;                                       % Noise figure
n_pdB = -174 +10*log10(W) + Q ;              % noise power
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


%j=1;
%% SINR Coverage probability

rc=100;
lambda = 1/(pi*rc^2);
syms r
p(r) = exp(-beta.*r);
BL = 1 - exp(-2*pi*lambda*int(r.*p(r),0,inf));
fL(x) = (2*pi*lambda.*x.*p(x).*exp(-2*pi*lambda*int(r.*p(r),0,x)))./BL;

BN = 1 - exp(-2*pi*lambda*int(r.*(1-p(r)),0,inf));
fN(x) = (2*pi*lambda.*x.*(1-p(x)).*exp(-2*pi*lambda*int(r.*(1-p(r)),0,x)))./BN;

psiL(x) = ((CN/CL).^(1/alphaN)).*(x.^(alphaL/alphaN));
psiN(x) = (CL/CN).^(1/alphaL).*(x.^(alphaN/alphaL));
%
%------ Association probability-----------------------
AL = vpa(BL.*int(exp(-2*pi*lambda*int(r.*(1-p(r)),0,psiL(x))).*fL(x),0,inf),6);
AN= 1 - AL;
%j=j+1;
%--------------- Probability density function--------------------

FL(x) = (BL.*fL(x)/AL).*(exp(-2*pi*lambda*int(r.*(1-p(r)),0,psiL(x))));
FN(x) = (BN.*fN(x)/AN).*(exp(-2*pi*lambda*int(r.*(1-p(r)),0,psiN(x))));

gnN = N_N*(factorial(N_N))^(-1/N_N);
gnL = N_L*(factorial(N_L))^(-1/N_L);

av = a./(Mt*Mr);
thresh_dB = -5:1:15;
thresh = 10.^(thresh_dB./10);

hwait = waitbar(0,'Please wait...');
for T =1:length(thresh)
sL = 0;
sN = 0;
  for n=1:N_L
      
        qn(x) = 0*x; vn(x) = 0*x;
       
        
        for k=1:4 
        syms x r
        q(x,r)=F(N_L,((n*gnL*av(k)*thresh(T).*(x.^alphaL))./(N_L.*(r.^alphaL)))).*p(r).*r;
        qn(x) = qn(x) +  b(k)*int(q(x,r),r,x,inf);
        v(x,r) = F(N_N,((n*CN*gnL*av(k)*thresh(T).*(x.^alphaL))./(CL*N_N.*(r.^alphaN)))).*(1-p(r)).*r;
        vn(x) = vn(x) +  b(k).*int(v(x,r),r,psiL(x),inf);
       q(x,r)=0*x;
       v(x,r)=0*x;
         end
        qn(x) = 2*pi*lambda*qn(x);
        vn(x) = 2*pi*lambda*vn(x);
         m(x) = exp((-n*gnL.*(x.^alphaL).*thresh(T).*n_p)./(CL*Mr*Mt)-qn(x)-vn(x)).*FL(x);
        sL = sL + ((-1)^(n+1)*factorial(N_L)/(factorial(n)*factorial(N_L-n)))*int(m(x),0,inf);
  end

  
    for n=1:N_N
        wn(x)=0*x; zn(x)=0*x;
        for k=1:4
        w(x,r) = F(N_L,(n*CL*gnN*av(k)*thresh(T)*x.^alphaN/(CN*N_L*r.^alphaL))).*p(r).*r;
        wn(x) = wn(x) + b(k).*int(w(x,r),r,psiN(x),inf);
        z(x,r) = F(N_N,(n*gnN*av(k)*thresh(T)*x.^alphaN/(N_N*r.^alphaN))).*(1-p(r)).*r;
        zn(x) = zn(x) + b(k).*int(z(x,r),r,x,inf);
        w(x,r) = 0;
        z(x,r) = 0;
        end
        wn(x) = 2*pi*lambda*wn(x);
        zn(x) = 2*pi*lambda*zn(x);
        
        c(x) = exp((-n*gnN.*(x.^alphaN).*thresh(T)*n_p)/(CN*Mr*Mt)-wn(x)-zn(x)).*FN(x);
        sN = sN + ((-1)^(n+1)*factorial(N_N)/(factorial(n)*factorial(N_N-n)))*int(c(x),0,inf); 
    end
    PL(T) = vpa(sL,3);
    PN(T) = vpa(sN,3);
    Pc(T) = AL*PL(T) + AN*PN(T);
    waitbar(T/length(thresh),hwait);
 end
 
close(hwait);
plot(thresh_dB,Pc)




