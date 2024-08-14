close all
clear
clc
format short

% parameter1
X = 2;
T = 3;
dx = 0.01;
dt = 0.00005;
M=round(X/dx); 
N=round(T/dt);
Nt1 = linspace(0, T, N+1);
Nt = linspace(0, T, N+1);
Nx = linspace(0, X, M+1);
[Nt, Nx] = meshgrid(Nt, Nx);

% parameter2
Du = 0.001;
Da = 0.001;
LA1 = 0.08;
LA2 = 0.02;
alpha2 = 0.6;
alpha1 = 0.4;
beta1 = 0.2;
beta2 = 0.5;
d1 = 0.0005;
d2 = 0.0005;
gamma1 = 0.1;
gamma3 = 0.3;
gamma4 = 0.3;


A1 =1;
A2 =1;
A3 =1;
A4 = 1;
A5 = 1;

B1 = 1;
B2 = 1;
B3 = 1;
B4 = 1;
B5 = 1;

C1 = 0.5;
C2 = 0.3;
C3 = 0.5;
C4 = 1.8;
C5 = 5.5;
C6 = 2.1;

Iud = 1.2;
Iad = 1.2;
LSud = 0.5;
LIud = 0.2;
LRud = 0.8;

J_sum = 0;
precision1 = 0.001;
precision2 = 10000;


K1 = zeros(M+1,1);
K2 = zeros(M+1,1);
K12 = zeros(M+1,1);
K21 = zeros(M+1,1);
Su = zeros(M+1,N+1);
Iu = zeros(M+1,N+1);
Ru = zeros(M+1,N+1);
LSu = zeros(M+1,N+1);
LIu = zeros(M+1,N+1);
LRu = zeros(M+1,N+1);
Sa = zeros(M+1,N+1);
Ia = zeros(M+1,N+1);
Ra = zeros(M+1,N+1);
Sus = zeros(M+1,N+1);

oldSu = zeros(M+1,N+1);
oldIu = zeros(M+1,N+1);
oldRu = zeros(M+1,N+1);
oldLSu = zeros(M+1,N+1);
oldLIu = zeros(M+1,N+1);
oldLRu = zeros(M+1,N+1);
oldSa = zeros(M+1,N+1);
oldIa = zeros(M+1,N+1);
oldRa = zeros(M+1,N+1);
oldSus = zeros(M+1,N+1);

%%%%ajoint variables
lamda1 = zeros(M+1, N+1);
lamda2 = zeros(M+1, N+1);
lamda3 = zeros(M+1, N+1);
lamda4 = zeros(M+1, N+1);
lamda5 = zeros(M+1, N+1);
lamda6 = zeros(M+1, N+1);
lamda7 = zeros(M+1, N+1);
lamda8 = zeros(M+1, N+1);
lamda9 = zeros(M+1, N+1);
lamda10 = zeros(M+1, N+1);
num = 0;

%%control strategies
alpha3 = zeros(M+1, N+1);%%u1
beta3 = zeros(M+1, N+1);%%u2
gamma2 = zeros(M+1, N+1);%%u3
u4 = zeros(M+1, N+1);
u5 = zeros(M+1, N+1);
u6 = zeros(M+1, N+1);

for i = 1:M+1
 x = i*dx;
 Su(i, 1) = 9 + 5*sin(2*x);
 Iu(i, 1) = 3 + 3*sin(20*x);
 Ru(i, 1) = 0;
 LSu(i, 1) = 0;
 LIu(i, 1) = 0;
 LRu(i, 1) = 0;
 Sa(i, 1) = 3 + 2*sin(2*x);
 Ia(i, 1) = 2 + 2*sin(20*x);
 Ra(i, 1) = 0;
 Sus(i, 1) = 0;

 K1(i,1) = 0.0515 + 0.0515*0.4*cos(x+pi);
 K2(i,1) = 0.061 + 0.061*0.5*cos(x+pi);
 K12(i,1) = 0.0034 + 0.0034*0.3*cos(x+pi);
 K21(i,1) =  0.00064 + 0.00064*0.5*cos(x+pi);
end

while 1
    for n = 1:N 
       for i = 2:M 
    
                LaplaceSu = (Su(i-1,n) - 2*Su(i, n) + Su(i+1, n))/(dx*dx); 
                LaplaceIu = (Iu(i-1,n) - 2*Iu(i, n) + Iu(i+1, n))/(dx*dx);
                LaplaceRu = (Ru(i-1,n) - 2*Ru(i, n) + Ru(i+1, n))/(dx*dx);
                LaplaceLSu = (LSu(i-1,n) - 2*LSu(i, n) + LSu(i+1, n))/(dx*dx);
                LaplaceLIu = (LIu(i-1,n) - 2*LIu(i, n) + LIu(i+1, n))/(dx*dx);
                LaplaceLRu = (LRu(i-1,n) - 2*LRu(i, n) + LRu(i+1, n))/(dx*dx);
                LaplaceSa = (Sa(i-1,n) - 2*Sa(i, n) + Sa(i+1, n))/(dx*dx);
                LaplaceIa = (Ia(i-1,n) - 2*Ia(i, n) + Ia(i+1, n))/(dx*dx);
                LaplaceRa = (Ra(i-1,n) - 2*Ra(i, n) + Ra(i+1, n))/(dx*dx);
                LaplaceSus = (Sus(i-1,n) - 2*Sus(i, n) + Sus(i+1, n))/(dx*dx);
    
                dSu1 = Du*LaplaceSu + LA1 - K1(i,1)*Su(i,n)*Iu(i,n) - K21(i,1)*Su(i,n)*Ia(i,n) - gamma1*Su(i,n) + alpha2*Ru(i,n) + gamma2(i,n)*LSu(i,n) - d1*Su(i,n) - alpha3(i,n)*Su(i,n) - gamma3*Su(i,n) - u6(i,n)*Su(i,n) + gamma4*Sus(i,n);
                dIu1 = Du*LaplaceIu + K1(i,1)*Su(i,n)*Iu(i,n) + K21(i,1)*Su(i,n)*Ia(i,n) - gamma1*Iu(i,n) + gamma2(i,n)*LIu(i,n) - alpha1*Iu(i,n) -u4(i,n)*Iu(i,n) - d1*Iu(i,n);
                dRu1 = Du*LaplaceRu + alpha1*Iu(i,n) + u4(i,n)*Iu(i,n) - alpha2*Ru(i,n) - gamma1*Ru(i,n) + gamma2(i,n)*LRu(i,n) - d1*Ru(i,n) + alpha3(i,n)*Su(i,n);
                dLSu1 = Du*LaplaceLSu + gamma1*Su(i,n) - gamma2(i,n)*LSu(i,n) - d1*LSu(i,n);
                dLIu1 = Du*LaplaceLIu + gamma1*Iu(i,n) - gamma2(i,n)*LIu(i,n) - d1*LIu(i,n);
                dLRu1 = Du*LaplaceLRu + gamma1*Ru(i,n) - gamma2(i,n)*LRu(i,n) - d1*LRu(i,n);
                dSa1 = Da*LaplaceSa + LA2 - K2(i,1)*Sa(i,n)*Ia(i,n) - K12(i,1)*Sa(i,n)*Iu(i,n) + beta2*Ra(i,n) - d2*Sa(i,n) - beta3(i,n)*Sa(i,n);
                dIa1 = Da*LaplaceIa + K2(i,1)*Sa(i,n)*Ia(i,n) + K12(i,1)*Sa(i,n)*Iu(i,n) - beta1*Ia(i,n) -u5(i,n)*Ia(i,n) - d2*Ia(i,n);
                dRa1 = Da*LaplaceRa + beta1*Ia(i,n) +u5(i,n)*Ia(i,n) - beta2*Ra(i,n) - d2*Ra(i,n) + beta3(i,n)*Sa(i,n);
                dSus1 = Du*LaplaceSus + gamma3*Su(i,n) + u6(i,n)*Su(i,n) -gamma4*Sus(i,n) -d1*Sus(i,n);
    
                dSu2 = Du*LaplaceSu + LA1 - K1(i,1)*(Su(i,n)+dSu1*dt/2)*(Iu(i,n)+dIu1*dt/2) - K21(i,1)*(Su(i,n)+dSu1*dt/2)*(Ia(i,n)+dIa1*dt/2) - gamma1*(Su(i,n)+dSu1*dt/2) + alpha2*(Ru(i,n)+dRu1*dt/2) + gamma2(i,n)*(LSu(i,n)+dLSu1*dt/2) - d1*(Su(i,n)+dSu1*dt/2) - alpha3(i,n)*(Su(i,n)+dSu1*dt/2) - gamma3*(Su(i,n)+dSu1*dt/2)- u6(i,n)*(Su(i,n)+dSu1*dt/2) + gamma4*(Sus(i,n)++dSus1*dt/2);
                dIu2 = Du*LaplaceIu + K1(i,1)*(Su(i,n)+dSu1*dt/2)*(Iu(i,n)+dIu1*dt/2) + K21(i,1)*(Su(i,n)+dSu1*dt/2)*(Ia(i,n)+dIa1*dt/2) - gamma1*(Iu(i,n)+dIu1*dt/2) + gamma2(i,n)*(LIu(i,n)+dLIu1*dt/2) - alpha1*(Iu(i,n)+dIu1*dt/2) -u4(i,n)*(Iu(i,n)+dIu1*dt/2)- d1*(Iu(i,n)+dIu1*dt/2);
                dRu2 = Du*LaplaceRu + alpha1*(Iu(i,n)+dIu1*dt/2) + u4(i,n)*(Iu(i,n)+dIu1*dt/2) - alpha2*(Ru(i,n)+dRu1*dt/2) - gamma1*(Ru(i,n)+dRu1*dt/2) + gamma2(i,n)*(LRu(i,n)+dLRu1*dt/2) - d1*(Ru(i,n)+dRu1*dt/2) + alpha3(i,n)*(Su(i,n)+dSu1*dt/2);
                dLSu2 = Du*LaplaceLSu + gamma1*(Su(i,n)+dSu1*dt/2) - gamma2(i,n)*(LSu(i,n)+dLSu1*dt/2) - d1*(LSu(i,n)+dLSu1*dt/2);
                dLIu2 = Du*LaplaceLIu + gamma1*(Iu(i,n)+dIu1*dt/2) - gamma2(i,n)*(LIu(i,n)+dLIu1*dt/2) - d1*(LIu(i,n)+dLIu1*dt/2);
                dLRu2 = Du*LaplaceLRu + gamma1*(Ru(i,n)+dRu1*dt/2) - gamma2(i,n)*(LRu(i,n)+dLRu1*dt/2) - d1*(LRu(i,n)+dLRu1*dt/2);
                dSa2 = Da*LaplaceSa + LA2 - K2(i,1)*(Sa(i,n)+dSa1*dt/2)*(Ia(i,n)+dIa1*dt/2) - K12(i,1)*(Sa(i,n)+dSa1*dt/2)*(Iu(i,n)+dIu1*dt/2) + beta2*(Ra(i,n)+dRa1*dt/2) - d2*(Sa(i,n)+dSa1*dt/2) - beta3(i,n)*(Sa(i,n)+dSa1*dt/2);
                dIa2 = Da*LaplaceIa + K2(i,1)*(Sa(i,n)+dSa1*dt/2)*(Ia(i,n)+dIa1*dt/2) + K12(i,1)*(Sa(i,n)+dSa1*dt/2)*(Iu(i,n)+dIu1*dt/2) - beta1*(Ia(i,n)+dIa1*dt/2) - u5(i,n)*(Ia(i,n)+dIa1*dt/2) - d2*(Ia(i,n)+dIa1*dt/2);
                dRa2 = Da*LaplaceRa + beta1*(Ia(i,n)+dIa1*dt/2) + u5(i,n)*(Ia(i,n)+dIa1*dt/2) - beta2*(Ra(i,n)+dRa1*dt/2) - d2*(Ra(i,n)+dRa1*dt/2) + beta3(i,n)*(Sa(i,n)+dSa1*dt/2);
                dSus2 = Du*LaplaceSus + gamma3*(Su(i,n) + dSu1*dt/2) + u6(i,n)*(Su(i,n) + dSu1*dt/2) -gamma4*(Sus(i,n)+dSus1*dt/2)-d1*(Sus(i,n)+dSus1*dt/2);
    
                dSu3 = Du*LaplaceSu + LA1 - K1(i,1)*(Su(i,n)+dSu2*dt/2)*(Iu(i,n)+dIu2*dt/2) - K21(i,1)*(Su(i,n)+dSu2*dt/2)*(Ia(i,n)+dIa2*dt/2) - gamma1*(Su(i,n)+dSu2*dt/2) + alpha2*(Ru(i,n)+dRu2*dt/2) + gamma2(i,n)*(LSu(i,n)+dLSu2*dt/2) - d1*(Su(i,n)+dSu2*dt/2) - alpha3(i,n)*(Su(i,n)+dSu2*dt/2)- gamma3*(Su(i,n)+dSu2*dt/2)- u6(i,n)*(Su(i,n)+dSu2*dt/2) + gamma4*(Sus(i,n)++dSus2*dt/2);
                dIu3 = Du*LaplaceIu + K1(i,1)*(Su(i,n)+dSu2*dt/2)*(Iu(i,n)+dIu2*dt/2) + K21(i,1)*(Su(i,n)+dSu2*dt/2)*(Ia(i,n)+dIa2*dt/2) - gamma1*(Iu(i,n)+dIu2*dt/2) + gamma2(i,n)*(LIu(i,n)+dLIu2*dt/2) - alpha1*(Iu(i,n)+dIu2*dt/2)  -u4(i,n)*(Iu(i,n)+dIu2*dt/2)- d1*(Iu(i,n)+dIu2*dt/2);
                dRu3 = Du*LaplaceRu + alpha1*(Iu(i,n)+dIu2*dt/2) + u4(i,n)*(Iu(i,n)+dIu2*dt/2) - alpha2*(Ru(i,n)+dRu2*dt/2) - gamma1*(Ru(i,n)+dRu2*dt/2) + gamma2(i,n)*(LRu(i,n)+dLRu2*dt/2) - d1*(Ru(i,n)+dRu2*dt/2) + alpha3(i,n)*(Su(i,n)+dSu2*dt/2);
                dLSu3 = Du*LaplaceLSu + gamma1*(Su(i,n)+dSu2*dt/2) - gamma2(i,n)*(LSu(i,n)+dLSu2*dt/2) - d1*(LSu(i,n)+dLSu2*dt/2);
                dLIu3 = Du*LaplaceLIu + gamma1*(Iu(i,n)+dIu2*dt/2) - gamma2(i,n)*(LIu(i,n)+dLIu2*dt/2) - d1*(LIu(i,n)+dLIu2*dt/2);
                dLRu3 = Du*LaplaceLRu + gamma1*(Ru(i,n)+dRu2*dt/2) - gamma2(i,n)*(LRu(i,n)+dLRu2*dt/2) - d1*(LRu(i,n)+dLRu2*dt/2);
                dSa3 = Da*LaplaceSa + LA2 - K2(i,1)*(Sa(i,n)+dSa2*dt/2)*(Ia(i,n)+dIa2*dt/2) - K12(i,1)*(Sa(i,n)+dSa2*dt/2)*(Iu(i,n)+dIu2*dt/2) + beta2*(Ra(i,n)+dRa2*dt/2) - d2*(Sa(i,n)+dSa2*dt/2) - beta3(i,n)*(Sa(i,n)+dSa2*dt/2);
                dIa3 = Da*LaplaceIa + K2(i,1)*(Sa(i,n)+dSa2*dt/2)*(Ia(i,n)+dIa2*dt/2) + K12(i,1)*(Sa(i,n)+dSa2*dt/2)*(Iu(i,n)+dIu2*dt/2) - beta1*(Ia(i,n)+dIa2*dt/2) - u5(i,n)*(Ia(i,n)+dIa2*dt/2) - d2*(Ia(i,n)+dIa2*dt/2);
                dRa3 = Da*LaplaceRa + beta1*(Ia(i,n)+dIa2*dt/2) + u5(i,n)*(Ia(i,n)+dIa2*dt/2) - beta2*(Ra(i,n)+dRa2*dt/2) - d2*(Ra(i,n)+dRa2*dt/2) + beta3(i,n)*(Sa(i,n)+dSa2*dt/2);
                dSus3 = Du*LaplaceSus + gamma3*(Su(i,n) + dSu2*dt/2) + u6(i,n)*(Su(i,n) + dSu2*dt/2) -gamma4*(Sus(i,n)+dSus2*dt/2)-d1*(Sus(i,n)+dSus2*dt/2);
    
                dSu = (dSu1 + 2*dSu2 + 2*dSu3 + Du*LaplaceSu + LA1 - K1(i,1)*(Su(i,n)+dSu3*dt)*(Iu(i,n)+dIu3*dt) - K21(i,1)*(Su(i,n)+dSu3*dt)*(Ia(i,n)+dIa3*dt) - gamma1*(Su(i,n)+dSu3*dt) + alpha2*(Ru(i,n)+dRu3*dt) + gamma2(i,n)*(LSu(i,n)+dLSu3*dt) - d1*(Su(i,n)+dSu3*dt) - alpha3(i,n)*(Su(i,n)+dSu3*dt) - gamma3*(Su(i,n)+dSu3*dt)- u6(i,n)*(Su(i,n)+dSu3*dt) + gamma4*(Sus(i,n)++dSus3*dt))/6;
                dIu = (dIu1 + 2*dIu2 + 2*dIu3 + Du*LaplaceIu + K1(i,1)*(Su(i,n)+dSu3*dt)*(Iu(i,n)+dIu3*dt) + K21(i,1)*(Su(i,n)+dSu3*dt)*(Ia(i,n)+dIa3*dt) - gamma1*(Iu(i,n)+dIu3*dt) + gamma2(i,n)*(LIu(i,n)+dLIu3*dt) - alpha1*(Iu(i,n)+dIu3*dt) -u4(i,n)*(Iu(i,n)+dIu3*dt) - d1*(Iu(i,n)+dIu3*dt))/6;
                dRu = (dRu1 + 2*dRu2 + 2*dRu3 + Du*LaplaceRu + alpha1*(Iu(i,n)+dIu3*dt) + u4(i,n)*(Iu(i,n)+dIu3*dt) - alpha2*(Ru(i,n)+dRu3*dt) - gamma1*(Ru(i,n)+dRu3*dt) + gamma2(i,n)*(LRu(i,n)+dLRu3*dt) - d1*(Ru(i,n)+dRu3*dt) + alpha3(i,n)*(Su(i,n)+dSu3*dt))/6;
                dLSu = (dLSu1 + 2*dLSu2 + 2*dLSu3 + Du*LaplaceLSu + gamma1*(Su(i,n)+dSu3*dt) - gamma2(i,n)*(LSu(i,n)+dLSu3*dt) - d1*(LSu(i,n)+dLSu3*dt))/6;
                dLIu = (dLIu1 + 2*dLIu2 + 2*dLIu3 + Du*LaplaceLIu + gamma1*(Iu(i,n)+dIu3*dt) - gamma2(i,n)*(LIu(i,n)+dLIu3*dt) - d1*(LIu(i,n)+dLIu3*dt))/6;
                dLRu = (dLRu1 + 2*dLRu2 + 2*dLRu3 + Du*LaplaceLRu + gamma1*(Ru(i,n)+dRu3*dt) - gamma2(i,n)*(LRu(i,n)+dLRu3*dt) - d1*(LRu(i,n)+dLRu3*dt))/6;
                dSa = (dSa1 + 2*dSa2 + 2*dSa3 + Da*LaplaceSa + LA2 - K2(i,1)*(Sa(i,n)+dSa3*dt)*(Ia(i,n)+dIa3*dt) - K12(i,1)*(Sa(i,n)+dSa3*dt)*(Iu(i,n)+dIu3*dt) + beta2*(Ra(i,n)+dRa3*dt) - d2*(Sa(i,n)+dSa3*dt) - beta3(i,n)*(Sa(i,n)+dSa3*dt))/6;
                dIa = (dIa1 + 2*dIa2 + 2*dIa3 + Da*LaplaceIa + K2(i,1)*(Sa(i,n)+dSa3*dt)*(Ia(i,n)+dIa3*dt) + K12(i,1)*(Sa(i,n)+dSa3*dt)*(Iu(i,n)+dIu3*dt) - beta1*(Ia(i,n)+dIa3*dt) - u5(i,n)*(Ia(i,n)+dIa3*dt) - d2*(Ia(i,n)+dIa3*dt))/6;
                dRa = (dRa1 + 2*dRa2 + 2*dRa3 + Da*LaplaceRa + beta1*(Ia(i,n)+dIa3*dt) + u5(i,n)*(Ia(i,n)+dIa3*dt) - beta2*(Ra(i,n)+dRa3*dt) - d2*(Ra(i,n)+dRa3*dt) + beta3(i,n)*(Sa(i,n)+dSa3*dt))/6;
                dSus = (dSus1 + 2*dSus2 + 2*dSus3 + Du*LaplaceSus + gamma3*(Su(i,n) + dSu3*dt) + u6(i,n)*(Su(i,n) + dSu3*dt) -gamma4*(Sus(i,n)+dSus3*dt)-d1*(Sus(i,n)+dSus3*dt))/6;
    
                Su(i, n+1) = Su(i, n) + dt*dSu;
                Iu(i, n+1) = Iu(i, n) + dt*dIu;
                Ru(i, n+1) = Ru(i, n) + dt*dRu;
                LSu(i, n+1) = LSu(i, n) + dt*dLSu;
                LIu(i, n+1) = LIu(i, n) + dt*dLIu;
                LRu(i, n+1) = LRu(i, n) + dt*dLRu;
                Sa(i, n+1) = Sa(i, n) + dt*dSa;
                Ia(i, n+1) = Ia(i, n) + dt*dIa;
                Ra(i, n+1) = Ra(i, n) + dt*dRa;
                Sus(i, n+1) = Sus(i, n) + dt*dSus;
         
    
                if i == 2
                    Su(i-1,n+1) = Su(i,n+1);
                    Iu(i-1,n+1) = Iu(i,n+1);
                    Ru(i-1,n+1) = Ru(i,n+1);
                    LSu(i-1,n+1) = LSu(i,n+1);
                    LIu(i-1,n+1) = LIu(i,n+1);
                    LRu(i-1,n+1) = LRu(i,n+1);
                    Sa(i-1,n+1) = Sa(i,n+1);
                    Ia(i-1,n+1) = Ia(i,n+1);
                    Ra(i-1,n+1) = Ra(i,n+1);
                    Sus(i-1,n+1) = Sus(i,n+1);
                end
                
                if i == M
                    Su(i+1,n+1) = Su(i,n+1);
                    Iu(i+1,n+1) = Iu(i,n+1);
                    Ru(i+1,n+1) = Ru(i,n+1);
                    LSu(i+1,n+1) = LSu(i,n+1);
                    LIu(i+1,n+1) = LIu(i,n+1);
                    LRu(i+1,n+1) = LRu(i,n+1);
                    Sa(i+1,n+1) = Sa(i,n+1);
                    Ia(i+1,n+1) = Ia(i,n+1);
                    Ra(i+1,n+1) = Ra(i,n+1);
                    Sus(i+1,n+1) = Sus(i,n+1);
                end
       end
    end
    
    Iu_sum = 0;
    LSu_sum = 0;
    LIu_sum = 0;
    LRu_sum = 0;
    Ia_sum = 0;
    
    Iu_T_sum = 0;
    LSu_T_sum = 0;
    LIu_T_sum = 0;
    LRu_T_sum = 0;
    Ia_T_sum = 0;
    
    u1_sum = 0;
    u2_sum = 0;
    u3_sum = 0;
    u4_sum = 0;
    u5_sum = 0;
    u6_sum = 0;
    %%%%%%
    for n = 1:N+1
        for i = 1:M+1
            Iu_sum = Iu_sum + (A1 * (Iu(i,n) - Iud)^2)*dt*dx;
            Ia_sum = Ia_sum + (A2 * (Ia(i,n) - Iad)^2)*dt*dx;
            LSu_sum = LSu_sum + (A3*(LSu(i,n) - LSud)^2)*dt*dx;
            LIu_sum = LIu_sum + (A4*(LIu(i,n) - LIud)^2)*dt*dx;
            LRu_sum = LRu_sum +(A5*(LRu(i,n) - LRud)^2)*dt*dx;
            if n == N+1
              Iu_T_sum = Iu_T_sum + (B1*(Iu(i,n) - Iud)^2)*dx;
              Ia_T_sum = Ia_T_sum + (B2*(Ia(i,n) - Iad)^2)*dx;
              LSu_T_sum =LSu_T_sum + (B3*(LSu(i,n) - LSud)^2)*dx;
              LIu_T_sum =LIu_T_sum + (B4*(LIu(i,n) - LIud)^2)*dx;
              LRu_T_sum =LRu_T_sum +  (B5*(LRu(i,n) - LRud)^2)*dx;
            end
        end
    end

    Iu_sum = roundn(Iu_sum,-2);
    Ia_sum = roundn(Ia_sum,-2);
    LSu_sum=roundn(LSu_sum,-2);
    LIu_sum=roundn(LIu_sum,-2);
    LRu_sum=roundn(LRu_sum,-2);
    Iu_T_sum = roundn(Iu_T_sum,-2);
    Ia_T_sum = roundn(Ia_T_sum,-2);
    LSu_T_sum=roundn(LSu_T_sum,-2);
    LIu_T_sum=roundn(LIu_T_sum,-2);
    LRu_T_sum=roundn(LRu_T_sum,-2);

    %%%cost function
    J_sum = Iu_sum + Ia_sum + LSu_sum + LIu_sum + LRu_sum + Iu_T_sum + Ia_T_sum + LSu_T_sum + LIu_T_sum + LRu_T_sum;


    break;
end

Iu_t = zeros(1,N+1);
Ia_t = zeros(1,N+1);
LSu_t = zeros(1,N+1);
LIu_t = zeros(1,N+1);
LRu_t = zeros(1,N+1);
LN_t = zeros(1,N+1);
u1_t = zeros(1,N+1);
u2_t = zeros(1,N+1);
u3_t = zeros(1,N+1);
u4_t = zeros(1,N+1);
u5_t = zeros(1,N+1);
u6_t = zeros(1,N+1);

for n = 1:N+1
    Iu_t(n) = sum(Iu(:,n)*dx);
    Ia_t(n) = sum(Ia(:,n)*dx);
    LSu_t(n) = sum(LSu(:,n)*dx);
    LIu_t(n) = sum(LIu(:,n)*dx);
    LRu_t(n) = sum(LRu(:,n)*dx);
    LN_t(n) = LSu_t(n) + LIu_t(n) + LRu_t(n);
end

Iu_t = Iu_t*100;
Ia_t = Ia_t*100;
LN_t = LN_t*100;

figure(1);
%%%%
set(gcf,'position',[250 300 560 140]);
subplot(1,3,1);
pcolor(Nt,Nx,Iu);
shading flat;
xlabel('$t$','interpreter','latex');
ylabel('$x$','interpreter','latex');
title('$I_{U}$','interpreter','latex');
colorbar;
subplot(1,3,2);
pcolor(Nt,Nx,Ia);
shading flat;
xlabel('$t$','interpreter','latex');
ylabel('$x$','interpreter','latex');
title('$I_{A}$','interpreter','latex');
colorbar;
subplot(1,3,3);
pcolor(Nt,Nx,LSu+LIu+LRu);
shading flat;
xlabel('$t$','interpreter','latex');
ylabel('$x$','interpreter','latex');
title('$L_{N}$','interpreter','latex');
colorbar;

figure(2);
set(gcf,'position',[250 300 560 140]);
subplot(1,3,1);
plot(Nt1,Iu_t);
%%%%
hold on;
plot(Nt1(1),Iu_t(1),'r.','MarkerSize', 10);
plot(Nt1(N+1),Iu_t(N+1),'r.','MarkerSize', 10);
Iu_t(1)=round(Iu_t(1));
Iu_t(N+1)=round(Iu_t(N+1));
sIu0=sprintf('(%.0f,%.0f)',Nt1(1),Iu_t(1));
sIutf=sprintf('(%.0f,%.0f)',Nt1(N+1),Iu_t(N+1));
text(Nt1(1),Iu_t(1),sIu0,'FontSize',8);
text(Nt1(N+1),Iu_t(N+1),sIutf,'FontSize',8);
%%%%%%
xlabel('$t$','interpreter','latex');
ylabel('$number$','interpreter','latex');
title('$I_{U}$','interpreter','latex');

subplot(1,3,2);
plot(Nt1,Ia_t);
%%%%
hold on;
plot(Nt1(1),Ia_t(1),'r.','MarkerSize', 10);
plot(Nt1(N+1),Ia_t(N+1),'r.','MarkerSize', 10);
Ia_t(1)=round(Ia_t(1));
Ia_t(N+1)=ceil(Ia_t(N+1));
sIa0=sprintf('(%.0f,%.0f)',Nt1(1),Ia_t(1));
sIatf=sprintf('(%.0f,%.0f)',Nt1(N+1),Ia_t(N+1));
text(Nt1(1),Ia_t(1),sIa0,'FontSize',8);
text(Nt1(N+1),Ia_t(N+1),sIatf,'FontSize',8);
%%%%%
xlabel('$t$','interpreter','latex');
ylabel('$number$','interpreter','latex');
title('$I_{A}$','interpreter','latex');
subplot(1,3,3);
plot(Nt1,LN_t);
%%%%%
hold on;
plot(Nt1(1),LN_t(1),'r.','MarkerSize', 10);
plot(Nt1(N+1),LN_t(N+1),'r.','MarkerSize', 10);
LN_t(1)=round(LN_t(1));
LN_t(N+1)=round(LN_t(N+1));
sLN0=sprintf('(%.0f,%.0f)',Nt1(1),LN_t(1));
sLNtf=sprintf('(%.0f,%.0f)',Nt1(N+1),LN_t(N+1));
text(Nt1(1),LN_t(1),sLN0,'FontSize',8);
text(Nt1(N+1),LN_t(N+1),sLNtf,'FontSize',8);
%%%%%
xlabel('$t$','interpreter','latex');
ylabel('$number$','interpreter','latex');
title('$L_{N}$','interpreter','latex');