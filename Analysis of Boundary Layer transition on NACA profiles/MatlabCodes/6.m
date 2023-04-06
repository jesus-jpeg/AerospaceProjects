
clc; close all; clear all;
%
%%*** Número de Reynolds, Pr, Ec, igas ****
Re=7*10^4; Pr=6.7; Ec=0; igas=0; %Cambiamos valores para agua o gas%
%Mallado en eta y xi
etamin=0; 
etamax=10/sqrt(Re); Neta=2000; 
heta=(etamax-etamin)/(Neta-1); 
eta(1:Neta)=etamin+((1:Neta)-1)*heta;
ximin=0.001; 
ximax=pi/2; Nxi=3000; 
hxi=(ximax-ximin)/(Nxi-1); 
xi(1:Nxi)=ximin+((1:Nxi)-1)*hxi;
%
% Matrices derivadas en eta 
D1=sparse(Neta,Neta); 
D2=sparse(Neta,Neta);
for j=2:Neta-1
  D1(j,j-1)=-1/2/heta;
  D1(j,j+1)=1/2/heta;
  D2(j,j-1)=1/heta^2; 
  D2(j,j)=-2/heta^2; 
  D2(j,j+1)=1/heta^2;
end
%
%%%%%%%%CORRIENTE EXTERIOR EXPERIMENTAL%%%%%%%%%%%%%%%%%%%
xiexp=[0 7.45 14.16 23.11 26.83 30.56 38.01 40.99 50 55 60 65 70 75 80 90 100 180]*pi/180/2;
cpexp=[1 0.97 0.82 0.53 0.40 0.28 -0.10 -0.25 -0.9 -1.1 -1.2 -1.3 -1.35 -1.3 -1.16 -1.17 -1.17 -1.17];
for k=1:Nxi, 
    cpe(k)=interp1(xiexp,cpexp,xi(k));
    ue(k)=sqrt(1-cpe(k));
end
duedxi(2:Nxi)=(ue(2:Nxi)-ue(1:Nxi-1))/hxi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if igas==1, 
    Thetae(1:Nxi)=1-Ec/2*(1-ue.^2);
else
    Thetae(1:Nxi)=1;
end
%Valores en estacion inicial n=1
% (es la n-1 para n=2) 
unm1(1:Neta)=ue(1);
vnm1(1:Neta)=0;
%
Thetanm1(1:Neta)=Thetae(1);
% Pintar solución en Nplot estaciones xiplot: 
X(1:Nxi)=1/2*(1-cos(2*xi));
Nplot=7;
Xplot=[0.06 0.14 0.21 0.29 0.35 0.4055 1]; 
icont=1;
%

%Avance en estaciones
for n=2:Nxi,
    % Matrices del sistema y r
    Mv=spdiags(vnm1',0,Neta,Neta);
    Mu=spdiags(unm1',0,Neta,Neta)/hxi; 
    A=-1/Re*D2+Mv*D1+Mu; 
    r=unm1.^2/hxi+ue(n)*duedxi(n);
%
%implementar condiciones de contorno
    A(1,:)=0; A(1,1)=1; r(1)=0;
    A(Neta,:)=0; A(Neta,Neta)=1;
    r(Neta)=ue(n);
% %Solución
    un=(A\r')';
    % Cálculo de cf(xi_n)
    cf(n)=2/Re*un(2)/heta; 

% Cálculo de vn
    dudxi=(un-unm1)/hxi;
    vn(1)=0;
    for j=2:Neta
    vn(j)=vn(j-1)-0.5*heta*(dudxi(j)+dudxi(j-1));
    end
   
%% 
Mat_un=spdiags(un',0,Neta,Neta)/hxi; 
Mat_vn=spdiags(vn',0,Neta,Neta); 
Atheta=Mat_un+Mat_vn*D1-D2/Re/Pr; 
rtheta=un.*Thetanm1/hxi-(D1*un')'.^2*Ec/Re+igas*Ec*un*ue(n)*duedxi(n); 
Atheta(1,:)=0; 
Atheta(1,1)=1; rtheta(1)=0;
Atheta(Neta,:)=0; 
Atheta(Neta,Neta)=1; 
rtheta(Neta)=Thetae(n); 
Thetan=(Atheta\rtheta')'; 
dThetadeta_0(n)=Thetan(2)/heta;
hold on 

%Pintamos perfil de temperaturas
if X(n)>= Xplot(icont) || un(2)<=0 
plot(1+(300-350)/350*Thetan,eta) 
legend('X=0.06','X=0.14','X=0.21','X=0.29','X=0.35','X_s=0.4055')
title('Perfil de temperaturas para el agua')
xlabel('T/T_p')
ylabel('\eta')
axis([0.8 1.1 etamin etamax]); 



% Actualizar u y T para próximaestación
 if un(2)<=0 
    ns=n;
    break 
 end
icont=icont+1
end
unm1=un; vnm1=vn;

 Thetanm1=Thetan;
end

%punto de separación%
[X(ns)]
%Cálculo de Nusselt%
sigma=1
Nu=trapz(xi(2:ns),dThetadeta_0(2:ns)); 
Nu_corr=(0.75-0.16*exp(- 0.018/sigma^3.1))*Re^0.5*Pr^(1/3);
[Nu Nu_corr]
%**********************
%
