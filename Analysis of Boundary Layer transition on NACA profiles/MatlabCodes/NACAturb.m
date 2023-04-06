clc; close all; clear all;
%
t=0.09; % es t/c
h=0; % es h/c
p=0;  % es p/c
alpha=1.2*pi/180; % ángulo de ataque

% Número de Reynolds
 Re=[5*10^4 1*10^5 5*10^5 1*10^6 1*10^7]
 for z=1:5 %%bucle
[Next, xiext, Xext, Yext, cpext, Nint, xiint, Xint, Yint, cpint]=Mpaneles(h,p,t,alpha);
iext=1; %cambiar para intradós o extradós%
if iext==1;
    Nxi=5000;
    ximin=0.0001; ximax=xiext(Next); xi=linspace(ximin,ximax,Nxi);
    X(1:Nxi)=interp1(xiext,Xext,xi(1:Nxi));
    Y(1:Nxi)=interp1(xiext,Yext,xi(1:Nxi));
    c_p(1:Nxi)=interp1(xiext,cpext,xi(1:Nxi));
else 
    Nxi=5000;
    ximin=0.0001; ximax=xiint(Nint); xi=linspace(ximin,ximax,Nxi);
    X(1:Nxi)=interp1(xiint,Xint,xi(1:Nxi));
    Y(1:Nxi)=interp1(xiint,Yint,xi(1:Nxi));
    c_p(1:Nxi)=interp1(xiint,cpint,xi(1:Nxi));
end

%Mallado
Neta=5000; etamin=0; etamax=40/sqrt(Re(z));
heta=(etamax-etamin)/(Neta-1);
eta(1:Neta)=etamin+((1:Neta)-1)*heta;

% hxi=(ximax-ximin)/(Nxi-1);
% xi(1:Nxi)=ximin+((1:Nxi)-1)*hxi;

d(1:Nxi-1)=sqrt((X(2:Nxi)-X(1:Nxi-1)).^2+(Y(2:Nxi)-Y(1:Nxi-1)).^2);
tX(1:Nxi-1)=(X(2:Nxi)-X(1:Nxi-1))./d;
tX(Nxi)=tX(Nxi-1);
tY(1:Nxi-1)=(Y(2:Nxi)-Y(1:Nxi-1))./d;
tY(Nxi)=tY(Nxi-1);
if iext==1
    nX=-tY; nY=tX;
else
    nX=tY; nY=-tX;
end

% Mallado en xi
xi(1)=0;
for n=2:Nxi
xi(n)=xi(n-1)+d(n-1);
end
hxi(2:Nxi)=xi(2:Nxi)-xi(1:Nxi-1); %
% Matrices derivadas en eta 
D1=sparse(Neta,Neta); D2=sparse(Neta,Neta);
for j=2:Neta-1
  D1(j,j-1)=-1/2/heta;
D1(j,j+1)=1/2/heta;
  D2(j,j-1)=1/heta^2; 
  D2(j,j)=-2/heta^2; D2(j,j+1)=1/heta^2;
end

% Corriente exterior 
ue(1:Nxi)=sqrt(1-c_p); 
duedxi(2:Nxi)=(ue(2:Nxi)-ue(1:Nxi-1))./hxi(2:Nxi);
%Valores en estacion inicial n=1
% (es la n-1 para n=2) 
unm1(1:Neta)=ue(1); vnm1(1:Neta)=0;
%
% Pintar solución en Nplot estaciones xiplot:
Nplot=7;
Xplot=[0.6,0.63,0.68,0.72,0.75,0.81,0.85]; icont=1;
%Avance en estaciones
%
%%%% Criterio Pohlhausen %%%%%%%%%%%%%%% 
iturb=0; ntr=0; Xtr=10;
Lambda_graf=[ -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
 Pohl_cri=[100 120 140 170 250 360 645 1150 2000 3500 5700 8000 9700 11000 11500];
delta(1)=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=2:Nxi,
     delta1 (n)=trapz ( eta ,1-unm1/ue (n-1)) ;
    Lambda(n)=Re(z)*(delta(n-1)^2)*duedxi(n); 
    Polh(n)=delta1(n)*Re(z)*ue(n);
if Polh(n)>=interp1(Lambda_graf,Pohl_cri,Lambda(n))
        if iturb==0,
           ntr=n; 
           Xtr=X(ntr);
           Ytr=Y(ntr);
        end
iturb=1;
end

if iturb==1 
    uast(n)=sqrt(1/Re(z)*unm1(2)/heta); 
    nuTn(1:Neta)=min(0.41*(Re(z)*eta*uast(n)-10.7*tanh(eta*uast(n)*Re(z)/10.7)),...
        0.018*delta1(n)*ue(n)*Re(z));
else
    nuTn(1:Neta)=0;
end
    f(1:Neta)=1+nuTn;
    fdot(1:Neta)=D1*f';
%
Mv=spdiags(vnm1'-fdot'/Re(z),0,Neta,Neta); Mu=spdiags(unm1',0,Neta,Neta)/hxi(n) ;
    Mf=spdiags(f'/Re(z),0,Neta,Neta);
r=unm1.^2/hxi(n)+ue(n)*duedxi(n); A=-Mf*D2+Mv*D1+Mu;
A(1,:)=0; A(1,1)=1; r(1)=0; A(Neta,:)=0; A(Neta,Neta)=1;
r(Neta)=ue(n);
    un=(A\r')';
    c_f(n)=2/Re(z)*un(2)/heta;

%     Perfil de velocidades
%     figure(1)
%    if X(n)>= Xplot(icont)
%     plot(un,eta)
%     legend('X=0.12','X=0.26','X=0.37','X=0.58','X=0.71','X=0.74','X_s=0.7816')
%     title(['Perfil de velocidades para Re=' num2str(Re)])
%     xlabel('u')
%     ylabel('\eta')
%     axis([0 2 0 etamax])
%  hold on
%  icont=icont+1
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%Ley de la pared:
% if iturb==1
%          figure(3)
%          if X(n)>= Xplot(icont)
%     plot(log10(Re*uast(n)*eta),un/uast(n))
%       legend('X=0.6','X=0.63','X=0.68','X=0.72','X=0.75','X=0.81','X=0.85')
%     title(['Ley de la pared para Re=' num2str(Re)])
%         axis([0 4 0 40]);
%        hold on
% 
%      %plot(log10(Re*uast(n)*eta),Re*uast(n)*eta,'g')
%  
%      %plot(log10(Re*uast(n)*eta),2.5*log(Re*uast(n)*eta)+5,'r')
%    icont=icont+1
%    end
% % % % %
%      end

    if un(2)<= 0
        ns=n;
        break
end

% Espesor de capa límite
% Cálculo de espesor de capa límite
j=1;
      while un(j)<=0.99*ue(n)
           j=j+1;
           delta(n)=eta(j);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Cálculo de vn
    dudxi=(un-unm1)./hxi(n);
    vn(1)=0;
    for j=2:Neta
    vn(j)=vn(j-1)- 0.5*heta*(dudxi(j)+dudxi(j-1));
    end
%
% Actualizar u y v para próxima estación
 unm1=un; vnm1=vn;

end

%Representación coeficiente de fricción en función de Re   
%     figure (4)
%         plot(Re(z).*xi(2:ns),c_f(2:ns),'DisplayName',['Re_c=' num2str(Re(z))])
%    title('Coeficiente de fricción para el extradós') 
%    legend('-DynamicLegend');
%    xlabel('Re_x')
%    ylabel('c_f')
%     axis([0 1.5*10^5 0 max(c_f)])
%     hold all

% Representación puntos sobre el perfil%
figure(2)
plot(Xext,Yext,'b',Xint,Yint,'r','LineWidth',2)
hold on
plot(X(ns),Y(ns),'o','MarkerSize',8,'MarkerFaceColor','[0 1 0]')
plot(Xtr,Ytr,'o','MarkerSize',8,'MarkerFaceColor','[1 0 0]')
legend( 'Extradós','Intradós','punto desprendimiento','punto transición')
hold off
axis([0 1 -0.4 0.4])

nalpha=nX*cos(alpha)+nY*sin(alpha) ; talpha=tX*cos(alpha)+tY*sin(alpha) ;
C_dp=trapz(xi(1:ns),- c_p(1:ns).*nalpha(1:ns))...
- c_p(ns)*trapz(xi(ns:Nxi),nalpha(ns:Nxi)); 
C_df=trapz(xi(2:ns),c_f(2:ns).*talpha(2:ns));
C_d=C_dp+C_df;

% Coeficiciente de resistencia de fricción en función de log10Re
% Recd=[4  4.7, 5, 5.7,  6, 6.7, 7];%Calcular para cada caso%
% Cdf=[0.0226 0.0144 0.0125 0.0094 0.0085 0.0068 0.0062];
% figure (5)
% plot(Recd,Cdf)
% title ([ 'Coeficiente de resistencia de fricción'])
% xlabel('log_{10}Re')
% ylabel('C_{df}')
%  [Xtr X(ns)]

% Caso simétrico
% sigma=2*max(Y(1:Nxi));
 2*[C_dp C_df C_d]
 
 %Pendiente curva de sustentación%
nperp=-sin(alpha)*nX+cos(alpha)*nY; 
C_l=trapz(xi(1:ns),-c_p(1:ns).*nperp(1:ns))-c_p(ns)*trapz(xi(ns:Nxi),nperp(ns:Nxi));
[C_l]
end %%%%bucle
