function Ap3

t=0.12; % es t/c
h=0; % es h/c
p=0;  % es p/c
alpha=0*pi/180; % ángulo de ataque
%Número de Reynolds
Re=5*10^5;

[Next, xiext, Xext, Yext, cpext, Nint, xiint, Xint, Yint, cpint]=Mpaneles(h,p,t,alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nxi=5000;
iext=1; %cambiar para intradós o extradós%
if iext==1;
    ximin=0.0001; ximax=xiext(Next); xi=linspace(ximin,ximax,Nxi);
    X(1:Nxi)=interp1(xiext,Xext,xi(1:Nxi));
    Y(1:Nxi)=interp1(xiext,Yext,xi(1:Nxi));
    c_p(1:Nxi)=interp1(xiext,cpext,xi(1:Nxi));
else 
    ximin=0.0001; ximax=xiint(Nint); xi=linspace(ximin,ximax,Nxi);
    X(1:Nxi)=interp1(xiint,Xint,xi(1:Nxi));
    Y(1:Nxi)=interp1(xiint,Yint,xi(1:Nxi));
    c_p(1:Nxi)=interp1(xiint,cpint,xi(1:Nxi));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Mallado:
Neta=5000; etamin=0; etamax=40/sqrt(Re);
heta=(etamax-etamin)/(Neta-1);
eta(1:Neta)=etamin+((1:Neta)-1)*heta;

hxi=(ximax-ximin)/(Nxi-1);
xi(1:Nxi)=ximin+((1:Nxi)-1)*hxi;

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

%Matrices derivadas en eta:
D1=sparse(Neta,Neta);
D2=sparse(Neta,Neta);
for j=2:Neta-1
    D1(j,j-1)=-1/2/heta;
    D1(j,j+1)=1/2/heta;
    D2(j,j-1)=1/heta^2;
    D2(j,j)=-2/heta^2;
    D2(j,j+1)=1/heta^2;
end

%Corriente exterior
ue(1:Nxi)=sqrt(1-c_p(1:Nxi));
duedxi(2:Nxi)=(ue(2:Nxi)-ue(1:Nxi-1))./hxi;

%Valores en estacion inicial n=l
% (es la n-1 para n=2)
unml(1:Neta)=ue(1);
vnml(1:Neta)=0;


Nplot=7;
Xplot=[0.2 0.3 0.5 0.6 0.7 0.8 1];
icont=1;

for n=2:Nxi,

% Matrices del sistema y r
Mv=spdiags(vnml',0,Neta,Neta);
Mu=spdiags(unml',0,Neta,Neta)/hxi;
A=-1/Re*D2+Mv*D1+Mu;
r=unml.^2/hxi+ue(n)*duedxi(n);

%implementación de las condiciones de contorno
A(1, :)=0; A(1,1)=1; r(1)=0;
A(Neta, :)=0; A(Neta,Neta)=1;
r(Neta)=ue(n);
%
%Solución
un=(A\r')';

% Cálculo de cf (xi n)
c_f(n)=2/Re*un(2)/heta;
% Perfil de velocidades
if X(n)>= Xplot(icont) || un(2)<=0 
      plot(un,eta)
     axis ( [0 2 etamin etamax])
     hold on
     if un(2)<=0
         ns=n;
         break
     end
     icont=icont+1;
end

% Cálculo de vn
dudxi=(un-unml)./hxi;
vn(1)=0;
for j=2:Neta
vn(j)=vn(j-1)-...
0.5*heta*(dudxi(j)+dudxi(j-1));
end

% Actualizar u y v para próxima
%estación
unml=un; vnml=vn;
end
%Punto de separacion:
[X(ns)]

nalpha=nX*cos(alpha)+nY*sin(alpha) ; 
talpha=tX*cos(alpha)+tY*sin(alpha) ;
C_dp=trapz(xi(1:ns),- c_p(1:ns).*nalpha(1:ns))-c_p(ns)*trapz(xi(ns:Nxi),...
    nalpha(ns:Nxi));
C_df=trapz(xi(2:ns),c_f(2:ns).*talpha(2:ns)); 
C_d=C_dp+C_df;
% Coeficientes basados en el espesor
sigma=2*max(Y(1:Nxi));
2*[C_dp C_df C_d]/sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
