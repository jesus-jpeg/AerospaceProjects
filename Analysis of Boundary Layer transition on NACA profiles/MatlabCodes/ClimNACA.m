function ClimNACA
clear all; close all; clc;
%Parámetros perfil NACA%
t=0.09; %t/c%
h=0; %h/c%
p=0; %p/c%
alpha=0*pi/180; %ángulo ataque%

[Next, xiext, Xext, Yext, cpext, Nint, xiint, Xint, Yint, cpint]=Mpaneles(h,p,t,alpha);

figure
plot(Xint,Yint,'b','LineWidth',2)
hold on
plot(Xext,Yext,'b','LineWidth',2)
axis([0 1 -0.5 0.5])

figure
plot(Xint(1:Nint),cpint(1:Nint),'b','LineWidth',2)
hold on
plot(Xext(1:Next),cpext(1:Next),'g','LineWidth',2)
xlabel('X/c')
ylabel('cp')
legend('cp_{int}','cp_{ext}')
axis([-0.1 1.1 -1 1])

[Xint(1:Nint)' Yint(1:Nint)' cpint(1:Nint)' xiint(1:Nint)']
[Xext(1:Next)' Yext(1:Next)' cpext(1:Next)' xiext(1:Next)']

function [Next,xiext,Xext,Yext,cpext,Nint,xiint,Xint,Yint,cpint]=Mpaneles(h,p,t,alpha)

if p==0
    p=10^-12
end

m=700;
% Número de paneles:
n = 2*m-2;
Xp(1:m) = 1-cos(pi/2*((1:m)-1)/(m-1));

% Ecuaciones perfil NACA:
zE= t/0.2*(0.2969*sqrt(Xp)-0.126*Xp-0.3516*Xp.^2+0.2843*Xp.^3-0.1036*Xp.^4);
for i=1:m,
    if Xp(i)<=p,
        zc(i)=h/p^2*(2*p*Xp(i)-Xp(i)^2);
        dzcdx(i)=h/p^2*(2*p-2*Xp(i));
    else
        zc(i)=h/(1-p)^2*(1-2*p+2*p*Xp(i)-Xp(i)^2);
        dzcdx(i)=h/(1-p)^2*(2*p-2*Xp(i));
    end
end
theta(1:m)=atan(dzcdx);

% Definición del extradós e intradós
xe(1:m)=Xp-zE.*sin(theta);
ze(1:m)=zc+zE.*cos(theta);
xi(1:m)=Xp+zE.*sin(theta);
zi(1:m)=zc-zE.*cos(theta);

% Vértices de paneles:
xv(1:m-1)=xe(1:m-1)
xv(m:n)=xi(m:-1:2);
zv(1:m-1)=ze(1:m-1);
zv(m:n)=zi(m:-1:2);
xv(n+1)=xv(1); zv(n+1)=zv(1);

% Vectores normales y tangentes
difxv(1:n)=xv(2:n+1)-xv(1:n);
difzv(1:n)=zv(2:n+1)-zv(1:n);
D(1:n)=sqrt(difxv.^2+difzv.^2);

vn_x=difzv./D;
vn_z=-difxv./D;

vt_x=difxv./D;
vt_z=difzv./D;

% Puntos medios de paneles
xm(1:n)=(xv(2:n+1)+xv(1:n))/2;
zm(1:n)=(zv(2:n+1)+zv(1:n))/2;
Dm(1:n-1)=sqrt((xm(2:n)-xm(1:n-1)).^2+(zm(2:n)-zm(1:n-1)).^2);
Dm(n)=sqrt((xm(n)-xm(1))^2+(zm(n)-zm(1))^2);

% Matriz método paneles
for j=1:n,
    bjkn=(xv(1:n)-xm(j)).*vn_x+(zv(1:n)-zm(j)).*vn_z;
    bjkt=(xv(1:n)-xm(j)).*vt_x+(zv(1:n)-zm(j)).*vt_z;
    I(j,1:n)=(atan((bjkt+D)./bjkn)-atan(bjkt./bjkn))./pi;
    I(j,j)=-1;
end

xbs=1; zbs=0; 
vne_x=sin(alpha); 
vne_z=-cos(alpha);
vte_x=cos(alpha); 
vte_z=sin(alpha);

for j=1:n,
bjen= (xbs-xm(j))*vne_x+(zbs-zm(j)).*vne_z;
bjet= (xbs-xm(j))*vte_x+(zbs-zm(j)).*vte_z;
I(j,n+1)=(pi/2*sign(bjen)-atan(bjet/bjen))/pi;
end
b(1:n)=2*(xm(1:n)*cos(alpha)+zm(1:n)*sin(alpha));

% Condición de Kutta
I(n+1,1:(n+1))=0; b(n+1)=0;
I(n+1,m-2)=-1/Dm(m-2);
I(n+1,m-1)=1/Dm(m-2);
I(n+1,m)=-1/Dm(m);
I(n+1,m+1)=1/Dm(m);
% Cálculo de potenciales
phi= (I\b')';
Gamma=phi(n+1);
Cl=2*Gamma;

% cp en extradós e intradós geométrico
v(1)=(phi(1)-phi(n))/Dm(n);
cp(1)=1-v(1)*v(1);
v(2:n)=(phi(2:n)-phi(1:n-1))./Dm(1:n-1);
cp(2:n)=1-v(2:n).^2;
cpe(1:m-1)=cp(1:+m-1);
cpi(2:m-1)=cp(n:-1:(m+1));
cpi(1)=cpe(1); cpi(m)=1; cpe(m)=1;

% Cálculos de cp para las corrientes de extradós e intradós 
%(ambas a partir del punto de remanso en la parte frontal) .
maxcpe=-2; 
for i=1:m-1,
    if cpe(i)> maxcpe,
        maxcpe=cpe(i); imaxcpe=i;
        xmaxcpe(i)=xe(i);
    end
end
maxcpi=-2; 
for i=1:m-1,
    if cpi(i)> maxcpi,
        maxcpi=cpi(i); imaxcpi=i;
xmaxcpi(i)=xi(i);
    end
end
if maxcpe >= maxcpi,
    X_PR=xe(imaxcpe);
    Xint(1:imaxcpe)=xe(imaxcpe:-1:1);
    Yint(1:imaxcpe)=ze(imaxcpe:-1:1);
    Xint((imaxcpe+1):(imaxcpe+m-1))=xi(2:m);  
    Yint((imaxcpe+1):(imaxcpe+m-1))=zi(2:m);
    cpint(1:imaxcpe)=cpe(imaxcpe:-1:1);
    cpint((imaxcpe+1):(imaxcpe+m-1))=cpi(2:m);
    Xext(1:(m-imaxcpe+1))=xe(imaxcpe:m);
    Yext(1:(m-imaxcpe+1))=ze(imaxcpe:m);
    cpext(1:(m-imaxcpe+1))=cpe(imaxcpe:m);
end
if maxcpi>= maxcpe,
    X_PR=xi(imaxcpi);
    Xext(1:imaxcpi)=xi(imaxcpi:-1:1);
    Yext(1:imaxcpi)=zi(imaxcpi:-1:1);
    Xext((imaxcpi+1):(imaxcpi+m-1))=xe(2:m);
    Yext((imaxcpi+1):(imaxcpi+m-1))=ze(2:m);
    cpext(1:imaxcpi)=cpi(imaxcpi:-1:1);
    cpext((imaxcpi+1):(imaxcpi+m-1))=cpe(2:m);
    Xint(1:(m-imaxcpi+1))=xi(imaxcpi:m);
    Yint(1:(m-imaxcpi+1))=zi(imaxcpi:m);
    cpint(1:(m-imaxcpi+1))=cpi(imaxcpi:m);
end

Nint=length(Xint);
Next=length(Xext);

xiint(1)=0;
for k=2:Nint,
    xiint(k)=xiint(k-1)+sqrt((Xint(k)-Xint(k-1))^2+(Yint(k)-Yint(k-1))^2);
end
xiext(1)=0
for k=2:Next,
    xiext(k)=xiext(k-1)+sqrt((Xext(k)-Xext(k-1))^2+(Yext(k)-Yext(k-1))^2);
end
   
