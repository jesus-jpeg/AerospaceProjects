clc; close all; clear all;
%
% Número de Reynolds
Re=1.347e6
%
%Mallado en eta y xi
etamin=0;
etamax=40/sqrt(Re);
Neta=2000; heta=(etamax-etamin)/(Neta-1); 
eta(1:Neta)=etamin+((1:Neta)-1)*heta;
ximin=0; ximax=1; Nxi=2000; 
hxi=(ximax-ximin)/(Nxi-1); 
xi(1:Nxi)=ximin+((1:Nxi)-1)*hxi;
% Matrices derivadas 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%Corriente exterior y valores en estacion inicial n=1
% (es la n-1 para n=2)
ue(1:Nxi)=1;
unm1(1:Neta)=ue(1);
vnm1(1:Neta)=0; 
%
%Avance en estaciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%
iturb=0; xitr=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pintar solución en Nplot
%estaciones xiplot:
Nplot=7;
xiplot=[0.24 0.36 0.49 0.67 0.75 0.89 0.99];
icont=1;

for n=2:Nxi,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo de Reichardt:
 delta1(n)=trapz(eta,1-unm1/ue(n-1));
    if Re*delta1(n)>=645
        if iturb==0,
        xitr=xi(n);
    end
    iturb=1;
end 
if iturb==1 
    uast(n)=sqrt(1/Re*unm1(2)/heta); 
    nuTn(1:Neta)=min(0.41*(Re*eta*uast(n)-10.7*tanh(eta*uast(n)*Re/10.7)),0.018*delta1(n)*ue(n)*Re);
else
 nuTn(1:Neta)=0;
end %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(1:Neta)=1+nuTn;
fdot(1:Neta)=D1*f';

%
% Matrices del sistema y r
Mu=spdiags(unm1',0,Neta,Neta)/hxi; 
Mv=spdiags(vnm1'- fdot'/Re,0,Neta,Neta); 
Mf=spdiags(f'/Re,0,Neta,Neta); r=unm1.^2/hxi;
 A=-Mf*D2+Mv*D1+Mu;
%
%%%%%%%%%%%%%%%%%%%%%%%
%
%implementar condiciones de contorno
 A(1,:)=0; A(1,1)=1; r(1)=0;
 A(Neta,:)=0; A(Neta,Neta)=1;
r(Neta)=1;
%
%Solución
    un=(A\r')';
    cf(n)=2/Re*un(2)/heta;

%
% Ley de la pared:
% figure(1)
% if iturb==1
%          
%          if xi(n)>= xiplot(icont)
% plot(log10(Re*uast(n)*eta),un/uast(n))
% legend('\xi=0.24','\xi=0.36','\xi=0.49','\xi=0.67','\xi=0.75','\xi=0.89','\xi=1')
% title(['Ley de la pared para Re=' num2str(Re)])
%         axis([0 4 0 30]);
%        hold on
% 
% %     plot(log10(Re*uast(n)*eta),Re*uast(n)*eta,'g')
% % 
% %     plot(log10(Re*uast(n)*eta),2.5*log(Re*uast(n)*eta)+5,'r')
%    icont=icont+1
%    end
% % % % % % % 
% end
    

%Perfil de velocidades
% figure(2)
%     if xi(n)>= xiplot(icont)
%     plot(un,eta)
%     legend('\xi=0.15','\xi=0.34','\xi=0.52','\xi=0.67','\xi=0.75','\xi=0.89','\xi=1')
%     title(['Perfil de velocidades para Re=' num2str(Re)])
%     axis([0 2 0 0.035])
%  hold on
%   icont=icont+1;
%     end
  
%
% Cálculo de vn
    dudxi=(un-unm1)/hxi;
    vn(1)=0;
    for j=2:Neta
    vn(j)=vn(j-1)-0.5*heta*(dudxi(j)+dudxi(j-1));
end
%
% Actualizar u y v para próxima estacion
 unm1=un; vnm1=vn;
 for i=1:Neta
    if round(un(i),2)==0.99 
        delta(n)=eta(i) ;
    end
 end

end
% C\'alculo de C
% figure(3)
% C=delta.*xi.^(-4/5)*Re^(1/5); 
% plot(xi,C)
% axis([0 1 0 0.5 ]) 
% title('Calculo de C')
% legend (['Re=' num2str(Re)])

%Cálculo de C1
% figure (4)
% loglog(Re.*xi,cf) 
% title('Calculo de C_1 y m') 
% xlabel('Re_x')
% ylabel('c_f')
% axis([1e3 1e7 1e-3 1e-1])
% hold on


%Representación coeficiente de fricción en función de Re   
%     figure (5)
%     plot(Re.*xi,cf)
%    title(['Coeficiente de friccion para Re=' num2str(Re)]) 
%    xlabel('Re_x')
%    ylabel('c_f')
%     axis([0 13e5 0 0.02])

%Comparación coeficientes de fricción  
%     figure (6)
%     plot(Re.*xi,cf)
%     hold on
%     plot(Re.*xi,0.0925./(Re.*xi).^0.23329)
%    title(['Comparación de los coeficientes de fricción']) 
%    xlabel('Re_x')
%    ylabel('c_f')
%    legend('c_f','c_fley')
%     axis([0 14e5 0 0.02])
    
%Representación coeficiente de fricción en función de xi
% figure (7)
% plot(xi,cf)
% title(['Coeficiente de friccion para Re=' num2str(Re)])
% xlabel('\xi')
% ylabel('c_f')
% axis ([0 1 0 0.05])

%Representación espesores de capa límite y de desplazamiento
%adimensionales
% figure (8)
% plot(xi,delta1,'g')
% hold on
% plot(xi , delta , 'b')
% axis ([0 1 0 0.02])
% legend([ '\delta_1/L'] ,[ '\delta/L'])
% title ([ 'Espesor de capa limite y de desplazamiento adimensional para Re_L = ' num2str(Re)])

%punto transición
[xitr]

%C_D frente a log_10Re
% figure(9)
% Recd=[4.1, 4.3, 4.6, 4.8, 4.9, 5.1, 5.2, 5.3, 5.4, 5.6, 5.8, 6.2, 6.3, 6.6, 7.1];
% Cd=[0.0127 0.0074 0.0052 0.0045 0.0042 0.0039 0.0037 0.0041 
%0.0043 0.0042 0.0040 0.00375 0.00355 0.0033 0.0029];
% %Cálculo C_D
% C_D=trapz(xi(2:Nxi),cf(2:Nxi));
% [C_D]
% figure (6)
% semilogx(Recd,Cd)
% xlabel('log_{10}Re_L')
% ylabel('C_D')
% hold on
% title('C_D para una placa plana')