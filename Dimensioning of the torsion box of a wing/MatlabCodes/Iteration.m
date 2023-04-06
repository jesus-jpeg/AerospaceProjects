a=1000;%mm
b=3700;%mm
Mx=11399300000;%N*mm
y=500;%mm
T=-3102.8*10^6;%N*mm
Sy=1489400;%N
rho=2590;%kg/m^3
g=9.81;%m/s^2

Ixx=(1/12)*a^(3)*(tba+tbs+ti)+(a^(2)*b*th)/2
A=2*b*th+a*(tba+tbs+ti)
A1=a*(b-d)
A2=a*d
d=b/2


sigmaz=Mx*y/Ixx
qs1=-((Sy/Ixx)*a*((a*tbs/8)+th*(b-d)/2))
qs2=-(Sy*a/Ixx)*(((tbs+ti)*a/(8))+b*th/4+th*d/2) %maximo en s1=b-d, s2=d
qt2=(1/(2*G*A1)*(2*(b-d)/th+a/tbs+a/ti)*T/(2*A1)+a/ti*T/(2*A1))/(1/(2*G*A1)*((2*(b-d)/th+a/tbs+a/ti)*2*A2/(2*A1)-a/ti)-1/(2*G*A2)*((2*d/th+a/tba+a/ti)-a/ti*2*A2/2A1))
qt1=(2*A2*qt2-T)/(2*A1)

%T=2*A1*qt1+2*A2*qt2
%(1/(2*G*A1))*((2*(b-d)/th+a/tbs+a/ti)*qt1-(a/ti)*qt2)=(1/(2*G*A2))*((2*d/th+a/tba+a/ti)*qt2-(a/ti)*qt1)
%de estas dos ec sacamos qt1 y qt2
q2=((Sy/Ixx)*(((tbs+ti)*a^2*d)/(4*th)+(b-d)*d+d^2+3*a^3/24+((tbs+ti)*a^3)/(8*tba)+th*b*a^2/2*tba)+(Sy/Ixx)*((tbs*(b-d)*a^2)/(8*th)+(a*(b-d)^2)/4+(a^3)/12)*(a/ti)/(2*(b-d)/th+a/tbs+a/ti))/((2*d/th+a/tba+a/ti)-(a/ti)^2/(2*(b-d)/th+a/tbs+a/ti))
q1=((a/ti)*q2+(Sy/Ixx)*((tbs*(b-d)*a^2)/(8*th)+(a*(b-d)^2)/4+(a^3)/12))/2*(b-d)/th+a/tbs+a/ti)


%(2*(b-d)/th+a/tbs+a/ti)*q1-(a/ti)*q2-(Sy/Ixx)*((tbs*(b-d)*a^2)/(8*th)+(a*(b-d)^2)/4+(a^3)/12)=0
%(2*d/th+a/tba+a/ti)*q2-(a/ti)*q1-(Sy/Ixx)*(((tbs+ti)*a^2*d)/(4*th)+(b-d)*d+d^2+3*a^3/24+((tbs+ti)*a^3)/(8*tba)+th*b*a^2/2*tba)=0

