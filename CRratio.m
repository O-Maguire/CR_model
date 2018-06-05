% Code writen by Oisin Maguire June 2015
% based on CR model developed by Colombant and Tonon:
% "X-ray emission in laser-produced plasma"
% CR_temp_ratio I/O structure is (atomic number, temperature in eV, wavelength in nm)
% returns the ratios of the ion stages including neutral and bear :-)
% Last update July 2015
function Ratios = CRratio(A, temperature, lambda)
Ratios=zeros(A+1,2);
for i=1:A+1
    Ratios(i,1)=i;
end
c=299792458;
me=9.1093897E-31;
esp0=8.854187817E-12;
ec=1.60217733E-19;
lambda = lambda*1e-9;
ne=4*(pi^2)*(c^2)*me*esp0/((ec^2)*(lambda^2))*1E-6;
Te=temperature;
potentialdata=load('Ionisation_5.csv');
IP=zeros(A,1);
for i=1:A
    IP(i)=potentialdata(i,A);
end
edata=load('electron_order2.dat');
outer_e_sn=zeros(A+1,1); i=0;
for fixing=A:-1:1
    i=i+1;
    outer_e_sn(i)=edata(fixing);
end
outer_e_sn(A+1)=0;
n(i,1)=1;
S(i,1)=(((9E-6)*outer_e_sn(1)*((Te/IP(1))^(1/2)))/((IP(1)^(3/2))*(4.88+(Te/IP(1)))))*exp(-IP(1)/Te);
alphaR(i,1)=(5.2E-14)*((IP(1)/Te)^(1/2))*0*(0.429+(0.5*log10(IP(1)/Te))+(0.469*((Te/IP(1))^(1/2))));
alpha3b(i,1)=((2.97E-27)*outer_e_sn(1))/((Te*(IP(1)^2))*(4.88+(Te/IP(1))));
for j=2:length(IP)
    S(i,j)=(((9E-6)*outer_e_sn(j)*((Te/IP(j))^(1/2)))/((IP(j)^(3/2))*(4.88+(Te/IP(j)))))*exp(-IP(j)/Te);
    alphaR(i,j)= (5.2E-14)*((IP(j)/Te)^(1/2))*(j-1)*(0.429+(0.5*log10(IP(j)/Te))+(0.469*((Te/IP(j))^(1/2))));
    alpha3b(i,j)=((2.97E-27)*outer_e_sn(j))/((Te*(IP(j)^2))*(4.88+(Te/IP(j))));
    n(i,j)=(S(i,j-1)/(alphaR(i,j)+ne*alpha3b(i,j)))*n(i,j-1);
end
Sbare=(((9E-6)*1*((Te/IP(A))^(1/2)))/((IP(A)^(3/2))*(4.88+(Te/IP(A)))))*exp(-IP(A)/Te);
alphaRbare=(5.2E-14)*(0.469)*A;
n(i,length(IP)+1)=(n(i,length(IP)).*Sbare)./(alphaRbare);
nt(i)=sum(n(i,:));
frac(i,:)=n(i,:)./nt(i);
Ratios(:,2)=frac(i,:);
end