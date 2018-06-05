% Code written by Oisin Maguire in 2015 \approx Feb 
% Last update July 2015
close all; clear all; clc;
c=299792458; me=9.1093897E-31; esp0=8.854187817E-12; ec=1.60217733E-19;     %Physical constants
options.Interpreter='tex';                                                  %turn on tex for dialog boxes
%get laser wavelength and electrom density values
title='Input laser wavelength and electrom density';
prompt={'Please enter the laser wavelength (in nm):',...
    ['Please enter the fraction of the critical density you wish to use or the density   ';...
    'you wish to use (<=1 for fraction, >1 for actual density (/cm^{-3}))               ']};
def={'1064','1'};
answer=inputdlg(prompt,title,1,def,options);
lambda = str2num(answer{1})*1e-9;                                           %laser wavelength
eden = str2num(answer{2});
if eden <= 1
    frac_ne=eden;
    ne=4*(pi^2)*(c^2)*me*esp0/((ec^2)*(lambda^2))*1E-6*frac_ne;             %calculate critical density and x fraction to be used
else ne=eden;end                                                            %user input electron density
title='Power density or temperature';
prompt={['Do you wish to work run for a specific power density, or work in temperature';...
    'Please enter 1 for temperature, or enter the power density (in Wcm^{-2})    ']};
def={'1'};
answer=inputdlg(prompt,title,1,def,options);
pd = str2num(answer{1});                                                    %power density
prompt1= 'what is the atomic number?';
A=str2num(cell2mat(inputdlg(prompt1)));                                     %atomic number
if pd>1
    Te=5.2e-6*(A^(1/5))*((((lambda*(10^6))^(2))*pd)^(3/5));
    file_name_temp=sprintf('Atomic_number_%d_T_%d.csv',A,Te);
    file_name_temp_im=sprintf('Atomic_number_%d_T_%d.jpg',A,Te);
else
    title1='Input temperature choice';
    prompt1={'Do you wish to use  a range of temperatures (enter ''0'') or a single temperature (enter temperature in eV)'};
    def1={'0'};
    answer1=eval(char(inputdlg(prompt1,title1,1,def1,options)));
    if answer1==0
        title2='Input temperature range';
        prompt2={'Please enter the start Temperature (eV)',...
            'Please enter the step size (eV) ',...
            'Please enter the end Temperature (eV)'};
        def2={'0','0.1','15'};
        answer2=inputdlg(prompt2,title2,1,def2,options);
        Te_start=str2num(answer2{1});
        Te_step=str2num(answer2{2});
        Te_end=str2num(answer2{3});
        Te=Te_start:Te_step:Te_end;
        file_name_temp=sprintf('Atomic_number_%d_Tstart_%d_Tend_%d.csv',A,Te_start,Te_end);
        file_name_temp_im=sprintf('Atomic_number_%d_Tstart_%d_Tend_%d.png',A,Te_start,Te_end);
        if Te(1)==0;Te(1)=0.00001;end                                       %can't use 0, so set to small number
    else Te=answer1;end
end
potentialdata=load('Ionisation_5.csv');                           %Ionisation table input
IP=zeros(A,1);
for i=1:A
    IP(i)=potentialdata(i,A);
end
edata=load('electron_order2.dat');                                %Electron remaining 
outer_e_sn=zeros(A+1,1); i=0;
for fixing=A:-1:1
    i=i+1;
    outer_e_sn(i)=edata(fixing);
end
outer_e_sn(A+1)=0;
for i = 1:length(Te)                                                        %Temperature loops
    n(i,1)=1;
    S(i,1)=(((9E-6)*outer_e_sn(1)*((Te(i)/IP(1))^(1/2)))/((IP(1)^(3/2))*(4.88+(Te(i)/IP(1)))))*exp(-IP(1)/Te(i));
    alphaR(i,1)=(5.2E-14)*((IP(1)/Te(i))^(1/2))*0*(0.429+(0.5*log10(IP(1)/Te(i)))+(0.469*((Te(i)/IP(1))^(1/2))));
    alpha3b(i,1)=((2.97E-27)*outer_e_sn(1))/((Te(i)*(IP(1)^2))*(4.88+(Te(i)/IP(1))));
    for j=2:length(IP)
        S(i,j)=(((9E-6)*outer_e_sn(j)*((Te(i)/IP(j))^(1/2)))/((IP(j)^(3/2))*(4.88+(Te(i)/IP(j)))))*exp(-IP(j)/Te(i));
        alphaR(i,j)= (5.2E-14)*((IP(j)/Te(i))^(1/2))*(j-1)*(0.429+(0.5*log10(IP(j)/Te(i)))+(0.469*((Te(i)/IP(j))^(1/2))));
        alpha3b(i,j)=((2.97E-27)*outer_e_sn(j))/((Te(i)*(IP(j)^2))*(4.88+(Te(i)/IP(j))));
        n(i,j)=(S(i,j-1)/(alphaR(i,j)+ne*alpha3b(i,j)))*n(i,j-1);
    end
    Sbare=(((9E-6)*1*((Te(i)/IP(A))^(1/2)))/((IP(A)^(3/2))*(4.88+(Te(i)/IP(A)))))*exp(-IP(A)/Te(i));
    alphaRbare=(5.2E-14)*(0.469)*A;
    n(i,length(IP)+1)=(n(i,length(IP)).*Sbare)./(alphaRbare);
    nt(i)=sum(n(i,:));
    frac(i,:)=n(i,:)./nt(i);
end
im=figure(1);
plot(Te,frac,'LineWidth',2)
set(gca,'FontSize',26)
set(gcb,'FontSize',26)
set(gcf, 'color', 'w');
xlabel('Te (eV)','FontName','Times','FontSize',26)
ylabel('n_{ion stage Z}/n_{total}','FontName','Times','FontSize',26)
output_sn(:,1)=Te.';                                                        %Set up output files
output_sn(:,2:length(frac(1,:))+1)=frac;
axis([Te_start Te_end 0.01 1])
maximize
saveas(im, file_name_temp_im)
dlmwrite(file_name_temp, output_sn, ',');