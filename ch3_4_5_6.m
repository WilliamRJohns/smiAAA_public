% ex4a.m             
%
% Fitting 1st column of the admittance matrix of 6-terminal system 
% (power system distribution network)
%
% -Reading frequency admittance matrix Y(s) from disk.
% -Extracting 1st column: f(s) (contains 6 elements)
% -Fitting f(s) using vectfit3.m 
%   -Initial poles: 25 linearly spaced complex pairs (N=50)
%   -5 iterations
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%
close all;clear all;
%clear all;
clc;

disp('Reading data from file ...') %--> s(1,Ns), bigY(Nc,Nc,Ns)
fid1=fopen('fdne.txt','r');
Nc=fscanf(fid1,'%f',1);
Ns=fscanf(fid1,'%f',1);
bigY=zeros(Nc,Nc,Ns); s=zeros(1,Ns);
for k=1:Ns
  s(k)=fscanf(fid1,'%e',1);
  for row=1:Nc
    for col=1:Nc
      dum1=fscanf(fid1,'%e',1);
      dum2=fscanf(fid1,'%e',1);   
      bigY(row,col,k)=dum1+j*dum2;
    end
  end
end
s=i*s;
fclose(fid1);

%Extracting first column:
for n=1:Nc
  f(n,:)=squeeze(bigY(n,1,:)).';
end  
freq=s/(2*pi*1i);
%Introduction Figure
figure()
loglog(freq,abs(f),'-*','Linewidth',1)
title('Six entries of Admitance Matrix')
xlabel('Frequency Hz')
ylabel('Abs(Entries)')

%Approximate the 6 functions with miAAA and 10 Lawson iternations
[aaafl,wjl,aaaf,zj,wj,fj,err] = miaaa(f,s,1e-6,true,10);
wden = wj(1:length(zj)); wnum = wj(length(zj)+1:end);
[poles_aaa,res_aaa,pfaaaf,~,bestpoly]=properrational(zj.',wnum,wden,fj.',aaaf,s);

%Figure miAAA 
figure()
semilogx(freq,abs(f-aaafl)+10^-13,'Linewidth',1.5)
title('miAAA approximation errors tol=10^{-6}')
xlabel('Frequency Hz')
ylabel('Error')
legend('|f_1(Z)-B_1(Z)|','|f_2(Z)-B_2(Z)|','|f_3(Z)-B_3(Z)|','|f_4(Z)-B_4(Z)|','|f_5(Z)-B_5(Z)|','|f_6(Z)-B_6(Z)|')

%Figure MIMO Lawson
figure()
plot(freq,max(abs(f-aaafl),[],1)+10^-13,'b',freq,max(abs(f-aaaf),[],1)+10^-13,'r','Linewidth',1.5)
title('miAAA-L Errors tol=10^{-6}')
xlabel('Frequency Hz')
ylabel('Max Error')
legend('miAAA-L','miAAA')

%Approximate the six functions with smiAAA with 10 Lawson iternations
[saaafl,swj,saaaf,szj,swj,sfj,serr] = smiaaa(f,s,1e-6,true,10,1);
swden = swj; swnum = swj;
[spoles_aaa,sres_aaa,spfaaaf,~,sbestpoly]=properrational(szj.',swnum,swden,sfj.',saaaf,s);

%Figure smiAAA 
figure()
semilogy(freq,abs(f-saaafl)+10^-13,'Linewidth',1.5)
title('smiAAA approximation errors tol=10^{-6}')
xlabel('Frequency Hz')
ylabel('Error')
legend('|f_1(Z)-B_1(Z)|','|f_2(Z)-B_2(Z)|','|f_3(Z)-B_3(Z)|','|f_4(Z)-B_4(Z)|','|f_5(Z)-B_5(Z)|','|f_6(Z)-B_6(Z)|') 

%Figure stable MIMO Lawson
figure()
plot(freq,max(abs(f-saaafl),[],1)+10^-13,'b',freq,max(abs(f-saaaf),[],1)+10^-13,'r','Linewidth',1.5)
title('smiAAA-L Errors tol=10^{-6}')
xlabel('Frequency Hz')
ylabel('Max Error')
legend('smiAAA-L','smiAAA')

%Figure Partial Fractions smiAAA
figure();
loglog(freq,abs(saaaf-spfaaaf),'Linewidth',1.5)
legend('|B_1(Z)-pf_1(Z)|','|B_2(Z)-pf_2(Z)|','|B_3(Z)-pf_3(Z)|','|B_4(Z)-pf_4(Z)|','|B_5(Z)-pf_5(Z)|','|B_6(Z)-pf_6(Z)|')
xlabel('Frequency Hz')
ylabel('Error')
title('Error in Partial Fraction Conversion smiAAA')

%Figure Pole Comparison
figure();
scatter(real(poles_aaa),imag(poles_aaa),'rx','Linewidth',2);hold on;
scatter(real(spoles_aaa),imag(spoles_aaa),'bo','Linewidth',2)
xline(0);
legend('miAAA poles','smiAAApoles')

%Convergence Plot for Different Values of ref
tol=1e-8;
figure();
title('Convergence of smiAAA')
set(gca,'yscale','log')
hold on;
for ref=[.25 1 2 6 100]
[lpaaaf,pwj,paaaf,pzj,~,pfj,err] = smiaaa(f,s,tol,false,10,ref);
semilogy([11:length(err)],err(11:end),'-*');
%xline(length(err),'-',num2str(ref),'LabelHorizontalAlignment','center');
disp(length(err))
end
hold off;
legend('ref=.25','ref=1','ref=2','ref=6','ref=100');
ylabel('Log error')
xlabel('Iteration')

