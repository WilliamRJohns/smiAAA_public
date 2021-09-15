function [pol]=przd(poles,residues)
%Use the g-eigenvalue problem from AAA along with the deflation idea from
%zpf to calculate the zeroes of a rational function
rv=residues;
pv=poles;
count=0;

while(abs(sum(rv))<10^-10)
%Preform a deflation
count=count+1;
rv(1) = [];             %Remove the first residue  
    first_pv = pv(1);   %Save the first pole
    pv(1) = [];         %Remove the first poles from the list
    
    %Recalculate the residues over the new set of poles
    fr_pv = zeros(length(pv),1);
    for i = 1:length(pv)
        fr_pv(i,1) = first_pv -pv(i);
    end
    rv = rv.*fr_pv;
    rv = rv/norm(rv);    
end

%Build and solve the general e-value problem
m=length(pv);
B = eye(m+1);
B(1,1) = 0;
E = [0 rv.'; ones(m, 1) diag(pv)];
pol = eig(E, B);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

%fprintf("%d deflations preformed\n",count)
       
end
