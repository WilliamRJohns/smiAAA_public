function [bestbcr,bestw,bcr,z,wj,fz,err] = symmetricsmiaaah2o(f,Z,tol,normalize,iter,ref)
%Accepts a matrix of function values over Z with each row being a single
%function's values
%Computes Multifunction AAA rational approximations over a common set of poles
%in barycentric and proper rational form such that each functions approximation 
%is with tol in inf-norm sense, where iter = number of Lawson optimization iterations
%and normalize=true, normalizes the functions before computing the approximation and renormalizes after

%Returns the Lawson optimized barycentric approximation bestbcr and its
%associated wieghts bestw. The pre optimization approximation bcr and its
%weights w and the support points z and function values at the support
%poitns fz.


mmax=111;        %Max number of Support Points
k=size(f,1);     %The number of functions

if(normalize)    %Nomalize Each function     for H2
    norms=norm(f,'fro')^2;
    f=f/norms;
end

% %domain normalization (might help cleanup)
% %Causes Numerical problems with reflection of poles
% a=2*1i/(max(Z)-min(Z));
% b=1i-a*max(Z);
% Z=a*Z+b;

%Initilize Variables
z=[];           %Support Points
fz=[];          %functions Values at support points
M=length(Z);    
bcr=mean(f,2);  %Initial "barycentric" "rational" 'approximation' for selection of first support point
J=1:M;          %Index vector of Samples not chosen as support points
Jz=[];          %Index vector of the Support Points
C=[];           %Cauchy Matrix
sm={};          %Spares Matrices for Loewner Matrix
err=[];         %Vector for error at each iteration
for i=1:k
    sm{i}=spdiags(f(i,:).',0,M,M);
end

%Primary Iteration to compute the barycentric rational approximations
for n=1:mmax+1
    %Compute the H_2 error
    [~,~,h2]=comp_error(f,bcr);
    if(h2<tol)
        disp('H_2 Under Tolerance')
        break
    end
    
    dev=(abs(f-bcr));
    err=[err max(dev,[],'all')];
    %Find the next support point and its conjugate for Hermitian Symetry
    totaldev=sum(dev,1);    
    [~,I]=max(totaldev);         %Find the Point with the largest sum of deviations
    z=[z Z(I) conj(Z(I))];       %Update the Support Points
    fz=[fz f(:,I) conj(f(:,I))]; %Update the function values at Support Points
    J(J==I)=[];                  %Update the Sample point indices 
    Jz=[Jz I];                   %Update the Support Point idices(used later in Lawson)
    
    %Build Loewner Matrices for each function
    C = [C 1./(Z.'-Z(I)) 1./(Z.'-conj(Z(I)))]; %Next row of Cauchy Matrix
    L=[];
    for i=1:k
        Li=sm{i}*C - C*diag(fz(i,:)); 
        L=[L;Li(J,:)];
    end
    
    Lh=[];
    %Compute the complex-real operator from fast AAA
    for ii=1:size(L,2)
        Lh= [Lh [real(L(:,ii)) -1*imag(L(:,ii)) ;imag(L(:,ii)) real(L(:,ii))]];
    end
    H=kron(eye(n),[1 0;0 1;1 0;0 -1]);
    
    [~, ~, W] = svd(Lh*H, 0);      %Solve for the common weights
    wht=W(:,end);
    wh=H*wht;
    w=wh(1:2:end)+1i*wh(2:2:end);
    if(true)
    %Compute the poles
    poles=przd(z.',w);
    %Find the non passive poles and replace them with passive poles
    np_poles=poles(real(poles)>0);
    while (~isempty(np_poles))

         %fprintf("%d Non-passive pole found at iteration %d\n",length(np_poles),n);
        
         p_poles=-1*ref*real(np_poles)+1i*imag(np_poles);
         for ii=1:length(np_poles)
             w = w.*(z.'-p_poles(ii))./(z.'-np_poles(ii));
         end
         poles=przd(z.',w);
         np_poles=poles(real(poles)>0);
    
    end
    end
    bcr = f;                    %Sets approximation values = functions values at support points
    D = C*w;                    %Calculate the Denominator
    for i=1:k
        Nk = C*(w.*fz(i,:).');  %Calculate the Numerator
        bcr(i,J)=Nk(J)./D(J);   %Sets approximation values at non support points
    end  
end

cleanup_flag=true;
if (cleanup_flag)
    [poles_aaa,res_aaa,pfaaaf,~,bestpoly]=properrational(z.',w,w,fz.',f,Z);
    [test]=micleanup(1,poles_aaa,res_aaa,1,z,fz,1,Z,f.',1e-13);
end


wj=w; % Store unoptimized weights
wh_den=wh;
%Lawson Optimization (iterative wighted least squares)
%With poles (den weights) fixed
if(iter>0)
%Initialize Variables of Lawson Iteration

m=length(z);
bestbcr=bcr;        %Best Approximation so far
bestw=repmat(w,2,1);        %w associated with best approximation;
lbcr=[];            %Approximation at current Lawson iteration
[~,~,maxerror]=comp_error(f,bestbcr);

L=[];           %Build Non-interpolatory Loewner Matrix
%replacement for the interpolation condition r(ti)~f(ti)
%Lsupp=[diag(ones(m/2,1)) diag(ones(m/2,1)) diag(ones(m/2,1)) diag(ones(m/2,1))];
Lsupp=zeros(m/2,2*m);
Lsupp(:,1:2:m)=diag(ones(m/2,1));
Lsupp(:,m+1:2:end)=diag(ones(m/2,1));
for i=1:k
    Li=[sm{i}*C C*diag(fz(i,:))];
    Li(Jz,:)=Lsupp;
    L=[L;Li];   
end

%Compute the complex-real operator from fast AAA
Lh=[];  %non interpolatory loewner matrix
for ii=1:size(L,2)
    Lh= [Lh [real(L(:,ii)) -1*imag(L(:,ii)) ;imag(L(:,ii)) real(L(:,ii))]];
end
A=Lh(:,2*m+1:end);
B=Lh(:,1:2*m);
H=kron(eye(m/2),[1 0;0 1;1 0;0 -1]);
for l=1:iter %Lawson Iterations
disp("iter= "+l);

%Fix den weights
b=B*wh_den;
wnum=A*H\b;
wh_num=H*wnum;
%optimized num weights
wnum=wh_num(1:2:end)+1i*wh_num(2:2:end);

%Fix num weights
a=A*wh_num;
wden=B*H\a;
wh_den=H*wden;
%optimized den weights
wden=wh_den(1:2:end)+1i*wh_den(2:2:end);

%reflect any unstable poles added in optimization
poles=przd(z.',wden);
%Find the unstable
np_poles=poles(real(poles)>0);
while (~isempty(np_poles))
         fprintf("%d unstable pole found at optimization iteration %d\n",length(np_poles),l);
         p_poles=-1*ref*real(np_poles)+1i*imag(np_poles);
         for ii=1:length(np_poles)
             wden = wden.*(z.'-p_poles(ii))./(z.'-np_poles(ii));
         end
         poles=przd(z.',w);
         np_poles=poles(real(poles)>0);    
end
w=[wden,wnum];
D = C(J,:)*w(1:m);                                      %Calculate the Denominator
    for i=1:k
        Nk = C(J,:)*(w(m+1:end).*fz(i,:).');            %Calculate the Numerator for each funtion
        lbcr(i,J)=Nk./D;                                %Set approximation values at non support points
        lbcr(i,Jz)=w(m+1:2:end)'.*f(i,Jz)./w(1:2:m)';   %Set approximation values at support points
    end
%Compare with the old approximation
[~,~,lmaxerror]=comp_error(f,lbcr);

if(maxerror>lmaxerror)
    fprintf('Optimized H_2 from %d to %d\n',maxerror,lmaxerror);
    maxerror=lmaxerror;
    bestbcr=lbcr;
    bestw=w;
end

end
else
    bestbcr=bcr;
    bestw=w;
end

cleanup_flag=false;
if (cleanup_flag)
    wden=wj;
    [poles_aaa,res_aaa,pfaaaf,~,bestpoly]=properrational(z.',wnum,wden,fz.',f,Z);
    [test]=micleanup(1,poles_aaa,res_aaa,1,z,fz,1,Z,f.',1e-13);
end


%Remove any support points with zero weight
m=length(w)/2;
io=find(w(1:m)==0);
io2=find(w(m+1:end)==0);
io=intersect(io,io2);
z(io)=[];
fz(:,io)=[];
w(io)=[];
if(~isempty(io))
    disp('zero wieght removed')
end


if(normalize)    %Re-nomalize Each function
    bcr=bcr*norms;
    bestbcr=bestbcr*norms;
    fz(1:k,:)=fz*norms;
end

end