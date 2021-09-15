function [bestbcr,bestw,bcr,z,wj,fz,err] = smiaaa(f,Z,tol,normalize,iter,ref)
%Accepts a matrix of function values over Z with each row being a single
%function's values
%The values of Z should be on the imaginary axis for stability enforcement
%to make sense

%Computes Multifunction smiAAA rational approximations over a common set of poles
%in barycentric form such that each functions approximation 
%is within tol in inf-norm sense, and all poles are stable. 

%normalize=true, normalizes the functions before computing the approximation and renormalizes after
%iter = number of Lawson optimization iterations with fixed poles to maintain stabiliety
%ref=a,a is real, reflects unstable poles p to stable poles q with real(q)=-a(real(p))

%Returns the Lawson optimized barycentric approximation bestbcr and its
%associated wieghts bestw. The pre optimization approximation bcr and its
%weights w and the support points z and function values at the support
%poitns fz.


mmax=100;        %Max number of Support Points
k=size(f,1);     %The number of functions

if(normalize)    %Nomalize Each function
    norms=vecnorm(f,2,2);
    f=f(1:k,:)./norms(1:k);
end

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
    %Compute the deviations at each points in Z
    dev=(abs(f-bcr));
    err=[err max(dev,[],'all')];
    if(max(dev)<tol)
        break
    end
    
    %Find the next support points
    totaldev=sum(dev,1);    
    [~,I]=max(totaldev);    %Find the Point with the largest sum of deviations
    z=[z Z(I)];             %Update the Support Points
    fz=[fz f(:,I)];         %Update the function values at Support Points
    J(J==I)=[];             %Update the Sample point indices 
    Jz=[Jz I];              %Update the Support Point idices(used later in Lawson)
    
    %Build Loewner Matrices for each function
    C = [C 1./(Z.' - Z(I))]; %Next row of Cauchy Matrix
    L=[];
    for i=1:k
        Li=sm{i}*C - C*diag(fz(i,:)); 
        L=[L;Li(J,:)];
    end
    
    [~, ~, W] = svd(L, 0);      %Solve for the common weights
    w=W(:,n);
    
    if(n>1)
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
wj=w; % Store unoptimized weights
%Lawson Optimization (iterative wighted least squares)
%With poles (den weights) fixed
if(iter>=0)
%Initialize Variables of Lawson Iteration
m=length(z);
gama=2;             %Initial Lawson Exponent
lw=ones(1,M);       %Initial Lawson Weigths
lw=lw/norm(lw);
bestbcr=bcr;        %Best Approximation so far
bestw=repmat(w,2,1);        %w associated with best approximation;
lbcr=[];            %Approximation at current Lawson iteration
maxerror=max(abs(f-bestbcr),[],'all');    %Max err of best approx

L=[];           %Build Non-interpolatory Loewner Matrix
%replacement for the interpolation condition r(ti)~f(ti)
Lsupp=[diag(ones(m,1)) -1*diag(ones(m,1))];
for i=1:k
    Li=[sm{i}*C -1*C*diag(fz(i,:))];
    Li(Jz,:)=Lsupp;
    L=[L;Li];   
end

for l=1:iter %Lawson Iterations
lws=repmat(lw,1,k);
d=spdiags(sqrt(lws).',0,size(L,1),size(L,1));
%Get the Lawson weighted approximation weight vector
%[~, ~, W] = svd(d*L,0);
wnum=(d*L(:,m+1:end))\(-1*d*L(:,1:m)*w(1:m));    %Resolve num weights
w=[w(1:m) ;wnum];
D = C(J,:)*w(1:m);                              %Calculate the Denominator
    for i=1:k
        Nk = C(J,:)*(w(m+1:end).*fz(i,:).');        %Calculate the Numerator for each funtion
        lbcr(i,J)=Nk./D;                            %Set approximation values at non support points
        lbcr(i,Jz)=w(m+1:end)'.*f(i,Jz)./w(1:m)';   %Set approximation values at support points
    end
%Compare with the old approximation
lmaxerror=max(abs(f-lbcr),[],'all');
poles=przd(z.',w(1:m));
np_poles=poles(real(poles)>0);
if(maxerror<=lmaxerror)
    gama=gama;
    disp('Changed gama')
else
    fprintf('Optimized from %d to %d\n',maxerror,lmaxerror);
    maxerror=lmaxerror;
    bestbcr=lbcr;
    bestw=w;
end

%Update the Lawson wieghts(I extended to multiple functions with max)
testlw=lw.*((max(abs(f-lbcr),[],1)));  
testlw(find(testlw==Inf))=max(testlw(testlw~=Inf));
testlw(isnan(testlw))=mean(testlw(~isnan(testlw)));
if(~isempty(find(testlw==NaN)) | ~isempty(find(testlw==inf)))    %If all values were NaN or Inf they can't be fixed
    disp('Lawson terminated,Weights could not be fixed')
    break;
end
testlw(find(testlw==0))=mean(testlw(testlw~=0));   %avoid any 0's from perfect interpolation
testlw=testlw/norm(testlw);
if(norm(testlw-lw,inf)<1e-8) %This Tolerance should be tested
    disp('Lawson Terminated on Weight fixed point')
    break;
end
lw=testlw;
 
end
end

%Remove any support points with zero weight
io=find(w(1:m)==0);
io2=find(w(M+1:end)==0);
io=intersect(io,io2);
z(io)=[];
fz(:,io)=[];
w(io)=[];
if(~isempty(io))
    disp('zero wieght removed')
end

if(normalize)    %Re-nomalize Each function
    bcr=bcr(1:k,:).*norms(1:k);
    bestbcr=bestbcr(1:k,:).*norms(1:k);
    fz(1:k,:)=fz(1:k,:).*norms(1:k);
end

end