function [p_poles,err] = FastAAAcompare(f,Z,tol,normalize,iter)
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


mmax=100;         %Max number of Support Points
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
spferr=[];         %Vector for the stable pf error at each iteration
for i=1:k
    sm{i}=spdiags(f(i,:).',0,M,M);
end

%Primary Iteration to compute the barycentric rational approximations
for n=1:300
    if(length(z)==length(Z));
        fprintf('FastAAA ran out of support points!\n');
        p_poles=[]; %Flag if stable partial fraction did no converge
        break;
    end
    
    %Compute the deviations at each points in Z
    dev=(abs(f-bcr));
    err=[err max(dev,[],'all')];
    
    %Find the next support point and its conjugate for Hermitian Symetry
    totaldev=max(dev,[],1);      %I think they are chossing support points like this
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
    
    poles=przd(z.',w);
    %Find the unstable and stable poles
    np_poles=poles(real(poles)>0);
    p_poles=poles(real(poles)<=0);
    %if(~isempty(np_poles))
    %     fprintf("%d Non-passive pole found at iteration %d\n",length(np_poles),n);
    %end
    
    %compute the stable partial fraction
    if(~isempty(p_poles))
        %Compute the residues
        Cnum=(poles-z).^(-1);
        Cden=-1*Cnum.^2;
        res=((Cnum*(w.*fz.'))./(Cden*w)).';
        
        %keep only stable poles/residue pairs
        %p_poles=poles(real(poles)<=0);
        p_res=res(:,real(poles)<=0);
        spf=[];pf=[];
        for ii=1:k
            spf=[spf ; pfeval(Z, p_poles, p_res(ii,:).', [0])];
            pf=[pf ; pfeval(Z, poles, res(ii,:).', [0])];
        end
    rem=f-spf;
    rem2=f-pf;
    else
        rem=f; spf=zeros(size(f)); rem2=f;pf=zeros(size(f));
    end
    
    %spf=[];
    for ii=1:k
        %Compute the polynomial part
        polypart=polyfit(Z,rem(ii,:),1);
        polypart2=polyfit(Z,rem2(ii,:),1);
        %Save in bcr because convinient 
        spf(ii,:)=spf(ii,:)+polyval(polypart,Z);
        pf(ii,:)=pf(ii,:)+polyval(polypart2,Z);
    end
  
    dev2=(abs(f-pf));
    dev=(abs(f-spf));
    spferr=[spferr max(dev,[],'all')];
    %fprintf("Partial Fraction error %d\n",max(dev2,[],'all'));
    %fprintf("stable pf error %d\n",max(dev,[],'all'));
    
    if(max(dev,[],'all')<tol)
        disp('Stable pf Under Tolerance')
        fprintf('with %d stable poles\n',length(p_poles));     
        break
    end
    
    bcr = f;                    %Sets approximation values = functions values at support points
    D = C*w;                    %Calculate the Denominator
    for i=1:k
        Nk = C*(w.*fz(i,:).');  %Calculate the Numerator
        bcr(i,J)=Nk(J)./D(J);   %Sets approximation values at non support points
    end 
    err=[err max((abs(f-bcr)),[],'all')];
    %fprintf("Barycentric error %d\n",err(end));
end


if(normalize)    %Re-nomalize Each function
    bcr=bcr(1:k,:).*norms(1:k);
    bestbcr=bestbcr(1:k,:).*norms(1:k);
    fz(1:k,:)=fz(1:k,:).*norms(1:k);
end
%Find the best stable partial fraction
[besterr,idx]=min(spferr);
fprintf('Best FastAAA error %d, with %d support points\n',besterr,2*idx);



end