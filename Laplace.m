stp = 0.01; % boundary step discretization
c = 0.5+0.5*i; % power sum centre

z = [0: stp :2]; z = [z 2+i *[0: stp :1]]; % L- shaped boundary definition
z = [z [2: - stp :1]+i ]; z = [ z 1+ i *[1: stp :2]];
z = [z [1: - stp :0]+2* i ]; z = [z i *[2: - stp :0]].';
p = polygon ([0 2 2+ i 1+ i 1+2* i 2* i ]) ;

u = real(z).^2; % Dirichlet boundary condition

N = 10+ ceil ( log ( length ( z))); % smooth part , nr. of terms
[a , H] = polyfitA (z -c ,u , N); % Vandermonde with Arnoldi
wz = @(x ) polyvalA (a ,H ,x (:) -c) ; % function handle , smooth part
kz = imag(wz(c));

up = u - real ( wz (z)) ; % singular part of b.c.
[r , pol,~,~,ZJ,~,~,~] = aaa ( up ,z ,'mmax' ,1000 , 'lawson' ,0) ; % AAA approximation
m = find ( isinpoly ( pol ,p ,1E-16) ==0) ; % find poles outside
polm = pol (m);


%lets do some smiaaa science 
[bestbcr,bestw,bcr,zj,wj,fz,err] = bmiaaa(up.',z.',1e-7,false,0,1,p);


B = 1./( z - polm.') ; % coefficient matrix
B (: , end +1) = 1;
B = [ real(B) -imag(B) ];
b = reshape (B\ up ,[] ,2) ; % LS approximation
b = b (:,1) +i*b (: ,2) ;
wp = @(x) [1./(x(:) - polm.') ones(size(x(:)))]*b ; % function handle , singular part
kp = imag(wp(c));

%aaa error with ext poles
disp('aaa error with all poles');
disp(norm(up-r(z),'inf'));

disp('aaa error with ext poles');
disp(norm(up-real(wp(z)),'inf'));

w = @(x) wp(x)+ wz(x(:)) -i *( kz + kp ); % function handle , solution
maxerr = norm (u - real(w(z)) ,'inf') % maximum error at the boundary