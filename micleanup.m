%% Cleanup
% Multi-function extension of Froissart doublet cleanup from Chebfun.
function [r, pol, res, zer, z, f, w] = ...
    micleanup(r, pol, res, zer, z, f, w, Z, F, cleanup_tol) 
% Remove spurious pole-zero pairs.

% Find negligible residues:
ii = find(max(abs(res./max(abs(F),[],1)),[],2) < cleanup_tol );
ni = length(ii);
if ( ni == 0 )
    return
elseif ( ni == 1 )
    warning('CHEBFUN:aaa:Froissart','1 Froissart doublet');
else
    warning('CHEBFUN:aaa:Froissart',[int2str(ni) ' Froissart doublets']);
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = [];
    f(jj) = [];
end

% Remove support points z from sample set:
for jj = 1:length(z)
    F(Z == z(jj)) = [];
    Z(Z == z(jj)) = [];
end
m = length(z);
M = length(Z);

% Build Multi Function Loewner Matrix
SF = spdiags(F, 0, M, M);
C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix
L=[];
for i=1:k
    sm{i}=spdiags(F(:,i).',0,M,M);
    Li=sm{i}*C - C*diag(fz(i,:)); 
    L=[L;Li(J,:)];
end


% % Build Loewner matrix:
% % SF = spdiags(F, 0, M, M);
% % Sf = diag(f);
% % C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
% % A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(L, 0);
w = V(:,m);

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w);
[pol, res, zer] = prz(r, z, f, w);

end % End of CLEANUP().
