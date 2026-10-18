function [beta,info]=tikhonov_svd(A,y,lambda,L)
% Solve min ||A*beta-y||_2^2 + lambda||L*beta||_2^2.
% Use augmented SVD for general L, no explicit inverse or normal eq.
if nargin<4 || isempty(L),L=eye(size(A,2));end
assert(lambda>=0 && all(isfinite(A(:))) && all(isfinite(y(:))));
B=[A; sqrt(lambda)*L]; b=[y(:); zeros(size(L,1),1)];
[U,S,V]=svd(B,'econ');s=diag(S);
if isempty(s), beta=zeros(size(A,2),1);info.rank=0;return;end
tol=max(size(B))*eps(max(s));active=s>tol;
proj=U(:,active)'*b;
beta=V(:,active)*(proj./s(active));
info.rank=sum(active);info.singular_values=svd(A,'econ');
info.residual_norm=norm(A*beta-y(:));
info.penalty_norm=norm(L*beta);
info.condition_estimate=cond(A);
end
