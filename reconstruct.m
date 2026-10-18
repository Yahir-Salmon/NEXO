function [z,diagnostics]=reconstruct(xy,y,q,method,lambda,width)
if nargin<4, method='RBF';end
if nargin<5, lambda=0.05;end
if nargin<6, width=15;end
switch upper(method)
    case 'IDW'
        z=idw(xy,y,q,2);
        diagnostics=struct('method','IDW','lambda',NaN,'condition_estimate',NaN);
    case 'RBF'
        A=rbf_design(xy,xy,width);
        [weights,diagnostics]=tikhonov_svd(A,y,lambda);
        z=rbf_design(q,xy,width)*weights;
        diagnostics.method='RBF_Tikhonov';diagnostics.lambda=lambda;
    otherwise, error('NEXO:InvalidMethod','Allowed: IDW, RBF');
end
end
