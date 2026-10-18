function [xy,y]=synthetic_stations(n,seed,noise,cluster)
% SYNTHTIC laboratory domain in local km: not actual stations.
if nargin<1,n=24;end
if nargin<2,seed=20261018;end
if nargin<3,noise=1;end
if nargin<4,cluster=false;end
rng(seed,'twister');
if cluster
    xy=[28,28]+randn(n,2).*[11,10];xy=max([1,1],min([59,49],xy));
else
    xy=[1,1]+rand(n,2).*[58,48];
end
y=synthetic_truth(xy(:,1),xy(:,2))+noise*randn(n,1);
end
