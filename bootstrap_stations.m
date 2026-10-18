function result=bootstrap_stations(xy,y,q,method,lambda,width,B,seed)
% Station-pairs bootstrap, SYNTHETIC DEMO ONLY until real spatial model calibrated.
% Bootstrap percentiles are NOT guaranteed prediction-interval coverage.
if nargin<4,method='RBF';end
if nargin<5,lambda=.05;end
if nargin<6,width=15;end
if nargin<7,B=200;end
if nargin<8,seed=20261018;end
assert(B>=10 && size(xy,1)>5);
rng(seed,'twister');
n=size(xy,1);samples=zeros(size(q,1),B);
for b=1:B
    ix=randi(n,n,1);
    samples(:,b)=reconstruct(xy(ix,:),y(ix),q,method,lambda,width);
end
result.mean=mean(samples,2);
result.sd=std(samples,0,2);
ordered=sort(samples,2);
lo=max(1,ceil((B+1)*.025));hi=min(B,ceil((B+1)*.975));
result.percentile_025=ordered(:,lo);
result.percentile_975=ordered(:,hi);
result.description='Station bootstrap variability, not calibrated prediction coverage';
end
