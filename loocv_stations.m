function report=loocv_stations(xy,y,method,lambda,width)
% Leave-one-station-out CV; no shared full-data fit (avoids leakage).
if nargin<4, lambda=.05;end
if nargin<5, width=15;end
n=numel(y);p=zeros(n,1);
for k=1:n
    train=true(n,1);train(k)=false;
    [estimate,~]=reconstruct(xy(train,:),y(train),xy(k,:),method,lambda,width);
    p(k)=estimate;
end
report.predictions=p; report.errors=p-y(:);
report.rmse=sqrt(mean(report.errors.^2));report.mae=mean(abs(report.errors));
end
