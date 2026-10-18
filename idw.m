function z = idw(xy, y, q, power)
% IDW. All coordinates must have consistent metric units; no toolboxes.
if nargin < 4, power=2; end
assert(size(xy,2)==2 && size(q,2)==2 && size(xy,1)==numel(y));
assert(power>0 && all(isfinite(xy(:))) && all(isfinite(y(:))));
z=zeros(size(q,1),1);y=y(:);
for i=1:size(q,1)
    d=sqrt(sum((xy-q(i,:)).^2,2));
    idx=find(d<1e-12,1);
    if ~isempty(idx), z(i)=y(idx);
    else
        w=max(d,1e-12).^(-power);
        z(i)=sum(w.*y)/sum(w);
    end
end
end
