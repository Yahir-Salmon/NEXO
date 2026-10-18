function A=rbf_design(xy, centers, width)
% Gaussian RBF design (dimensionless). XY, centers share km or m.
assert(width>0 && size(xy,2)==2 && size(centers,2)==2);
A=zeros(size(xy,1),size(centers,1));
for j=1:size(centers,1)
    dsq=sum((xy-centers(j,:)).^2,2);
    A(:,j)=exp(-dsq/(2*width^2));
end
end
