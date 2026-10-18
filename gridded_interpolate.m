function z=gridded_interpolate(xAxis,yAxis,Z,xyQuery,method)
% Interpolation for true RECTILINEAR grids only (not irregular sensor locations).
% Standard MATLAB interp2; no special toolbox needed.
% Grid rows indexed by yAxis and columns by xAxis.
if nargin<5,method='linear';end
assert(isvector(xAxis)&&isvector(yAxis)&&isnumeric(Z));
assert(size(Z,1)==numel(yAxis)&&size(Z,2)==numel(xAxis));
assert(all(diff(xAxis)>0)&&all(diff(yAxis)>0));
assert(size(xyQuery,2)==2);
assert(ismember(method,{'linear','spline'}));
z=interp2(xAxis(:)',yAxis(:),Z,xyQuery(:,1),xyQuery(:,2),method,NaN);
end
