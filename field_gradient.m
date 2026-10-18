function [Cx,Cy,magnitude]=field_gradient(C,dx,dy)
% Spatial gradient; spacing MUST be meters for concentration gradients in ug/m4.
assert(dx>0 && dy>0 && ismatrix(C));
% MATLAB gradient: first output is derivative along columns (x).
[Cx,Cy]=gradient(C,dx,dy);
magnitude=hypot(Cx,Cy);
end
