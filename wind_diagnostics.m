function [divergence,rotation]=wind_diagnostics(U,V,dx,dy)
% div v=dU/dx+dV/dy; curl_z=dV/dx-dU/dy.
[Ux,Uy]=gradient(U,dx,dy);[Vx,Vy]=gradient(V,dx,dy);
divergence=Ux+Vy;rotation=Vx-Uy;
end
