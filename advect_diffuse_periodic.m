function [next,cfl]=advect_diffuse_periodic(C,u,v,D,dx,dy,dt,source,k)
% Conservative finite-volume periodic transport, 2D, SI units.
% dC/dt + div(v*C) = D Laplacian C + source - k*C.
% First-order upwind flux and central diffusion. Stability constraint checked.
if nargin<9, source=0;end
if nargin<10,k=0;end
assert(dx>0 && dy>0 && dt>0 && D>=0 && k>=0);
cfl=abs(u)*dt/dx+abs(v)*dt/dy+2*D*dt*(1/dx^2+1/dy^2)+k*dt;
if cfl>1+1e-12
    error('NEXO:Unstable','Explicit CFL/stability bound violated: %g >1',cfl);
end
right=C(:,[2:end,1]);down=C([2:end,1],:);
if u>=0,advX=u*C;else,advX=u*right;end
if v>=0,advY=v*C;else,advY=v*down;end
Fx=advX-D*(right-C)/dx;
Fy=advY-D*(down-C)/dy;
next=C-dt*((Fx-Fx(:,[end,1:end-1]))/dx+...
    (Fy-Fy([end,1:end-1],:))/dy)+dt*(source-k*C);
end
