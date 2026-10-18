function results=run_periodic_experiment(n,u,D,T)
% Closed-form Fourier-mode benchmark; all variables synthetic.
if nargin<1,n=64;end
if nargin<2,u=1.5;end
if nargin<3,D=10;end
if nargin<4,T=500;end
L=10000;h=L/n;gridx=(0:n-1)*h;
[X,Y]=meshgrid(gridx,gridx);
C=20+5*sin(2*pi*X/L).*sin(2*pi*Y/L);
initialMass=sum(C(:))*h*h;
maxDt=.6/(abs(u)/h +4*D/h^2);
nSteps=max(1,ceil(T/maxDt));dt=T/nSteps;
for it=1:nSteps
    C=advect_diffuse_periodic(C,u,0,D,h,h,dt);
end
Cexact=20+5*exp(-8*pi^2*D*T/L^2)*sin(2*pi*(X-u*T)/L).*sin(2*pi*Y/L);
results=struct('n',n,'u',u,'D',D,'T',T,'C',C,...
    'Cexact',Cexact,'L',L,'rmse',sqrt(mean((C(:)-Cexact(:)).^2)),...
    'relative_mass_error',abs(sum(C(:))*h*h-initialMass)/initialMass,...
    'dt',dt,'steps',nSteps);
end
