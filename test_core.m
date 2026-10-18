% Standard MATLAB no third-party toolbox. Print errors and assert failures.
here=fileparts(mfilename('fullpath'));addpath(genpath(here));
[xy,y]=synthetic_stations(16,123,0);
assert(norm(idw(xy,y,xy)-y)<1e-9);
for k=1:size(xy,1)
    tr=true(size(y));tr(k)=false;
    assert(all(isfinite(reconstruct(xy(tr,:),y(tr),xy(k,:),'RBF',.05,15))));
end
[x,ygrid]=meshgrid(0:100:1000,0:100:900);[gx,gy]=field_gradient(3*x+4*ygrid,100,100);
assert(max(abs(gx(:)-3))<1e-9);assert(max(abs(gy(:)-4))<1e-9);
C=rand(40,40)+3;N=advect_diffuse_periodic(C,2,-1,20,100,100,1);
assert(abs(sum(C(:))-sum(N(:)))<1e-9);
a=run_periodic_experiment(24,1.5,10,500);b=run_periodic_experiment(48,1.5,10,500);
assert(b.rmse<a.rmse);
fprintf('MATLAB core tests: passed (on current MATLAB runtime)\n');
