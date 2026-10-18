function tab=convergence_study(resolutions,u,D,T)
% Runs independent synthetic PDE benchmarks.
if nargin<1,resolutions=[24 48 96];end
if nargin<2,u=1.5;end
if nargin<3,D=10;end
if nargin<4,T=500;end
rms=zeros(numel(resolutions),1);balance=rms;
for k=1:numel(resolutions)
    out=run_periodic_experiment(resolutions(k),u,D,T);
    rms(k)=out.rmse;balance(k)=out.relative_mass_error;
end
tab=table(resolutions(:),rms,balance,'VariableNames',...
    {'cells_per_axis','reference_rmse_ug_m3','relative_mass_error'});
end
