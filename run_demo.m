% Execute FROM anywhere with: run('matlab/run_demo.m')
here=fileparts(mfilename('fullpath'));root=fileparts(here);
addpath(genpath(here));
[xy,y]=synthetic_stations(24,20261018,1,false);
[xx,yy]=meshgrid(linspace(1,59,96),linspace(1,49,80));
[z,info]=reconstruct(xy,y,[xx(:),yy(:)],'RBF',0.05,15);
z=reshape(z,size(xx));
cv=loocv_stations(xy,y,'RBF',0.05,15);
benchmark=run_periodic_experiment(48,1.5,10,500);
f=figure('Name','NEXO: reconstrucción sintética');contourf(xx,yy,z,15,'LineStyle','none');hold on;
scatter(xy(:,1),xy(:,2),36,y,'filled','MarkerEdgeColor','w');hold off;
axis equal tight;colorbar;title('SINTÉTICO: reconstrucción RBF regularizada');
xlabel('x local (km)');ylabel('y local (km)');
exportgraphics(f,fullfile(root,'results','figures','matlab_reconstruction.png'));
fprintf('Sintético LOOCV-RMSE: %.4f ug/m3\n',cv.rmse);
fprintf('Sintético PDE RMSE: %.4f ug/m3; balance masa %.3g\n',benchmark.rmse,benchmark.relative_mass_error);
