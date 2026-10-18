"""Executed reference tests + generated figures; no output is a real-air-quality result."""
import json, pathlib, sys, platform
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'/'integration'))
from nexo_numerics import *
from scipy.spatial import ConvexHull
from matplotlib.path import Path


def run():
    rng_seed=20261018
    xy,y=stations(24,seed=rng_seed,noise=1)
    gridx=np.linspace(1,59,116);gridy=np.linspace(1,49,96)
    X,Y=np.meshgrid(gridx,gridy); q=np.column_stack((X.ravel(),Y.ravel()))
    true=synthetic_field(X,Y)
    w,A,s=rbf_fit(xy,y,lam=.05);zrbf=rbf_predict(q,xy,w).reshape(X.shape)
    zidw=idw(xy,y,q).reshape(X.shape)
    cv_idw=loocv(xy,y,'idw');cv_rbf=loocv(xy,y,'rbf',lam=.05)
    # Convex hull mask -> report interior / exterior separately
    inside=Path(xy[ConvexHull(xy).vertices]).contains_points(q).reshape(X.shape)
    summaries={}
    for method,z in [('IDW',zidw),('RBF_Tikhonov',zrbf)]:
        summaries[method]={'grid_hull_rmse_synthetic':float(np.sqrt(np.mean((z[inside]-true[inside])**2))),
                           'grid_extrapolation_rmse_synthetic':float(np.sqrt(np.mean((z[~inside]-true[~inside])**2))),
                           'loo_station_rmse_synthetic':cv_idw['rmse'] if method=='IDW' else cv_rbf['rmse']}
    fig,axs=plt.subplots(1,3,figsize=(15,4.6),constrained_layout=True)
    for ax,z,t in zip(axs,[true,zidw,zrbf],['Verdad conocida','Distancia inversa (IDW)','RBF y Tikhonov']):
        im=ax.pcolormesh(X,Y,np.where(inside,z,np.nan) if t!='Verdad conocida' else z,shading='auto',cmap='viridis',vmin=15,vmax=44)
        ax.scatter(xy[:,0],xy[:,1],s=22,c='white',edgecolors='#15202b',linewidths=.65,zorder=3)
        ax.set(xlabel='x local (km)',ylabel='y local (km)',title=t,aspect='equal')
    fig.colorbar(im,ax=axs,label='Concentración SINTÉTICA (µg/m³)',shrink=.8)
    fig.suptitle('NEXO · Experimento sintético · no son estaciones reales',size=14)
    fig.savefig(ROOT/'results/figures/compare_synthetic.png',dpi=165);plt.close(fig)
    U,s,Vt=np.linalg.svd(A,full_matrices=False)
    fig,ax=plt.subplots(figsize=(7,4));ax.semilogy(np.arange(1,len(s)+1),s,'o-',label='singulares σᵢ');
    ax.axhline(np.sqrt(.05),ls='--',label='escala de filtro √λ',c='tab:red')
    ax.set(xlabel='Índice i',ylabel='σᵢ (adimensional)',title='Espectro SVD · diseño RBF · SINTÉTICO');ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(ROOT/'results/figures/svd_synthetic.png',dpi=175);plt.close(fig)
    lambdas=[0.0001,0.001,0.01,0.05,0.1,0.3,1,3,10]
    err=[]
    for lam in lambdas:err.append(loocv(xy,y,'rbf',lam=lam)['rmse'])
    fig,ax=plt.subplots(figsize=(7,4));ax.semilogx(lambdas,err,'o-',label='RMSE dejar-una-estación-fuera');ax.axhline(cv_idw['rmse'],ls='--',label='IDW (LOOCV)');ax.set(xlabel='λ (adimensional; penalización de pesos)',ylabel='RMSE sintético (µg/m³)',title='Selección regularización fuera de muestra · SINTÉTICO');ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(ROOT/'results/figures/cv_synthetic.png',dpi=175);plt.close(fig)
    # mass and mesh convergence
    n=48;L=10000.;xx=np.arange(n)*L/n; initial=analytic_periodic(xx,xx,0,1.5,10.,L)
    next_step=transport_step(initial,1.5,0.,10.,L/n,L/n,10)
    conserved=float(abs(next_step.sum()-initial.sum())/initial.sum())
    orders=[24,48,96]; conv=[convergence_run(z) for z in orders]
    fig,ax=plt.subplots(figsize=(7,4));ax.loglog(orders,conv,'o-',label='Error L2 contra solución analítica');ax.set(xlabel='Número de celdas por eje',ylabel='Error RMS (µg/m³)',title='Convergencia de malla · transporte-difusión sintético');ax.grid(alpha=.25,which='both');ax.legend();fig.tight_layout();fig.savefig(ROOT/'results/figures/convergence_synthetic.png',dpi=175);plt.close(fig)
    summary={'all_results_are':'SYNTHETIC_NUMERICAL_EXPERIMENTS_NOT_REAL_MONTERREY_OBSERVATIONS',
       'seed':rng_seed,'stations_synthetic':len(xy),'measurement_noise_sd_ug_m3':1,
       'reconstruction':summaries,'rbf_matrix_rank_tol':int(np.linalg.matrix_rank(A)),
       'rbf_condition_number':float(np.linalg.cond(A)),'singular_values':[float(v) for v in s],
       'lambda_LOOCV':[{'lambda':k,'rmse':v} for k,v in zip(lambdas,err)],
       'mass_relative_error_single_periodic_step':conserved,
       'mesh_convergence':[{'n':n,'rms':e} for n,e in zip(orders,conv)],
       'executed_with':{'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__},
       'notes':['Known smooth synthetic concentration; all x and y are local km, NOT georeferenced station positions.',
                'Convex hull used to mask interpolation maps; grid outside hull reports extrapolation error separately.',
                'LOOCV uses stationary synthetic sensor coordinates only, distinct from independent grid truth.',
                'Periodic PDE solver represents a laboratory benchmark, not calibrated real dispersion.']}
    (ROOT/'results/metrics/synthetic_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    print(json.dumps({k:v for k,v in summary.items() if k in ['reconstruction','mass_relative_error_single_periodic_step','mesh_convergence','rbf_condition_number']},indent=2))
    return summary

if __name__=='__main__': run()
