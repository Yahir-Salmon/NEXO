import sys, pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'src'/'integration'))
import numpy as np
import pytest
from nexo_numerics import *


def test_synthetic_reproducibility():
    a,b=stations(seed=42);c,d=stations(seed=42)
    assert np.array_equal(a,c) and np.array_equal(b,d)

def test_idw_exact_at_stations():
    x,y=stations(n=15); assert np.allclose(idw(x,y,x),y)

def test_reconstruction_stays_finite_when_station_removed():
    x,y=stations(n=21); q=np.array([[21,12],[30,33]])
    for idx in range(len(y)):
        assert np.isfinite(idw(np.delete(x,idx,0),np.delete(y,idx),q)).all()
        w,_,_=rbf_fit(np.delete(x,idx,0),np.delete(y,idx))
        assert np.isfinite(w).all()

def test_svd_tikhonov_stabilizes_small_singular_directions():
    x,y=stations(n=25,noise=2)
    _,A,s=rbf_fit(x,y,lam=0)
    w_reg,_,_=rbf_fit(x,y,lam=1)
    w_none,_,_=rbf_fit(x,y,lam=1e-11)
    assert np.linalg.norm(w_reg)<np.linalg.norm(w_none)
    assert s.max() > s.min()

def test_rank_known():
    A=np.array([[1,2,3],[2,4,6],[0,1,1]],float)
    assert np.linalg.matrix_rank(A)==2

def test_gradient_analytic():
    x=np.linspace(0,4,100);y=np.linspace(0,3,90);X,Y=np.meshgrid(x,y)
    gx,gy=gradient_2d(3*X+4*Y+7,x[1]-x[0],y[1]-y[0])
    assert np.max(abs(gx-3))<1e-10 and np.max(abs(gy-4))<1e-10

def test_mass_conservation_periodic():
    rng=np.random.default_rng(2);c=rng.uniform(1,40,(50,60))
    after=transport_step(c,2,-1,20,100,100,2)
    assert abs(np.sum(after)-np.sum(c))<1e-9

def test_solver_rejects_unstable_dt():
    with pytest.raises(ValueError,match='Unstable'):
        transport_step(np.ones((10,10)),3,1,30,10,10,20)

def test_mesh_convergence():
    err=[convergence_run(n) for n in [24,48,96]]
    assert err[2]<err[1]<err[0],err

def test_loocv_metrics_are_valid():
    xy,y=stations(n=12)
    metrics=loocv(xy,y)
    assert metrics['rmse']>=metrics['mae'] and len(metrics['predictions'])==12
