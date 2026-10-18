"""Independent numerical reference implementation for reproducibility tests.
The research reconstruction and PDE code live in matlab/. Synthetic only here.
Coordinates in km (reconstruction) / m (transport). No real station claims.
"""
import numpy as np
from scipy.linalg import svdvals
from scipy.spatial.distance import cdist


def synthetic_field(x, y):
    """Smooth known function, micrograms/m3; ALL SYNTHETIC."""
    x, y = np.broadcast_arrays(np.asarray(x), np.asarray(y))
    return (18 + 20*np.exp(-((x-22)**2/(2*9**2)+(y-28)**2/(2*11**2)))
            + 8*np.exp(-((x-49)**2/(2*6**2)+(y-13)**2/(2*8**2)))
            + 3*np.sin(x/12)*np.cos(y/10))


def stations(n=24, seed=20261018, cluster=False, noise=1.0):
    rng = np.random.default_rng(seed)
    if cluster:
        xy = np.clip(rng.normal([28,28], [11,10], (n,2)), [1,1], [59,49])
    else:
        xy = rng.uniform([1,1], [59,49], (n,2))
    values = synthetic_field(xy[:,0],xy[:,1]) + rng.normal(0, noise, n)
    return xy, values


def idw(train_xy, train_y, query_xy, power=2):
    d = cdist(np.asarray(query_xy, float), np.asarray(train_xy,float))
    out = np.zeros(len(d))
    exact = d < 1e-12
    sel = exact.any(axis=1)
    out[sel] = np.asarray(train_y)[exact[sel].argmax(axis=1)]
    d = np.maximum(d[~sel], 1e-12)
    w = d**(-power)
    out[~sel] = (w @ np.asarray(train_y))/(w.sum(axis=1))
    return out


def gaussian_design(xy, centers, scale=15):
    return np.exp(-0.5*(cdist(np.asarray(xy), np.asarray(centers))/scale)**2)


def rbf_fit(xy, y, lam=0.05, scale=15):
    A = gaussian_design(xy,xy,scale)
    U,s,Vt = np.linalg.svd(A, full_matrices=False)
    # Direct SVD Tikhonov filter: minimizes ||Aw-y||^2 + lam||w||^2
    w = Vt.T @ ((s/(s*s+lam)) * (U.T @ np.asarray(y)))
    return w, A, s


def rbf_predict(query, centers, weights, scale=15):
    return gaussian_design(query, centers, scale) @ weights


def loocv(xy, values, kind='idw', lam=0.05, scale=15):
    p = np.zeros(len(values))
    for i in range(len(values)):
        take = np.arange(len(values)) != i
        if kind=='idw':
            p[i] = idw(xy[take],values[take],xy[i:i+1])[0]
        elif kind=='rbf':
            w,_,_ = rbf_fit(xy[take],values[take],lam,scale)
            p[i] = rbf_predict(xy[i:i+1],xy[take],w,scale)[0]
        else: raise ValueError(kind)
    err = values - p
    return {'rmse':float(np.sqrt(np.mean(err*err))),
            'mae':float(np.mean(np.abs(err))), 'predictions':p}


def gradient_2d(z, dx, dy):
    gy, gx = np.gradient(z, dy, dx, edge_order=2)
    return gx,gy


def transport_step(c, u, v, D, dx, dy, dt, source=0, decay=0):
    """Conservative upwind flux + central diffusion, doubly periodic; C in ug/m3.
    No-clipping. Check explicit timestep stability before advance.
    """
    if dx<=0 or dy<=0 or dt<=0 or D<0 or decay<0: raise ValueError('invalid physics')
    c=np.asarray(c,float)
    cfl = abs(u)*dt/dx + abs(v)*dt/dy + 2*D*dt*(dx**-2+dy**-2)+decay*dt
    if cfl>1+1e-12: raise ValueError(f'Unstable explicit scheme: CFL={cfl:.3g}>1')
    fx = u*(c if u>=0 else np.roll(c,-1,axis=1)) - D*(np.roll(c,-1,axis=1)-c)/dx
    fy = v*(c if v>=0 else np.roll(c,-1,axis=0)) - D*(np.roll(c,-1,axis=0)-c)/dy
    return c - dt*(fx-np.roll(fx,1,axis=1))/dx - dt*(fy-np.roll(fy,1,axis=0))/dy + dt*(source-decay*c)


def analytic_periodic(x,y,t,u,D,L,base=20,amplitude=5):
    X,Y=np.meshgrid(x,y)
    return base+amplitude*np.exp(-8*np.pi**2*D*t/L**2)*np.sin(2*np.pi*(X-u*t)/L)*np.sin(2*np.pi*Y/L)


def convergence_run(n, L=10000., t_final=500., u=1.5, D=10.):
    h = L/n
    x=np.arange(n)*h; y=x
    c=analytic_periodic(x,y,0,u,D,L)
    dt_max=.60/(abs(u)/h+4*D/h**2)
    steps=int(np.ceil(t_final/dt_max)); dt=t_final/steps
    for _ in range(steps): c=transport_step(c,u,0,D,h,h,dt)
    expected=analytic_periodic(x,y,t_final,u,D,L)
    return float(np.sqrt(np.mean((c-expected)**2)))
