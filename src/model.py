"""Independent negative-binomial Bayesian regression by importance sampling.

The posterior is evaluated exactly. A multivariate Student-t proposal, located
at the MAP with curvature from finite differences, is corrected by likelihood /
proposal weights. Four independent batches permit Monte Carlo stability checks.
This is not a Gaussian/Laplace posterior approximation.
"""
import numpy as np
from scipy.special import gammaln, logsumexp, digamma, expit
from scipy.optimize import minimize
from scipy.stats import multivariate_t
import arviz as az

def design(x, center=None, scale=None):
    x=np.asarray(x,float)
    if x.ndim==1: x=x[:,None]
    if center is None: center=x.mean(axis=0)
    if scale is None: scale=x.std(axis=0,ddof=1)
    if np.any(scale<=0): raise ValueError('Constant predictor')
    return np.column_stack([np.ones(len(x)),(x-center)/scale]),center,scale

def log_posterior(theta,X,y):
    t=np.atleast_2d(theta); p=X.shape[1]
    # The last coordinate is log(phi); include its Jacobian.
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        eta=t[:,:p]@X.T; z=t[:,p,None];phi=np.exp(z)
        total=np.logaddexp(z,eta)
        ll=(gammaln(y+phi)-gammaln(phi)-gammaln(y+1)+phi*(z-total)+y*(eta-total)).sum(axis=1)
        lp=-.5*np.sum((t[:,:p]/2.5)**2,axis=1)-phi[:,0]+z[:,0]
        result=ll+lp
    return np.where(np.isfinite(result),result,-np.inf)

def log_gradient(theta,X,y):
    """Analytic score, avoiding cancellation in finite-difference optimization."""
    beta=theta[:-1];z=theta[-1];phi=np.exp(z)
    eta=X@beta;total=np.logaddexp(z,eta)
    score_eta=y-(y+phi)*expit(eta-z)
    grad_beta=X.T@score_eta-beta/2.5**2
    score_phi=digamma(y+phi)-digamma(phi)+z+1-total-(phi+y)*np.exp(-total)
    grad_z=phi*np.sum(score_phi)+1-phi
    return np.r_[grad_beta,grad_z]

def proposal(X,y):
    p=X.shape[1]; start=np.r_[np.log(np.mean(y)+.1),np.zeros(p-1),0.]
    f=lambda t:-log_posterior(t,X,y)[0]
    opt=minimize(f,start,jac=lambda t:-log_gradient(t,X,y),method='BFGS',options={'gtol':1e-5,'maxiter':1000})
    mode=opt.x; h=1e-3;d=len(mode);eye=np.eye(d)*h
    H=np.empty((d,d));f0=f(mode)
    for i in range(d):
        H[i,i]=(f(mode+eye[i])-2*f0+f(mode-eye[i]))/h**2
        for j in range(i):
            H[i,j]=H[j,i]=(f(mode+eye[i]+eye[j])-f(mode+eye[i]-eye[j])-f(mode-eye[i]+eye[j])+f(mode-eye[i]-eye[j]))/(4*h*h)
    if np.linalg.eigvalsh(H).min()<=0: raise RuntimeError('MAP curvature not positive definite')
    return mode,np.linalg.inv(H),bool(opt.success)

def fit(x,y,seed=0,n_per_batch=8192):
    y=np.asarray(y,float)
    if np.any(y<0) or not np.all(y==np.floor(y)): raise ValueError('NB observations must be nonnegative integers')
    X,center,scale=design(x); mode,cov,success=proposal(X,y)
    dist=multivariate_t(loc=mode,shape=1.2*cov,df=7)
    rng=np.random.default_rng(seed)
    for attempt in range(4):
        batches=[]; logs=[]; batch_means=[]
        n=n_per_batch*2**attempt
        for b in range(4):
            t=dist.rvs(size=n,random_state=rng)
            lw=log_posterior(t,X,y)-dist.logpdf(t)
            w=np.exp(lw-logsumexp(lw))
            batches.append(t);logs.append(lw);batch_means.append(w@t)
        theta=np.concatenate(batches);lw=np.concatenate(logs)
        weights=np.exp(lw-logsumexp(lw));ess=1/np.sum(weights**2)
        _,k=az.psislw(lw)  # diagnostic only; actual integration uses raw weights
        if ess>=4000 and k<.7: break
    means=weights@theta
    sd=np.sqrt(weights@((theta-means)**2))
    mcse=np.std(batch_means,axis=0,ddof=1)/2
    diagnostics=dict(importance_ess=float(ess),pareto_k=float(k),proposals=len(theta),
                     max_mcse_over_sd=float(np.max(mcse/sd)),map_optimizer_success=success,
                     numerical_pass=bool(ess>=4000 and k<.7 and np.max(mcse/sd)<.05))
    return dict(theta=theta,weights=weights,center=center,scale=scale,diagnostics=diagnostics,mode=mode,cov=cov)

def predict(fit,x,seed=0,size=30000):
    X,_,_=design(x,fit['center'],fit['scale'])
    rng=np.random.default_rng(seed)
    index=rng.choice(len(fit['theta']),size=size,p=fit['weights'])
    t=fit['theta'][index];phi=np.exp(t[:,-1])
    mu=np.exp(t[:,:-1]@X.T)
    return rng.negative_binomial(phi[:,None],phi[:,None]/(phi[:,None]+mu))

def parameter_summary(fit):
    t=fit['theta'].copy();w=fit['weights'];p=t.shape[1]-1
    if p>1:
        t[:,1:p]/=fit['scale']
        t[:,0]-=t[:,1:p]@fit['center']
    t[:,-1]=np.exp(t[:,-1])
    rows=[]
    for i in range(t.shape[1]):
        s=np.argsort(t[:,i]);cw=np.cumsum(w[s]);qs=np.interp([.025,.5,.975],cw,t[s,i])
        mean=w@t[:,i];sd=np.sqrt(w@((t[:,i]-mean)**2))
        rows.append(dict(parameter='intercept' if i==0 else ('phi' if i==p else f'beta{i}'),mean=mean,sd=sd,lower=qs[0],median=qs[1],upper=qs[2]))
    return rows
