import torch
import torch.nn.functional as F
from torch.autograd import Function
from scipy.special import kv
import scipy.io as sio
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
from scipy.io import savemat

#try delay embeddings to state space

class BesselKv(Function):
    @staticmethod
    def forward(ctx, v, x):
        v_np = float(v.detach().cpu().numpy())
        x_np = x.detach().cpu().numpy().astype(np.float64)
        y = torch.as_tensor(kv(v_np, x_np), dtype=torch.float64, device=x.device)
        ctx.save_for_backward(v, x, y)
        return y

    @staticmethod
    def backward(ctx, grad_output):
        v, x, y = ctx.saved_tensors
        v_np = float(v.detach().cpu().numpy())
        x_np = x.detach().cpu().numpy().astype(np.float64)

        dK_dx_np = np.array(-0.5 * (kv(v_np - 1, x_np) + kv(v_np + 1, x_np)), dtype=np.float64)

        eps = 1e-4
        kv_plus = kv(v_np + eps, x_np)
        kv_minus = kv(v_np - eps, x_np)
        dK_dv_np = (kv_plus - kv_minus) / (2 * eps)

        dK_dx = torch.tensor(dK_dx_np, dtype=x.dtype, device=x.device)
        dK_dv = torch.tensor(dK_dv_np, dtype=x.dtype, device=x.device)

        grad_x = grad_output * dK_dx
        grad_v = (grad_output * dK_dv).sum()
        return grad_v, grad_x

def bessel_kv(v, x):
    if not torch.is_tensor(v):
        v = torch.tensor(v, dtype=x.dtype, device=x.device)
    if not torch.is_tensor(x):
        x = torch.tensor(x, dtype=torch.float64, device=device)
    return BesselKv.apply(v, x)


def compute_error(s, v):
    PX = ((s * r1) ** v) * bessel_kv(v, (s*r1))
    PY = ((s * r2) ** v) * bessel_kv(v, (s*r2))

    G = (PX.T @ PX) / M
    G=(G+G.T)/2
    A = (PX.T @ PY) / M
    L = (PY.T @ PY) / M
    L=(L+L.T)/2

    eigvals, U = torch.linalg.eigh(G)
    idx = torch.argsort(torch.abs(eigvals), descending=True)
    eigvals = eigvals[idx[:r]]
    U = U[:, idx[:r]]
    U = U / torch.linalg.norm(U, dim=0, keepdim=True)
    Sig = torch.sqrt(torch.diag(eigvals))
    S2 = torch.diag(1.0 / torch.diag(Sig))

    A2 = S2 @ (U.mT @ A @ U) @ S2
    L2 = S2 @ (U.mT @ L @ U) @ S2
    L2 = 0.5 * (L2 + L2.mT) 
    K2 = A2

    all_error_bounds = torch.tensor([0.0], device=device).repeat(steps+1)

    for k in range(d):
        coefs_svd = torch.linalg.lstsq(PX, x[k, :].unsqueeze(1)).solution

        coefs = Sig @ U.T @ coefs_svd

        coefs_svd = U @ S2 @ coefs

        x_step = DATA[k,0:maxx-1]
        delta = torch.sqrt(x_step @ x_step.T/M -
                             (2/M) * ((PX @ coefs_svd).T @ x_step).real +
                             coefs.T @ coefs)

        Q = L2 - A2.mT @ K2
        Q = 0.5 * (Q + Q.mT)

        K_ops = [torch.eye(r, device=device, dtype=dtype)]
        for _ in range(steps):
            K_ops.append(K2 @ K_ops[-1])
        K_norm = torch.stack([torch.linalg.norm(K_ops[i],2) for i in range(steps + 1)])
        
        Kg = [coefs.squeeze()]
        KgQKg = [torch.sqrt(torch.real(Kg[0].T @ Q @ Kg[0]))]
        for i in range(1, steps):
            next_Kg = K2 @ Kg[-1]
            Kg.append(next_Kg)
            KgQKg.append(torch.sqrt(torch.real(next_Kg.T @ Q @ next_Kg)))
        KgQKg = torch.stack(KgQKg)

        error_list = []
        for i in range(steps):
            error_i = sum(K_norm[j] * KgQKg[i-j] for j in range(i+1)) + delta*K_norm[i+1]
            error_i=error_i.squeeze()
            error_list.append(error_i)
        error = torch.tensor([0.0], device=device).repeat(steps+1)
        error[0] = delta
        error[1:] = torch.stack(error_list)
        all_error_bounds=all_error_bounds+error
    all_error_bounds = all_error_bounds/d
    return all_error_bounds


def compute_exact_errors(s, v):
    PX = ((s * r1) ** v) * bessel_kv(v, s*r1) 
    PY = ((s * r2) ** v) * bessel_kv(v, s*r2)
    PX_unseen = ((s * r3) ** v) * bessel_kv(v,s*r3)

    G = PX.T @ PX / M
    G = (G+G.T)/2
    A = PX.T @ PY / M
    L = PY.T @ PY / M
    L = (L+L.T)/2

    eigvals, U = torch.linalg.eigh(G)
    idx = torch.argsort(torch.abs(eigvals), descending=True)
    eigvals = eigvals[idx[:r]]
    U = U[:, idx[:r]]
    U = U / torch.linalg.norm(U, dim=0, keepdim=True)
    Sig = torch.sqrt(torch.diag(eigvals))
    S2 = torch.diag(1.0 / torch.diag(Sig))

    K2 = S2 @ (U.T @ A @ U) @ S2

    all_exact_errors = torch.tensor([0.0], device=device).repeat(steps+1)

    G = PX_unseen.T @ PX_unseen / M

    for k in range(d):
        coefs_svd = torch.linalg.lstsq(PX, x[k, :].unsqueeze(1)).solution
        coefs = Sig @ U.T @ coefs_svd
        coefs_svd = U @ S2 @ coefs

        fut_coefs = coefs
        
        x_step = DATA[k,maxx:maxx+maxx2]
        err = torch.sqrt(x_step @ x_step.T/M -
                             (2/M) * ((PX_unseen @ coefs_svd).T @ x_step).real +
                             coefs_svd.T @ (G @ coefs_svd))
        exact_errors = [err.squeeze()]
        for i in range(steps):
            x_step = DATA[k,maxx+1+i:maxx+maxx2+1+i]
            fut_coefs = K2 @ fut_coefs
            fut_coefs_svd = U @ S2 @ fut_coefs
            err = torch.sqrt(x_step @ x_step.T/M -
                             (2/M) * ((PX_unseen @ fut_coefs_svd).T @ x_step).real +
                             fut_coefs_svd.T @ (G @ fut_coefs_svd))
            exact_errors.append(err.squeeze())
        exact_errors=torch.stack(exact_errors)
        all_exact_errors=all_exact_errors+exact_errors
    all_exact_errors = all_exact_errors/d
    return all_exact_errors

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dtype = torch.float64
torch.set_default_dtype(dtype)

data_dir = r"C:\Users" #add desired directory

# load data from .mat files
DATAp = sio.loadmat(os.path.join(data_dir, 'states_pluto2.mat'))['states_pluto2']
DATAc = sio.loadmat(os.path.join(data_dir, 'states_charon2.mat'))['states_charon2']
DATAh = sio.loadmat(os.path.join(data_dir, 'states_hydra2.mat'))['states_hydra2']
DATAs = sio.loadmat(os.path.join(data_dir, 'states_styx2.mat'))['states_styx2']
DATAn = sio.loadmat(os.path.join(data_dir, 'states_nix2.mat'))['states_nix2']
DATAk = sio.loadmat(os.path.join(data_dir, 'states_kerberos2.mat'))['states_kerberos2']

DATA = torch.tensor(
    np.hstack([DATAp, DATAc, DATAh, DATAs, DATAn, DATAk]),
    dtype=torch.float64
).to(device)

DATA=DATA.T

steps, r = 20, 50 
maxx = 10*365
maxx2 = 10*365

means = (DATA[:, :maxx]).mean(dim=1, keepdim=True) 
stds = (DATA[:, :maxx]).std(dim=1, keepdim=True)
DATA = (DATA - means) / stds

x = DATA[:, :maxx-1].to(device)
y = DATA[:, 1:maxx].to(device)
M = x.shape[1]
d = x.shape[0]

N = 200
torch.manual_seed(0)
cent = torch.rand(36, N, device=device)
x_min = x.min(dim=1).values.unsqueeze(1)
x_max = x.max(dim=1).values.unsqueeze(1)
cent = cent * (x_max - x_min) + x_min
cent = cent.T

x2 = DATA[:,maxx:maxx+maxx2]

r1 = torch.cdist(x.T, cent, p=2)
r2 = torch.cdist(y.T, cent, p=2)
r3 = torch.cdist(x2.T, cent, p=2)

error_bounds_init=compute_error(1/5,1)
print(error_bounds_init)
exact_errors_init=compute_exact_errors(1/5,1)
print(exact_errors_init)
error_bounds_best=compute_error(0.0522,4.73)
print(error_bounds_best)
exact_errors_best=compute_exact_errors(0.0522,4.73)
print(exact_errors_best)

data = {
    'error_bounds_init': error_bounds_init.detach().cpu().numpy(),
    'exact_errors_init': exact_errors_init.detach().cpu().numpy(),
    'error_bounds_best': error_bounds_best.detach().cpu().numpy(),
    'exact_errors_best': exact_errors_best.detach().cpu().numpy()
}

results_dir = os.path.join(data_dir, "results")
os.makedirs(results_dir, exist_ok=True)

savemat(os.path.join(results_dir, "hyperparameter_fixed_results.mat"), data)









