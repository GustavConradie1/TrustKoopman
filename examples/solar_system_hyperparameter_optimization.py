import torch
import torch.nn.functional as F
from torch.autograd import Function
from scipy.special import kv
import scipy.io as sio
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import os

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

    all_error_bounds = []

    for k in range(d):
        coefs_svd = torch.linalg.lstsq(PX, x[k, :].unsqueeze(1)).solution

        coefs = Sig @ U.T @ coefs_svd

        coefs_svd2 = U @ S2 @ coefs
        
        x_step = DATA[k,0:maxx-1]
        delta = torch.sqrt(x_step @ x_step.T/M -
                             (2/M) * ((PX @ coefs_svd2).T @ x_step).real +
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
        all_error_bounds.append(torch.mean(error))
    all_error_bounds = torch.stack(all_error_bounds)
    print('error bounds')
    print(torch.mean(all_error_bounds))
    return torch.mean(all_error_bounds)

def compute_exact_errors(s, v):
    """Compute exact Koopman errors for all steps."""
    PX = ((s * r1) ** v) * bessel_kv(v, s*r1)  
    PY = ((s * r2) ** v) * bessel_kv(v, s*r2) 
    PX_unseen = ((s * r3) ** v) * bessel_kv(v,s*r3)

    G = PX.T @ PX / M
    G = (G+G.T)/2
    A = PX.T @ PY / M
    L = PY.T @ PY / M
    L = (L+L.T)/2
    #print(G)

    eigvals, U = torch.linalg.eigh(G)
    idx = torch.argsort(torch.abs(eigvals), descending=True)
    eigvals = eigvals[idx[:r]]
    U = U[:, idx[:r]]
    U = U / torch.linalg.norm(U, dim=0, keepdim=True)
    Sig = torch.sqrt(torch.diag(eigvals))
    S2 = torch.diag(1.0 / torch.diag(Sig))

    K2 = S2 @ (U.T @ A @ U) @ S2

    all_exact_errors = []

    G = PX_unseen.T @ PX_unseen / M
    #print(G)

    for k in range(d):
        coefs_svd = torch.linalg.lstsq(PX, x[k, :].unsqueeze(1)).solution
        coefs = Sig @ U.T @ coefs_svd
        coefs_svd2 = U @ S2 @ coefs
        fut_coefs = coefs
        x_step = DATA[k,maxx:maxx+maxx2]
        err = torch.sqrt(x_step @ x_step.T/M -
                             (2/M) * ((PX_unseen @ coefs_svd2).T @ x_step).real +
                             coefs_svd2.T @ (G @ coefs_svd2))
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
        all_exact_errors.append(torch.mean(exact_errors))
    all_exact_errors = torch.stack(all_exact_errors)
    print('exact errors')
    print(torch.mean(all_exact_errors))
    return (torch.mean(all_exact_errors))


def compute_pointwise_errors(s, v):
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

    all_pointwise_errors = []

    for k in range(d):
        coefs_svd = torch.linalg.lstsq(PX, x[k, :].unsqueeze(1)).solution
        coefs = Sig @ U.T @ coefs_svd
        coefs_svd2 = U @ S2 @ coefs
        fut_coefs = coefs
        exact_val=DATA[k,maxx-1]
        pred_val=PY[-1,:] @ coefs_svd2
        err = torch.abs(exact_val-pred_val)
        pointwise_errors = [err.squeeze()]
        for i in range(steps):
            exact_val=DATA[k,maxx+i]
            fut_coefs = K2 @ fut_coefs
            fut_coefs_svd = U @ S2 @ fut_coefs
            pred_val = PY[-1,:] @ fut_coefs_svd
            err = torch.abs(exact_val-pred_val)
            pointwise_errors.append(err.squeeze())
        pointwise_errors=torch.stack(pointwise_errors)
        all_pointwise_errors.append(torch.mean(pointwise_errors))
    all_pointwise_errors=torch.stack(all_pointwise_errors)
    print('pointwise errors')
    print(torch.mean(all_pointwise_errors))
    return (torch.mean(all_pointwise_errors))


def plot_loss_landscape(s_min=1/20000, s_max=1/1000, ns=20,
                        v_min=0.1, v_max=5.0, nv=20):
    s_vals = np.logspace(np.log10(s_min), np.log10(s_max), ns)
    v_vals = np.linspace(v_min, v_max, nv)

    S, V = np.meshgrid(s_vals, v_vals)

    bounds_err = np.zeros((nv, ns))
    exact_err  = np.zeros((nv, ns))

    total_iters = ns * nv
    pbar = tqdm(total=total_iters, desc="Evaluating Loss Landscape")

    for i in range(nv):
        for j in range(ns):

            s0 = torch.tensor(S[i, j], device=device)
            v0 = torch.tensor(V[i, j], device=device)
            print(s0)
            print(v0)

            with torch.no_grad():
                bounds_err[i, j] = compute_error(s0, v0).item()
                exact_err[i, j]  = compute_exact_errors(s0, v0).item()

            pbar.update(1)

    pbar.close()

    fig = plt.figure(figsize=(14, 6))

    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    surf1 = ax1.plot_surface(S, V, bounds_err,
                             cmap='viridis', alpha=0.85, edgecolor='none')
    ax1.set_title("Koopman Error Bound")
    ax1.set_xlabel("s parameter")
    ax1.set_ylabel("v parameter")
    ax1.set_zlabel("Bound")
    fig.colorbar(surf1, ax=ax1, shrink=0.6)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf2 = ax2.plot_surface(S, V, exact_err,
                             cmap='plasma', alpha=0.85, edgecolor='none')
    ax2.set_title("Exact Prediction Error")
    ax2.set_xlabel("s parameter")
    ax2.set_ylabel("v parameter")
    ax2.set_zlabel("Exact Error")
    fig.colorbar(surf2, ax=ax2, shrink=0.6)

    plt.tight_layout()
    plt.show()

    return S, V, bounds_err, exact_err

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

#loss landscape - for testing
#S, V, bound_surface, exact_surface = plot_loss_landscape(
#    s_min=1/100, s_max=1/5, ns=10,
#    v_min=0.5, v_max=1.5, nv=10
#)

s = torch.tensor(np.log(1/5), device=device, dtype=dtype, requires_grad=True)
v = torch.tensor(np.log(1), device=device, dtype=dtype, requires_grad=True)
optimizer = torch.optim.Adam([s, v], lr=1e-2)
torch.autograd.set_detect_anomaly(True)

epochs = 501
loss_history = []
best_loss = float('inf')
best_params = (None, None)
best_epoch = 0
exact_error_history = []
pointwise_error_history = []

for epoch in tqdm(range(epochs), desc="Training Progress"):
    optimizer.zero_grad()
    s0 = torch.exp(s)
    v0 = torch.exp(v)
    loss = compute_error(s0, v0)
    loss.backward()
    optimizer.step()
    loss_val = loss.item()
    loss_history.append(loss_val)

    if loss_val < best_loss:
        best_loss = loss_val
        best_params = (torch.exp(s).item(), torch.exp(v).item())
        best_epoch = epoch

    with torch.no_grad():
        exact_errors = compute_exact_errors(s0, v0)
        exact_error_history.append(exact_errors.cpu().numpy())
        pointwise_errors = compute_pointwise_errors(s0,v0)
        pointwise_error_history.append(pointwise_errors.cpu().numpy())

    print(f"\n Best at epoch {best_epoch}: loss={best_loss:.6e}, "
      f"s={best_params[0]:.6f}, v={best_params[1]:.6f}")

#save results
results_dir = os.path.join(data_dir, "results")
os.makedirs(results_dir, exist_ok=True)

loss_history_np = np.array(loss_history)
exact_error_history_np = np.array(exact_error_history)
pointwise_error_history_np = np.array(pointwise_error_history)

np.savez(os.path.join(results_dir, "optimization_results.npz"),
         loss_history=loss_history_np,
         exact_error_history=exact_error_history_np,
         pointwise_error_history=pointwise_error_history_np,
         cent=cent.detach().cpu().numpy(),
         best_epoch=best_epoch,
         best_loss=best_loss,
         best_s=best_params[0],
         best_v=best_params[1])

#plot results
plt.figure(figsize=(8, 5))
plt.plot(loss_history, 'b-o', label='Loss (L2 error bound)')
plt.plot(exact_error_history, 'g-s', label='Exact L2 error')
plt.plot(pointwise_error_history, 'r-s', label='Exact pointwise error')
plt.axvline(best_epoch, color='r', linestyle='--', alpha=0.7)
plt.scatter(best_epoch, best_loss, color='red', zorder=5)
plt.text(best_epoch + 0.5, best_loss,
         f'Best: epoch={best_epoch}\nloss={best_loss:.3e}\ns={best_params[0]:.3f}, v={best_params[1]:.3f}',
         color='red', fontsize=9, va='bottom')

plt.xlabel('Epoch')
plt.ylabel('Error')
plt.title('Loss vs Exact Error during Optimization')
plt.grid(True)
plt.yscale('log')
plt.legend()
plt.tight_layout()
final_fig_path = os.path.join(results_dir, "loss_vs_epoch.png")
plt.savefig(final_fig_path, dpi=300, bbox_inches='tight')
plt.show()


