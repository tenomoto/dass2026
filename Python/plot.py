import numpy as np
import matplotlib.pyplot as plt
import eqocean as eo

def print_varnames(ax, x, dx=0.5, ffact=10):
    varnames = ["q0", "q2", "q4", "h", f"{ffact} x f"]
    nc = len(varnames)
    ax.set_xlim(0, nc * (2 * np.pi + dx) - dx)
    ax.set_ylim(0, 1)
    x_shift = x.copy()
    for v in varnames:
        ax.text(np.mean(x_shift), 0.5, v, ha='center', va='center', fontsize=11)
        x_shift = x_shift + 2 * np.pi + dx
        
    ax.axis('off')

def plot_waves_row(ax, x, q0, q2, q4, f, n, dx=0.5, ffact=10, add=False):
    h = eo.calc_h(q0, q2, q4)
    ymat = [q0, q2, q4, h, ffact * f]
    nc = len(ymat)
    if not add:
        ax.set_xlim(0, nc * (2 * np.pi + dx) - dx)
        ymin, ymax = np.min(ymat), np.max(ymat)
        
        if ymin == ymax:
            ax.set_ylim(ymin - 1, ymax + 1)
        else:
            ax.set_ylim(ymin, ymax)
            
    fmt = 'r--' if add else 'k-'
    x_shift = x.copy()
    for j in range(nc):
        if (n > 0) or (j < nc - 1):
            ax.plot(x_shift, ymat[j], fmt, linewidth=1)
            
            if not add:
                ax.plot([x_shift[0], x_shift[-1]], [0, 0], 'k-', linewidth=0.5)
            
        x_shift = x_shift + 2 * np.pi + dx
    if not add:
        ax.text(-0.3, 0, str(n), ha='right', va='center', fontsize=9)
        ax.axis('off')

def plot_waves(x, hist, thist=None):
    tmax = hist['q0'].shape[0]
    fig, axes = plt.subplots(nrows=tmax + 2, ncols=1, figsize=(8, 6))
    print_varnames(axes[0], x)
    for idx, n in enumerate(range(tmax - 1, -1, -1)):
        ax = axes[idx + 1]
        plot_waves_row(
            ax, x, 
            hist['q0'][n, :], hist['q2'][n, :], hist['q4'][n, :], hist['f'][n, :], 
            n, add=False
        )
        if thist is not None:
            plot_waves_row(
                ax, x, 
                thist['q0'][n, :], thist['q2'][n, :], thist['q4'][n, :], thist['f'][n, :], 
                n, add=True
            )
    axes[-1].axis('off')
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

