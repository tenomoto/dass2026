import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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

def plot_rmse(hist, thist, ghist=None):
    h = eo.calc_h(hist['q0'], hist['q2'], hist['q4'])
    h_t = eo.calc_h(thist['q0'], thist['q2'], thist['q4'])
    rmse_q0 = np.sqrt(np.mean((hist['q0'] - thist['q0'])**2, axis=1))
    rmse_q2 = np.sqrt(np.mean((hist['q2'] - thist['q2'])**2, axis=1))
    rmse_q4 = np.sqrt(np.mean((hist['q4'] - thist['q4'])**2, axis=1))
    rmse_h = np.sqrt(np.mean((h - h_t)**2, axis=1))
    rmse = dict(
        q0 = rmse_q0, q2 = rmse_q2, q4 = rmse_q4, h = rmse_h
    )
    fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
    colors = ['k', 'r', 'b', 'g']
    for c, name in zip(colors, ['q0', 'q2', 'q4', 'h']):
        ax.plot(rmse[name], c=c, ls='-', lw=2.0, label=name+' DA')
    if ghist is not None:
        h_g = eo.calc_h(ghist['q0'], ghist['q2'], ghist['q4'])
        grmse_q0 = np.sqrt(np.mean((ghist['q0'] - thist['q0'])**2, axis=1))
        grmse_q2 = np.sqrt(np.mean((ghist['q2'] - thist['q2'])**2, axis=1))
        grmse_q4 = np.sqrt(np.mean((ghist['q4'] - thist['q4'])**2, axis=1))
        grmse_h = np.sqrt(np.mean((h_g - h_t)**2, axis=1))
        rmse['gq0'] = grmse_q0
        rmse['gq2'] = grmse_q2
        rmse['gq4'] = grmse_q4
        rmse['gh'] = grmse_h
        for c, name in zip(colors, ['q0', 'q2', 'q4', 'h']):
            ax.plot(rmse['g'+name], c=c, ls='--', lw=2.0, label=name+' guess')
    if 'q0s' in hist.keys():
        sprd_q0 = np.mean(hist['q0s'], axis=1)
        sprd_q2 = np.mean(hist['q2s'], axis=1)
        sprd_q4 = np.mean(hist['q4s'], axis=1)
        rmse['q0s'] = sprd_q0
        rmse['q2s'] = sprd_q2
        rmse['q4s'] = sprd_q4
        for c, name in zip(colors, ['q0', 'q2', 'q4']):
            ax.plot(rmse[name+'s'], c=c, ls=':', lw=2.0, label=name+' sprd')
        ax.set_ylabel('RMSE, Spread', fontsize=12)
    else:
        ax.set_ylabel('RMSE', fontsize=12)
    ax.set_xlabel('cycle')
    ax.legend(ncol=3)
    plt.show()
    plt.close()
    return rmse

def plot_state(hist, plot_f=False):
    fig = plt.figure(figsize=(8, 6))
    if plot_f:
        ncol = 4
    else:
        ncol = 3
    icol = 1
    axs = []
    for name in ['q0', 'q2', 'q4']:
        ax = fig.add_subplot(1, ncol, icol)
        axs.append(ax)
        if icol > 1:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('Time')
        var = hist[name]
        t = np.arange(var.shape[0])
        x = np.arange(var.shape[1])
        p = ax.pcolormesh(x, t, var, shading='auto', vmin=-1.5, vmax=1.5, cmap='RdBu_r')
        ax.set_title(name)
        ax.set_xlabel('Location')
        icol += 1
    fig.colorbar(p, ax=axs, shrink=0.6, pad=0.01)
    if plot_f:
        ax = fig.add_subplot(1, ncol, icol)
        ax.set_yticklabels([])
        var = hist['f']
        t = np.arange(var.shape[0])
        x = np.arange(var.shape[1])
        p = ax.pcolormesh(x, t, var, shading='auto', cmap='PuOr', vmin=-0.15, vmax=0.15)
        fig.colorbar(p, ax=ax, shrink=0.6, pad=0.01)
        ax.set_title('f')
        ax.set_xlabel('Location')
    if 'q0s' in hist.keys():
        fig_s = plt.figure(figsize=(8, 6))
        icol = 1
        axs = []
        for name in ['q0s', 'q2s', 'q4s']:
            ax = fig_s.add_subplot(1, ncol, icol)
            axs.append(ax)
            if icol > 1:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Time')
            var = hist[name]
            t = np.arange(var.shape[0])
            x = np.arange(var.shape[1])
            p = ax.pcolormesh(x, t, var, shading='auto', vmin=0.0, vmax=0.4, cmap='viridis')
            ax.set_title(name)
            ax.set_xlabel('Location')
            icol += 1
        fig_s.colorbar(p, ax=axs, shrink=0.6, pad=0.01)
        if plot_f:
            ax = fig_s.add_subplot(1, ncol, icol)
            ax.set_yticklabels([])
            var = hist['fs']
            t = np.arange(var.shape[0])
            x = np.arange(var.shape[1])
            p = ax.pcolormesh(x, t, var, shading='auto', cmap='viridis', vmin=0.0, vmax=0.05)
            fig_s.colorbar(p, ax=ax, shrink=0.6, pad=0.01)
            ax.set_title('f')
            ax.set_xlabel('Location')
        return fig, fig_s
    elif 'q0_pf' in hist.keys():
        fig_s = plt.figure(figsize=(8, 6))
        icol = 1
        axs = []
        for name in ['q0_pf', 'q2_pf', 'q4_pf']:
            ax = fig_s.add_subplot(1, ncol, icol)
            axs.append(ax)
            if icol > 1:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Time')
            var = hist[name]
            t = np.arange(var.shape[0])
            x = np.arange(var.shape[1])
            p = ax.pcolormesh(x, t, var, shading='auto', cmap='coolwarm', 
            norm=mcolors.SymLogNorm(linthresh=0.03, linscale=0.03,
            vmin=-1.0, vmax=1.0, base=2))
            ax.set_title(name)
            ax.set_xlabel('Location')
            icol += 1
        fig_s.colorbar(p, ax=axs, shrink=0.8, pad=0.01)
        return fig, fig_s
    else:
        return fig