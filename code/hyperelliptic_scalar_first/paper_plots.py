"""Page-sized paper layouts, reusable from saved plot data without solving again."""
import json
from pathlib import Path
import matplotlib.pyplot as plt


def save_paper_plots(example_dir, plot_times, finest_states, analytic_times,
                     analytic_values, beta, rk4_steps):
    folder = Path(example_dir)
    data = dict(times=[float(t) for t in plot_times],
                states=[[float(v.real) for v in s] for s in finest_states],
                # Compute powers before discarding imaginary components.
                powers=[[float((s[j]*s[j+len(beta)]).real) for j in range(len(beta))]
                        for s in finest_states],
                analytic_times=[float(t) for t in analytic_times],
                analytic_states=[[float(v.real) for v in item['state']] for item in analytic_values],
                analytic_powers=[[float((item['state'][j]*item['state'][j+len(beta)]).real)
                                  for j in range(len(beta))] for item in analytic_values],
                analytic_rho=[float(item['inversion'].real) for item in analytic_values],
                beta=[int(b) for b in beta], steps=int(rk4_steps[-1]))
    (folder/'plot_data.json').write_text(json.dumps(data)+'\n')
    return render_paper_plots(folder, data)


def render_paper_plots(folder, data=None):
    folder = Path(folder)
    if data is None:
        data = json.loads((folder/'plot_data.json').read_text())
    n = len(data['beta'])
    def panel(ax, line, markers, name):
        ax.plot(data['times'], line, color='black', linestyle='--', linewidth=.8,
                label=f"RK4 ({data['steps']:,} steps)")
        ax.plot(data['analytic_times'], markers, color='black', marker='o',
                markersize=3, linestyle='None', label=f"sigma ({len(markers)} samples)")
        ax.set_xlabel(r'$z-z_0$', fontsize=15)
        ax.tick_params(labelsize=14)
        ax.legend([], [], title=name, title_fontsize=17, loc='best', framealpha=.9)
        ax.grid(alpha=.15)
    rows = 4 if n == 5 else 5
    fig, axes = plt.subplots(rows, 3, figsize=(11, 14.9), constrained_layout=True)
    # Complete u group followed by complete v group, with each panel named.
    for k, ax in enumerate(axes.flat):
        if k >= 2*n:
            ax.set_axis_off()
            continue
        tag, j = ('u', k) if k < n else ('v', k-n)
        panel(ax, [s[k] for s in data['states']],
              [s[k] for s in data['analytic_states']], rf'$\mathrm{{Re}}\,{tag}_{{{j+1}}}$')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles, labels, loc='center', fontsize=14, frameon=False)
    fig.savefig(folder/'modes_combined.png', dpi=200)
    figures = [fig]
    if n == 7:
        fig = plt.figure(figsize=(11, 14.9), constrained_layout=True)
        grid = fig.add_gridspec(4, 3)
        ax = fig.add_subplot(grid[0, :])
        panel(ax, [-sum(s)/sum(data['beta']) for s in data['powers']],
              data['analytic_rho'], r'$\rho(z)$')
        for j in range(9):
            ax = fig.add_subplot(grid[1+j//3, j%3])
            if j >= n:
                ax.set_axis_off()
                if j == 8:
                    ax.legend(handles, labels, loc='center', fontsize=14, frameon=False)
                continue
            panel(ax, [s[j] for s in data['powers']],
                  [s[j] for s in data['analytic_powers']], rf'$I_{{{j+1}}}$')
        fig.savefig(folder/'scalar_powers_combined.png', dpi=200)
        figures.append(fig)
    return figures
