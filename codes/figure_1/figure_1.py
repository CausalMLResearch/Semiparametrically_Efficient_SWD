import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
    "text.usetex": False,
})

# read simulation results
df = pd.read_csv('../../results/figure_1/simulation_results_20250828_045609.csv')
global_min = min(df['mse_aipw1'].min(), df['mse_or'].min())
global_max = max(df['mse_aipw1'].max(), df['mse_or'].max())

def reshape_data(df_filtered, pi_1_unique, pi_2_unique, column):
    """reshape data into grid"""
    Z = np.full((len(pi_1_unique), len(pi_2_unique)), np.nan)
    for i, p1 in enumerate(pi_1_unique):
        for j, p2 in enumerate(pi_2_unique):
            mask = (df_filtered['pi_1'] == p1) & (df_filtered['pi_2'] == p2)
            if mask.any():
                Z[i, j] = df_filtered.loc[mask, column].values[0]
    return Z

def setup_plot(ax, pi_1_range, pi_2_range, z_min, z_max, label, title):
    """plot settings"""
    ax.text2D(0.5, 1.02, f'{label} {title}', transform=ax.transAxes,
              fontsize=20, va='bottom', ha='center', fontfamily='serif')

    ax.set_xlabel(r'$\pi_2$', fontsize=20, labelpad=20, fontfamily='serif')
    ax.set_ylabel(r'$\pi_1$', fontsize=20, labelpad=20, fontfamily='serif')
    ax.set_zlabel('MSE', fontsize=20, labelpad=20, rotation=90, fontfamily='serif')
    ax.tick_params(axis='x', labelsize=18, pad=8)
    ax.tick_params(axis='y', labelsize=18, pad=8)
    ax.tick_params(axis='z', labelsize=18, pad=8)
    ax.view_init(elev=25, azim=-45)
    ax.set_xlim(pi_2_range[0], pi_2_range[1])
    ax.set_ylim(pi_1_range[0], pi_1_range[1])
    ax.set_zlim(z_min, z_max)
    ax.set_xticks(np.linspace(pi_2_range[0], pi_2_range[1], 5))
    ax.set_yticks(np.linspace(pi_1_range[0], pi_1_range[1], 5))

    z_ticks = np.linspace(max(0, z_min), 0.24, 4)
    z_ticks[0] = max(0, z_min)
    z_ticks[-1] = 0.24
    ax.set_zticks(z_ticks)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1f}'))
    ax.zaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2f}'))

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontfamily('serif')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontfamily('serif')
    for tick in ax.zaxis.get_major_ticks():
        tick.label1.set_fontfamily('serif')

pi_1_range, pi_2_range = [0, 1], [0, 1]
df_filtered = df[(df['pi_1'].between(*pi_1_range)) & (df['pi_2'].between(*pi_2_range))]
pi_1_unique = np.sort(df_filtered['pi_1'].unique())
pi_2_unique = np.sort(df_filtered['pi_2'].unique())
PI_2, PI_1 = np.meshgrid(pi_2_unique, pi_1_unique)
Z_or = reshape_data(df_filtered, pi_1_unique, pi_2_unique, 'mse_or')
Z_aipw = reshape_data(df_filtered, pi_1_unique, pi_2_unique, 'mse_aipw1')

fig = plt.figure(figsize=(12, 4))
norm = Normalize(vmin=global_min, vmax=global_max)
sm = ScalarMappable(norm=norm, cmap='turbo')

# or estimator
ax1 = fig.add_subplot(121, projection='3d')
surf1 = ax1.plot_surface(PI_2, PI_1, Z_or,
                         cmap='turbo',
                         vmin=global_min,
                         vmax=global_max,
                         rstride=1,
                         cstride=1,
                         linewidth=0,
                         antialiased=True)
contour1 = ax1.contourf(PI_2, PI_1, Z_or, 20, zdir='z', offset=-0.1, cmap='turbo', vmin=global_min, vmax=global_max)
setup_plot(ax1, pi_1_range, pi_2_range, -0.1, global_max * 1.1, '(a)', 'MSE of an OR estimator')

# aipw estimator
ax2 = fig.add_subplot(122, projection='3d')
surf2 = ax2.plot_surface(PI_2, PI_1, Z_aipw,
                         cmap='turbo',
                         vmin=global_min,
                         vmax=global_max,
                         rstride=1,
                         cstride=1,
                         linewidth=0,
                         antialiased=True)
contour2 = ax2.contourf(PI_2, PI_1, Z_aipw, 20, zdir='z', offset=-0.1, cmap='turbo', vmin=global_min, vmax=global_max)
setup_plot(ax2, pi_1_range, pi_2_range, -0.1, global_max * 1.1, '(b)', 'MSE of an AIPW estimator')

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)

n_colorbar_ticks = 4
cbar_ticks = np.linspace(global_min, global_max, n_colorbar_ticks)
cbar.set_ticks(cbar_ticks)
cbar.ax.set_yticklabels([f'{v:.2f}' for v in cbar_ticks])

for label in cbar.ax.get_yticklabels():
    label.set_fontfamily('serif')

plt.subplots_adjust(left=0.02, right=0.90, top=0.88, bottom=0.05, wspace=0.05)
filename = '2_period_mse.pdf'
plt.savefig(filename, bbox_inches='tight', dpi=1000, pad_inches=0.2)
print(f'saved: {filename}')
plt.show()
plt.close()
