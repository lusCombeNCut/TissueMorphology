#!/usr/bin/env python3
"""Plot circularity timeseries with error bars for selected stiffnesses."""

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'mathtext.fontset': 'dejavuserif',
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.size': 11,
    'axes.labelsize': 13,
    'legend.fontsize': 9,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
})

SELECTED_STIFFNESSES = [300.0, 1300.0]
COLOURS = {300.0: '#1f77b4', 1300.0: '#d62728'}

csv_path = '2d_sim_04042026/plots_circularity_timeseries/vertex2d_circularity_timeseries.csv'
out_path = '2d_sim_04042026/plots_circularity_timeseries/vertex2d_circularity_errorbar.png'

data = {}
with open(csv_path) as f:
    for row in csv.DictReader(f):
        s = float(row['stiffness_Pa'])
        if s not in SELECTED_STIFFNESSES:
            continue
        data.setdefault(s, {'t': [], 'mean': [], 'std': []})
        data[s]['t'].append(float(row['time_hours']))
        data[s]['mean'].append(float(row['mean_circularity']))
        data[s]['std'].append(float(row['std_circularity']))

fig, ax = plt.subplots(figsize=(8, 5))

for s in SELECTED_STIFFNESSES:
    d = data[s]
    t = np.array(d['t'])
    mean = np.array(d['mean'])
    std = np.array(d['std'])
    ax.errorbar(t, mean, yerr=std, fmt='x', color=COLOURS[s],
                label=f'{s:.0f} Pa', capsize=3, markersize=5, linewidth=0)

ax.set_xlabel('Time (hours)')
ax.set_ylabel(r'Circularity ($4\pi A / P^2$)')
ax.set_ylim(0, 1.05)
ax.legend(title='ECM Shear\nModulus', loc='best')

fig.tight_layout()
fig.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out_path}")
