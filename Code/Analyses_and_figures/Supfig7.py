import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt


no_native_merged = pd.read_table(path_to_sup7_data)

fig = plt.figure(figsize=(10, 10), constrained_layout=True)
gs1 = plt.GridSpec(nrows=2, ncols=2)
# ax_array = fig.subplots(2, 2,  )#sharex=True, sharey='row')
ax_array = []

sets = ['motif', 'non_motif']
methods = ['avg_lddt', 'norm_to_max_model']
arg_dict = {
	'motif': {
		'linker': {'c': 'b', 'marker': 'x', 'label':'linker'},
		'sep': {'c': 'deepskyblue', 'marker': '+', 'label':'sep chains'},
	},
	'non_motif': {
		'linker': {'c': 'r', 'marker': 'x', 'label':'linker'},
		'sep': {'c': 'coral', 'marker': '+', 'label':'sep chains'},
	}
}
for i, set in enumerate(sets):
	for j, method in enumerate(methods):
		plt_tab = no_native_merged.query('set == @set')
		if i==1 and j==0:
			ax = fig.add_subplot(gs1[j,i], sharey=ax_array[0])
		else:
			ax = fig.add_subplot(gs1[j,i])
		for subset in ['linker', 'sep']:
			# ax_array[j, i].scatter(
			ax.scatter(
				x=plt_tab.query('run_type == @subset')['dockq'],
				y=plt_tab.query('run_type == @subset')[method],
				**arg_dict[set][subset]
			)
		ax_array.append(ax)

font_size = 12
ax_array[0].set_title('Motif Set', fontsize=font_size)
ax_array[2].set_title('Non-motif Set', fontsize=font_size)
ax_array[0].set_ylabel('Average pLDDT (peptide)', fontsize=font_size)
ax_array[1].set_ylabel('Buried Surface Area [Å²]', fontsize=font_size)
ax_array[0].legend()
ax_array[2].legend()
ax_array[0].axhline(y=0.7, xmin=0, xmax=1, linestyle='--', lw=0.9, c='black')
ax_array[2].axhline(y=0.7, xmin=0, xmax=1, linestyle='--', lw=0.9, c='black')
ax_array[1].set_xlabel("DockQ", fontsize=font_size)
ax_array[3].set_xlabel("DockQ", fontsize=font_size)

plt.show()

