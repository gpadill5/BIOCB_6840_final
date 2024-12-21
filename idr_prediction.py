import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import gaussian_filter1d
import os

plt.rcParams.update({
    'font.size': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'lines.linewidth': 2,
    'font.family': 'serif',
    'font.serif': 'Helvetica',
    'figure.dpi': 1000
})

parser = argparse.ArgumentParser()
parser.add_argument("--i", type=str, required=True)
parser.add_argument("--o", type=str, required= True)
parser.add_argument("--s", type=str, required= True)
args = parser.parse_args()

if not os.path.exists(args.i):
    raise FileNotFoundError(f"File not found: {args.i}")

data_file = args.i
df = pd.read_csv(data_file, sep="\t")

pLDDT_threshold = 10
RSA_threshold = 0.581
sigma = 7
df['lddt_smoothed'] = gaussian_filter1d(df['lddt'], sigma)
df['rsa_smoothed'] = gaussian_filter1d(df['rsa'], sigma)

df['is_disordered'] = (df['lddt_smoothed'] < pLDDT_threshold) & (df['rsa_smoothed'] > RSA_threshold)
df['group'] = (df['is_disordered'] != df['is_disordered'].shift()).cumsum()

groups = df.groupby(['group', 'is_disordered'])['pos'].agg(['min', 'max']).reset_index()

length = 1
groups['length'] = groups['max'] - groups['min'] + 1
filtered_disordered_groups = groups[(groups['is_disordered']) & (groups['length'] > length)]

max_residue_position = df['pos'].max()

plt.figure(figsize=(12, 4))

plt.axhspan(0, 1, xmin=0, xmax=1, facecolor='blue', alpha=0.5, label="Ordered")

for _, row in filtered_disordered_groups.iterrows():
    start = row['min'] - 0.5
    end = row['max'] + 0.5
    plt.axvspan(start, end, facecolor='red', alpha=0.5, label="Disordered")

plt.title(f"Disordered Regions in Prediction Based on Smoothed RSA", pad=20)
plt.xlabel("Residue Position")
plt.ylabel("")
plt.ylim(0, 1)
plt.xlim(0, max_residue_position + 1)
plt.yticks([])
plt.legend(handles=[
    plt.Line2D([0], [0], color='red', lw=4, label='Disordered'),
    plt.Line2D([0], [0], color='blue', lw=4, label='Ordered')
])

ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=20))

plt.tight_layout()
plt.savefig(args.o)

print(f"Filtered disordered regions (length > {length}):")
for _, row in filtered_disordered_groups.iterrows():
    print(f"Positions {int(row['min'])} to {int(row['max'])}, Length: {int(row['length'])}")
