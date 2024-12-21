from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import argparse

plt.rcParams.update({
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'lines.linewidth': 2,
    'font.family': 'serif',
    'font.serif': 'Helvetica',
    'figure.dpi': 1000
})

parser = argparse.ArgumentParser()
parser.add_argument("--i", type=str, required=True)
parser.add_argument("--o", type=str, required=True)
args = parser.parse_args()

tree = Phylo.read(args.i, "newick")

def format_ticks(x, pos):
    return f"{x:.7f}"

fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(1, 1, 1)

Phylo.draw(tree, axes=ax, do_show=False, branch_labels=lambda c: f"{c.branch_length:.7f}" if c.branch_length else "")

ax.xaxis.set_major_formatter(FuncFormatter(format_ticks))

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

for text in ax.texts:
    text.set_fontsize(20)
    x, y = text.get_position()
    text.set_position((x, y + 0.01))

ax.set_title("TDP-43", fontsize=40)
ax.set_xlabel("Branch Length", fontsize=25)
ax.set_ylabel("Taxa", fontsize=25)

plt.tight_layout()
plt.savefig(args.o)