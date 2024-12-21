import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'lines.linewidth': 2,
    'font.family': 'serif',
    'font.serif': 'Helvetica',
    'figure.dpi': 1000
})

# Use alphafold_disorder_plot.py to generate IDR ranges for different \sigma values
data = {
    "TDP-43_human": [(76, 99), (261, 414)],
    "TDP-43_human_7": [(78, 100), (262, 414)],
    "TDP-43_human_5": [(79, 100), (182, 185), (262, 414)],
    "TDP-43_human_3": [(79, 100), (180, 187), (262, 414)],
    "TDP-43_human_1": [(1, 2), (11, 14), (22, 23), (47, 48), (65, 66), (78, 101),
                       (139, 142), (179, 187), (224, 225), (261, 414)],
    "FUS_human": [(1, 279), (368, 417), (437, 518)],
    "FUS_human_7": [(1, 279), (370, 425), (452, 526)],
    "FUS_human_5": [(1, 280), (371, 424), (452, 526)],
    "FUS_human_3": [(1, 280), (315, 317), (328, 330), (371, 423), (452, 526)],
    "FUS_human_1": [(1, 163), (165, 185), (187, 280), (314, 318), (328, 331),
                    (371, 423), (431, 432), (453, 512), (514, 526)],
    "Axin-1_human": [(1, 78), (216, 393), (402, 782)],
    "Axin-1_human_7": [(1, 78), (216, 392), (403, 780)],
    "Axin-1_human_5": [(1, 78), (216, 390), (404, 730), (735, 779)],
    "Axin-1_human_3": [(1, 79), (216, 388), (404, 471), (477, 635), (642, 729),
                     (737, 779), (830, 833)],
    "Axin-1_human_1": [(1, 78), (127, 132), (174, 175), (214, 216), (219, 357),
                     (360, 387), (394, 396), (406, 467), (470, 472), (478, 627),
                     (629, 636), (644, 708), (711, 714), (716, 729), (737, 780),
                     (792, 794), (816, 818), (829, 833), (842, 844), (861, 862)]
}

# Function to calculate the Jaccard index
def calculate_jaccard_index(reference_intervals, comparison_intervals):
    reference_set = set()
    for start, end in reference_intervals:
        reference_set.update(range(start, end))

    comparison_set = set()
    for start, end in comparison_intervals:
        comparison_set.update(range(start, end))
    intersection = reference_set & comparison_set
    union = reference_set | comparison_set

 
    return len(intersection) / len(union) if union else 0


def plot_intervals_with_jaccard(reference_label, reference_data, comparisons):
    fig, ax = plt.subplots(figsize=(12, 6))

    y_labels = []
    y_ticks = []
    y_pos = 0 

    # Plot reference intervals
    y_labels.append('MobiDB AlphaFold')
    y_ticks.append(y_pos + 5)
    ax.broken_barh([(start, end - start) for start, end in reference_data],
                   (y_pos, 10), facecolors='blue', edgecolors='black', label='Reference')
    y_pos += 15 

    for label, intervals in comparisons.items():
        # Extract "i" from "n_human_i"
        sigma_value = label.split("_")[-1] if label.split("_")[-1].isdigit() else "Ref"
        # Calculate Jaccard index
        jaccard_index = calculate_jaccard_index(reference_data, intervals)

        # Update label with \sigma=i and Jaccard index
        label_with_index = f"σ={sigma_value} (J={jaccard_index:.2f})"
        y_labels.append(label_with_index)
        y_ticks.append(y_pos + 5)
        ax.broken_barh([(start, end - start) for start, end in intervals],
                       (y_pos, 10), facecolors='black', edgecolors='black')
        y_pos += 15

    reference_patch = mpatches.Patch(color='blue', label='Reference')
    comparison_patch = mpatches.Patch(color='black', label='Gaussian Smoothed RSA')


    ax.legend(
        handles=[reference_patch, comparison_patch],
        loc='upper left',        
        bbox_to_anchor=(1.05, 1),  
        fontsize='small',      
        framealpha=0.8           
    )


    ax.set_xlabel('Residue Position')
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_title(f'IDR Labeling Comparison for $H. sapiens$ {reference_label}')
    plt.tight_layout()
    plt.savefig(f'{reference_label}_benchmark.svg')


grouped_data = {}
for key, value in data.items():
    base_key = key.split("_")[0]
    if base_key not in grouped_data:
        grouped_data[base_key] = {"reference": value, "comparisons": {}}
    elif key == base_key:
        grouped_data[base_key]["reference"] = value
    else:
        grouped_data[base_key]["comparisons"][key] = value

# Plot each group
for group, datasets in grouped_data.items():
    plot_intervals_with_jaccard(group, datasets["reference"], datasets["comparisons"])
