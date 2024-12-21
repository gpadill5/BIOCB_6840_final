import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Amino acid color map 
aa_colors = {
    "A": "#FF7F0E", "R": "#1F77B4", "N": "#2CA02C", "D": "#D62728",
    "C": "#9467BD", "Q": "#8C564B", "E": "#E377C2", "G": "#7F7F7F",
    "H": "#BCBD22", "I": "#17BECF", "L": "#FFBB78", "K": "#98DF8A",
    "M": "#C5B0D5", "F": "#DBDB8D", "P": "#9EDAE5", "S": "#AD494A",
    "T": "#E7969C", "W": "#C49C94", "Y": "#F7B6D2", "V": "#C7C7C7"
}

# Input your known/predicted IDR range(s)
# String must be label in alignment_file!
# int input sourced from sequences folder
disordered_regions = {
    "FUS_HUMAN": [(1, 278), (368, 426), (451, 526)],
    "FUS_MOUSE": [(1, 271), (361, 419), (444, 518)],
    "FUS_XENLA": [(1, 291), (383, 439), (463, 536)],
    "FUS_PANTR": [(1, 278), (370, 426), (450, 526)],
    "FUS_BOVIN": [(1, 265), (355, 414), (438, 513)]
}

# Parse clustal alignment file and ensure padding
def parse_clustal_alignment(file_path):
    species_sequences = {}
    max_length = 0  # Track the longest sequence 
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("CLUSTAL") or line.startswith("STOCKHOLM"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                species_id, sequence_fragment = parts[0], parts[1]
                # Skip lines where species_id contains no alphabetic characters (e.g., consensus lines)
                if not any(c.isalpha() for c in species_id):
                    continue
                species_sequences[species_id] = species_sequences.get(species_id, "") + sequence_fragment
                max_length = max(max_length, len(species_sequences[species_id]))
    # Pad all sequences to the same length with gaps
    for species_id in species_sequences:
        sequence = species_sequences[species_id]
        if len(sequence) < max_length:
            species_sequences[species_id] += "-" * (max_length - len(sequence))
    return species_sequences

# Calculate consensus line
def calculate_consensus(species_sequences):
    consensus = []
    num_species = len(species_sequences)
    alignment_length = len(next(iter(species_sequences.values())))  # Assume all sequences are the same length
    sequences = list(species_sequences.values())

    for col in range(alignment_length):
        column_residues = [seq[col] for seq in sequences if seq[col] != "-"]
        if len(column_residues) == 0:
            consensus.append(" ")
            continue
        most_common_residue = max(set(column_residues), key=column_residues.count)
        count = column_residues.count(most_common_residue)
        if count == num_species:
            consensus.append("*")
        elif count > num_species / 2:
            consensus.append(":")
        elif count > 1:
            consensus.append(".")
        else:
            consensus.append(" ")
    return "".join(consensus)

def lighten_color(color, amount=0.5):
    import matplotlib.colors as mc
    try:
        c = mc.cnames[color]
    except KeyError:
        c = color
    c = mc.to_rgb(c)
    c = [c[i] + (1 - c[i]) * amount for i in range(3)]
    return c

def plot_alignment_with_consensus(species_sequences, disordered_regions, wrap_length=60):
    consensus_line = calculate_consensus(species_sequences)

    boundary_counts = {species: 0 for species in species_sequences}

    species_prev_state = {species: None for species in species_sequences}
    current_y = 0
    min_y = current_y 
    num_sequences = len(species_sequences)
    num_wraps = (len(consensus_line) + wrap_length - 1) // wrap_length  
    fig_height = max(6, num_sequences * num_wraps * 0.6)  
    fig, ax = plt.subplots(figsize=(15, fig_height))

    alignment_length = len(consensus_line)

    aa_indices = {}  # species -> list of (aa_pos, idx_in_alignment)
    for species, seq in species_sequences.items():
        indices = []
        aa_pos = 0  # amino acid position (excluding gaps)
        for idx, res in enumerate(seq):
            if res != '-':
                aa_pos += 1
                indices.append((aa_pos, idx))  # (aa_pos, alignment_idx)
            else:
                indices.append((None, idx))  
        aa_indices[species] = indices

    for wrap_start in range(0, alignment_length, wrap_length):
        wrapped_consensus = consensus_line[wrap_start:wrap_start + wrap_length]
        wrapped_sequences = {
            species: seq[wrap_start:wrap_start + wrap_length]
            for species, seq in species_sequences.items()
        }

        for j, res in enumerate(wrapped_consensus):
            ax.text(j + 0.5, current_y + 0.5, res, ha="center", va="center", color="#000000", fontsize=12, family="monospace")
        current_y -= 0.75  
        min_y = min(min_y, current_y)

        for species, sequence in wrapped_sequences.items():
            indices = aa_indices[species][wrap_start:wrap_start + wrap_length]
            is_disordered = []
            aa_pos_list = []  
            for aa_pos, idx_in_alignment in indices:
                if aa_pos is not None:
                    disordered = False
                    if species in disordered_regions:
                        for start_aa, end_aa in disordered_regions[species]:
                            if start_aa <= aa_pos <= end_aa:
                                disordered = True
                                break
                    is_disordered.append(disordered)
                else:
                    disordered = False
                    if species in disordered_regions:
                        prev_aa_pos = aa_pos_list[-1] if aa_pos_list else None
                        next_idx = wrap_start + len(aa_pos_list) + 1
                        next_aa_pos = None
                        if next_idx < len(aa_indices[species]):
                            for future_aa_pos, future_idx in aa_indices[species][next_idx:]:
                                if future_aa_pos is not None:
                                    next_aa_pos = future_aa_pos
                                    break
                        for start_aa, end_aa in disordered_regions[species]:
                            if (prev_aa_pos is not None and start_aa <= prev_aa_pos <= end_aa) or \
                               (next_aa_pos is not None and start_aa <= next_aa_pos <= end_aa):
                                disordered = True
                                break
                    is_disordered.append(disordered)
                aa_pos_list.append(aa_pos)

            # Plot the sequence
            for j, res in enumerate(sequence):
                if res != '-':
                    color = aa_colors.get(res, "#000000") 
                    if not is_disordered[j]:
                        color = lighten_color(color, amount=0.5)
                    ax.add_patch(Rectangle((j, current_y), 1, 1, color=color))
                    ax.text(j + 0.5, current_y + 0.5, res, ha="center", va="center", fontsize=10, family="monospace")
                else:
                    # Just plot the '-' character without background for gaps
                    ax.text(j + 0.5, current_y + 0.5, res, ha="center", va="center", fontsize=10, family="monospace")
            ax.text(-1, current_y + 0.5, species, ha="right", va="center", fontsize=10, family="monospace")

            for j in range(len(is_disordered)):
                current_state = is_disordered[j]
                if species_prev_state[species] is not None and current_state != species_prev_state[species]:
                    boundary_counts[species] += 1
                species_prev_state[species] = current_state

            for j in range(1, len(is_disordered)):
                if is_disordered[j] != is_disordered[j - 1]:
                    # Draw a vertical bold black line at this position, use photoshop to put horizontal lines
                    x = j
                    y = current_y
                    ax.plot([x, x], [y, y + 1], color='black', linewidth=2)
            current_y -= 1
            min_y = min(min_y, current_y)

        current_y -= 1.5 
        min_y = min(min_y, current_y)

    ax.set_xlim(-1, wrap_length)
    ax.set_ylim(min_y - 1, 1)
    ax.axis("off") 
    plt.title("FUS Multiple Sequence Alignment (MSA)", fontsize=14)
    plt.tight_layout()
    plt.savefig('fus_MSA.svg')

    # Boarder debug
    print("Number of IDR boarders for each species:")
    for species, count in boundary_counts.items():
        print(f"{species}: {count} boarders")

alignment_file = "/clustal-omega-1.2.4/fus_aligned.clu"
species_sequences = parse_clustal_alignment(alignment_file)
try:
    plot_alignment_with_consensus(species_sequences, disordered_regions)
except Exception as e:
    print(f"An error occurred: {e}")


