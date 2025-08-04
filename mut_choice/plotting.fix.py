import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker

# Define file paths
abonly_file_path = "/home/luna.kuleuven.be/u0147150/Desktop/PhD/WP2_Overexpression/1_Insilico/2_FoldX/Analysis/STAT1GOF/Outcom_DNA_Int/merged.dna.txt"
gof_mutations_file_path = "/home/luna.kuleuven.be/u0147150/Desktop/PhD/WP2_Overexpression/1_Insilico/2_FoldX/KnownMutations/known_GOFmut"

# Define correct column headers for abonly data
correct_headers_abonly = [
    "Pdb", "Group1", "Group2", "IntraclashesGroup1", "IntraclashesGroup2",
    "Interaction Energy", "StabilityGroup1", "StabilityGroup2"
]

# Load the abonly data while skipping redundant headers
data_abonly = pd.read_csv(abonly_file_path, sep="\t", header=0, names=correct_headers_abonly, comment='P')

# Drop potential NaN rows and reset index
data_abonly = data_abonly.dropna().reset_index(drop=True)

# Extract MutationIndex from Pdb filenames (e.g., ./1BF5_Dimer_Repair_10_1.pdb → 10)
data_abonly["MutationIndex"] = data_abonly["Pdb"].str.extract(r'_(\d+)_\d+\.pdb')[0].astype(float)

# Convert "Interaction Energy" to numeric
data_abonly["Interaction Energy"] = pd.to_numeric(data_abonly["Interaction Energy"], errors='coerce')

# Load known GOF mutations data without headers
updated_gof_mutations = pd.read_csv(gof_mutations_file_path, sep="\t", header=None)

# Assign column names manually
updated_gof_mutations.columns = ["MutationIndex", "Mutation", "Domain"]

# Extract ResidueIndex from Mutation names (e.g., "T419R" -> 419)
updated_gof_mutations["ResidueIndex"] = updated_gof_mutations["Mutation"].str.extract(r'(\d+)')[0].astype(float)

# Merge extracted MutationIndex with ResidueIndex from known mutations
data_annotated_abonly = data_abonly.merge(updated_gof_mutations[["MutationIndex", "ResidueIndex", "Mutation", "Domain"]], on="MutationIndex", how="left")

# Sort data by ResidueIndex for correct scaling
data_annotated_abonly = data_annotated_abonly.sort_values(by="ResidueIndex")

# Define residue index ranges for domains
domain_ranges = {
    "ND": (1, 133),
    "CCD": (134, 316),
    "DBD": (317, 487),
    "LD": (488, 575),
    "SH2": (576, 682),
    "TAD": (683, 750)
}

# Define colors for domains
domain_colors = {
    "ND": "lightblue",
    "CCD": "blue",
    "DBD": "green",
    "LD": "orange",
    "SH2": "red",
    "TAD": "purple"
}

# Create the properly scaled plot
fig, ax = plt.subplots(figsize=(16, 6))

# Set y-axis limits with padding
y_min = data_annotated_abonly["Interaction Energy"].min() - 0.5
y_max = data_annotated_abonly["Interaction Energy"].max() + 0.5
ax.set_ylim(y_min, y_max)

# Set x-axis limits using true numerical scale of ResidueIndex
x_min = data_annotated_abonly["ResidueIndex"].min() - 5
x_max = data_annotated_abonly["ResidueIndex"].max() + 5
ax.set_xlim(x_min, x_max)

# Add background colors based on domain ranges
for domain, (start, end) in domain_ranges.items():
    ax.add_patch(patches.Rectangle(
        (start, y_min),  # Start x position and y lower bound
        end - start,  # Width
        y_max - y_min,  # Height
        color=domain_colors[domain], alpha=0.3, label=domain if domain not in ax.get_legend_handles_labels()[1] else ""
    ))

# Plot interaction energy data points **with proper ResidueIndex scaling**
ax.scatter(data_annotated_abonly["ResidueIndex"], data_annotated_abonly["Interaction Energy"], color="black", label="Interaction Energy", s=10)

# Ensure **all ResidueIndex labels align correctly**
unique_residues = data_annotated_abonly.drop_duplicates(subset=["ResidueIndex"])
ax.set_xticks(unique_residues["ResidueIndex"])
ax.set_xticklabels(unique_residues["Mutation"], rotation=90, fontsize=8)

# Label axes
ax.set_xlabel("Residue Index (Numerically Scaled)", fontsize=12)
ax.set_ylabel("Interaction Energy (kcal/mol)", fontsize=12)
ax.set_title("Interaction Energy Across Mutations - Residue Index Scaled X-Axis", fontsize=14, fontweight='bold')
ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)  # Light grid for clarity

# Add legend
ax.legend()

# Show the plot
plt.show()

