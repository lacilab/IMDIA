import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Arguments
hbonds_ext = sys.argv[1]
hbonds_prot_domains = sys.argv[2]
hbonds_prot_monomer = sys.argv[3]
hbonds_dna = sys.argv[4]
hbonds_dna_bb = sys.argv[5]
hbonds_dna_base = sys.argv[6]

hbonds_df = pd.read_csv(hbonds_ext)

##########  Protein acceptors and donors in the whole dimer
donor_df = hbonds_df[hbonds_df['donor_type']=='protein']
acceptor_df = hbonds_df[hbonds_df['acceptor_type']=='protein']

plt.figure(figsize=(10, 6))

# Plot donor residues in red
plt.scatter(donor_df['donor_res_num_adjusted'], donor_df['frame%'], 
           color='red', label='Donor', alpha=0.6)

# Plot acceptor residues in blue (if any have 3-letter codes)
plt.scatter(acceptor_df['acceptor_res_num_adjusted'], acceptor_df['frame%'], 
           color='blue', label='Acceptor', alpha=0.6)

plt.xlabel('Residue Number (Adjusted)')
plt.ylabel('Frame %')
plt.title('STAT1 Donor and Acceptor Residues')
plt.grid(True, alpha=0.3)
# plt.axvspan(0, 180, zorder=0, alpha=0.2, color='blue', label='CCD')
# plt.axvspan(181, 352, zorder=0, alpha=0.2, color='red', label='DBD')
# plt.axvspan(353, 440, zorder=0, alpha=0.2, color='green', label='LD')
# plt.axvspan(441, 548, zorder=0, alpha=0.2, color='yellow', label='SH2')
# plt.axvspan(549, 574, zorder=0, alpha=0.2, color='pink', label='Tail')
plt.axvspan(136, 317, zorder=0, alpha=0.2, color='blue', label='CCD')
plt.axvspan(318, 488, zorder=0, alpha=0.2, color='red', label='DBD')
plt.axvspan(489, 576, zorder=0, alpha=0.2, color='green', label='LD')
plt.axvspan(577, 683, zorder=0, alpha=0.2, color='yellow', label='SH2')
plt.axvspan(684, 708, zorder=0, alpha=0.2, color='pink', label='Tail')
plt.legend()
plt.savefig(hbonds_prot_domains)
plt.close()


##########  Protein acceptors and donors per monomer
plt.figure(figsize=(10, 6))

# Plot donor residues in red
plt.scatter(donor_df['donor'], donor_df['frame%'], 
           color='red', label='Donor', alpha=0.6)

# Plot acceptor residues in blue
plt.scatter(acceptor_df['acceptor'], acceptor_df['frame%'], 
           color='blue', label='Acceptor', alpha=0.6)

plt.xlabel('Atom Number')
plt.ylabel('Frame %')
plt.title('STAT1 Donor and Acceptor Residues')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig(hbonds_prot_monomer)
plt.close()


##########  dsDNA acceptors and donors 
# Filter for DNA donors and acceptors
dna_donor_df = hbonds_df[hbonds_df['donor_type'] == 'DNA']
dna_acceptor_df = hbonds_df[hbonds_df['acceptor_type'] == 'DNA']

plt.figure(figsize=(10, 6))

# Plot DNA donor residues in red using donor atom numbers
plt.scatter(dna_donor_df['donor'], dna_donor_df['frame%'], 
           color='red', label='Donor', alpha=0.6)

# Plot DNA acceptor residues in blue using acceptor atom numbers
plt.scatter(dna_acceptor_df['acceptor'], dna_acceptor_df['frame%'], 
           color='blue', label='Acceptor', alpha=0.6)

plt.xlabel('DNA Residue Number')
plt.ylabel('Frame %')
plt.title('DNA Donor and Acceptor Residues')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig(hbonds_dna)
plt.close()


##########  dsDNA backbone acceptors and donors
# Filter for DNA backbone donors and acceptors
donor_df_backbone = hbonds_df[hbonds_df['donor_section'] == 'backbone']
acceptor_df_backbone = hbonds_df[hbonds_df['acceptor_section'] == 'backbone']

# Create the plot
plt.figure(figsize=(10, 6))

# Plot donor residues in red (backbone only)
plt.scatter(donor_df_backbone['donor'], donor_df_backbone['frame%'], 
            color='red', label='Donor', alpha=0.6)

# Plot acceptor residues in blue (backbone only)
plt.scatter(acceptor_df_backbone['acceptor'], acceptor_df_backbone['frame%'], 
            color='blue', label='Acceptor', alpha=0.6)

# Customize the plot
plt.xlabel('Atom Number')
plt.ylabel('Frame %')
plt.title('STAT1 Donor and Acceptor Residues (DNA Backbone Only)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig(hbonds_dna_bb)
plt.close()


##########  dsDNA bases acceptors and donors
# Filter for DNA base donors and acceptors
donor_df_backbone = hbonds_df[hbonds_df['donor_section'] == 'base']
acceptor_df_backbone = hbonds_df[hbonds_df['acceptor_section'] == 'base']

# Create the plot
plt.figure(figsize=(10, 6))

# Plot donor residues in red (backbone only)
plt.scatter(donor_df_backbone['donor'], donor_df_backbone['frame%'], 
            color='red', label='Donor', alpha=0.6)

# Plot acceptor residues in blue (backbone only)
plt.scatter(acceptor_df_backbone['acceptor'], acceptor_df_backbone['frame%'], 
            color='blue', label='Acceptor', alpha=0.6)

# Customize the plot
plt.xlabel('Atom Number')
plt.ylabel('Frame %')
plt.title('STAT1 Donor and Acceptor Residues (DNA Bases Only)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig(hbonds_dna_base)
plt.close()
