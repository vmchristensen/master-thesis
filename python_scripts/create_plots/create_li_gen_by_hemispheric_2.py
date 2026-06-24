import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# PATHS
csv_path = 'C:\\UPC\\Thesis\\code\\master-thesis\\csv_files\\'
img_path  = 'C:\\UPC\\Thesis\\code\\master-thesis\\images\\'

# LOAD
df = pd.read_csv(csv_path + 'demographics_cleaned_filtered_258.csv')

# ── Color scheme ──────────────────────────────────────────────
LEFT_COLOR      = '#D64E4E'
RIGHT_COLOR     = '#2B81AC'
BILATERAL_COLOR = '#7A51A0'

# ── Plot ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

bins = np.linspace(df['LI_GEN'].min(), df['LI_GEN'].max(), 30)

for group, color, label in [('Left', LEFT_COLOR, 'Left'),
                              ('Right', RIGHT_COLOR, 'Right'),
                              ('Bilateral', BILATERAL_COLOR, 'Bilateral')]:
    subset = df[df['HEMISPHERIC_2'] == group]['LI_GEN'].dropna()
    ax.hist(subset, bins=bins, color=color, alpha=0.75, label=label,
            edgecolor='white', linewidth=0.6)

ax.axvline(0, color='#888888', linestyle='--', linewidth=1.2, alpha=0.7)

ax.set_title('LI_GEN Distribution by Lateralization Group',
             fontsize=13, fontweight='bold', pad=12)
ax.set_xlabel('LI_GEN', fontsize=11)
ax.set_ylabel('Participants (n)', fontsize=11)
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(axis='both', labelsize=10)
ax.legend(fontsize=10, frameon=False)
ax.invert_xaxis()

plt.tight_layout()
plt.savefig(img_path + 'LI_GEN_distribution_258.png', dpi=150, bbox_inches='tight')
plt.show()