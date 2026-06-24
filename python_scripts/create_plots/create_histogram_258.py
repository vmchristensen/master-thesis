import pandas as pd
import matplotlib.pyplot as plt

# PATHS
csv_path = 'C:\\UPC\\Thesis\\code\\master-thesis\\csv_files\\'
img_path = 'C:\\UPC\\Thesis\\code\\master-thesis\\images\\'

# LOAD
df = pd.read_csv(csv_path + 'demographics_cleaned_filtered_258.csv')

# ── Color scheme ──────────────────────────────────────────────
LEFT_COLOR      = '#D64E4E'  # red
RIGHT_COLOR     = "#2B81AC"  # dark blue
BILATERAL_COLOR = "#7A51A0"  # dark purple
TYPICAL_COLOR   = '#D64E4E'  # red
ATYPICAL_COLOR  = '#5369A0'  # dark periwinkle

df = df[df['ID'] != 'EZM30214']

# ── Data ──────────────────────────────────────────────────────
typical_counts    = df['HEMISPHERIC'].value_counts()
laterality_counts = df['HEMISPHERIC_2'].value_counts()

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
fig.patch.set_facecolor('white')

# ── Left subplot: Typical vs Atypical ─────────────────────────
ax1 = axes[0]
ax1.set_facecolor('white')

categories1 = ['Typical', 'Atypical']
counts1     = [typical_counts.get('Typical', 0), typical_counts.get('Atypical', 0)]
colors1     = [TYPICAL_COLOR, ATYPICAL_COLOR]
total1      = sum(counts1)

bars1 = ax1.bar(categories1, counts1, color=colors1, width=0.5, edgecolor='white', linewidth=1.2)

for bar, count in zip(bars1, counts1):
    ax1.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + 3,
             f'{count}\n({count/total1*100:.1f}%)',
             ha='center', va='bottom', fontsize=11, fontweight='bold', color='#2D2D2D')

ax1.set_title('Language Organization\n(Typical vs Atypical)', fontsize=13, fontweight='bold', pad=12)
ax1.set_ylabel('Count', fontsize=11)
ax1.set_ylabel('Participants (n)', fontsize=11)
ax1.set_ylim(0, max(counts1) * 1.2)
ax1.spines[['top', 'right']].set_visible(False)
ax1.tick_params(axis='both', labelsize=10)

# ── Right subplot: Left / Right / Bilateral ───────────────────
ax2 = axes[1]
ax2.set_facecolor('white')

categories2 = ['Left', 'Bilateral', 'Right']
counts2     = [laterality_counts.get('Left', 0),
               laterality_counts.get('Bilateral', 0),
               laterality_counts.get('Right', 0)]
colors2     = [LEFT_COLOR, BILATERAL_COLOR, RIGHT_COLOR]
bars2 = ax2.bar(categories2, counts2, color=colors2, width=0.5, edgecolor='white', linewidth=1.2)

for bar, count in zip(bars2, counts2):
    ax2.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + 3,
             f'{count}\n({count/total1*100:.1f}%)',
             ha='center', va='bottom', fontsize=11, fontweight='bold', color='#2D2D2D')

ax2.set_title('Language Lateralization\n(Left / Bilateral / Right)', fontsize=13, fontweight='bold', pad=12)
ax2.set_ylabel('Participants (n)', fontsize=11)
ax2.set_ylim(0, max(counts2) * 1.2)
ax2.spines[['top', 'right']].set_visible(False)
ax2.tick_params(axis='both', labelsize=10)

# ── Overall title ─────────────────────────────────────────────
fig.suptitle('Distribution of Language Organization and Lateralization',
             fontsize=15, fontweight='bold', y=1.01)

plt.tight_layout()
plt.savefig(img_path + 'hemispheric_dominance_258.png', dpi=150, bbox_inches='tight')
plt.show()