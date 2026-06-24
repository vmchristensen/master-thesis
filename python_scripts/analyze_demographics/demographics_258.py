# new demographics on only the 259 subjects


import pandas as pd

# ── PATHS ──────────────────────────────────────────────────────────────────────
csv_path     = 'C:\\UPC\\Thesis\\code\\master-thesis\\csv_files\\'
results_path = 'C:\\UPC\\Thesis\\code\\master-thesis\\demo_results\\'

# ── LOAD PRE-FILTERED CSV (259 subjects with images) ──────────────────────────
merged = pd.read_csv(csv_path + 'demographics_cleaned_filtered_258.csv')

print(f"Total subjects loaded: {len(merged)}")
print(f"EZM2: {(merged['cohort'] == 2).sum()}  |  EZM3: {(merged['cohort'] == 3).sum()}")

df2 = merged[merged['cohort'] == 2].reset_index(drop=True)
df3 = merged[merged['cohort'] == 3].reset_index(drop=True)


# ── METRICS FUNCTION ──────────────────────────────────────────────────────────
def compute_metrics(df):
    metrics = {
        # LI_GEN
        'LI_GEN_min':   df['LI_GEN'].min(),
        'LI_GEN_max':   df['LI_GEN'].max(),
        'LI_GEN_range': f"{df['LI_GEN'].min()}, {df['LI_GEN'].max()}",
        'LI_GEN_mean':  df['LI_GEN'].mean(),
        'LI_GEN_std':   df['LI_GEN'].std(),

        # Laterality Quotient (Edinburgh)
        'LateralityQ_min':   df['Laterality_Quotient'].min(),
        'LateralityQ_max':   df['Laterality_Quotient'].max(),
        'LateralityQ_range': f"{df['Laterality_Quotient'].min()}, {df['Laterality_Quotient'].max()}",
        'LateralityQ_mean':  df['Laterality_Quotient'].mean(),
        'LateralityQ_std':   df['Laterality_Quotient'].std(),

        # Edinburgh score
        'Edinburgh_min':   df['Puntuacion.Edinburgh'].min(),
        'Edinburgh_max':   df['Puntuacion.Edinburgh'].max(),
        'Edinburgh_range': f"{df['Puntuacion.Edinburgh'].min()}, {df['Puntuacion.Edinburgh'].max()}",
        'Edinburgh_mean':  df['Puntuacion.Edinburgh'].mean(),
        'Edinburgh_std':   df['Puntuacion.Edinburgh'].std(),

        # Age
        'Age_min':   df['Edad'].min(),
        'Age_max':   df['Edad'].max(),
        'Age_range': f"{df['Edad'].min()}, {df['Edad'].max()}",
        'Age_mean':  df['Edad'].mean(),
        'Age_std':   df['Edad'].std(),

        # Typical vs Atypical
        'Typical_count':      (df['HEMISPHERIC'] == 'Typical').sum(),
        'Atypical_count':     (df['HEMISPHERIC'] == 'Atypical').sum(),
        'Typical_percentage': f"{(df['HEMISPHERIC'] == 'Typical').mean() * 100:.2f}%",
        'Atypical_percentage':f"{(df['HEMISPHERIC'] == 'Atypical').mean() * 100:.2f}%",

        # Left / Bilateral / Right
        'Left_count':        (df['HEMISPHERIC_2'] == 'Left').sum(),
        'Bilateral_count':   (df['HEMISPHERIC_2'] == 'Bilateral').sum(),
        'Right_count':       (df['HEMISPHERIC_2'] == 'Right').sum(),
        'Percentage_Left':   f"{(df['HEMISPHERIC_2'] == 'Left').mean() * 100:.2f}%",
        'Percentage_Bilateral': f"{(df['HEMISPHERIC_2'] == 'Bilateral').mean() * 100:.2f}%",
        'Percentage_Right':  f"{(df['HEMISPHERIC_2'] == 'Right').mean() * 100:.2f}%",

        # Sex
        'Male_count':             (df['Sexo'] == 'Hombre').sum(),
        'Female_count':           (df['Sexo'] == 'Mujer').sum(),
        'Male_percentage':        f"{(df['Sexo'] == 'Hombre').mean() * 100:.2f}%",
        'Female_percentage':      f"{(df['Sexo'] == 'Mujer').mean() * 100:.2f}%",

        # Handedness subgroups (Bryden 1977 cutoffs on Puntuacion.Edinburgh)
        # 24–37 = mixed-handed, 38–50 = left-handed
        'Lefthanded_count':       ((df['Puntuacion.Edinburgh'] >= 38) & (df['Puntuacion.Edinburgh'] <= 50)).sum(),
        'Mixedhanded_count':      ((df['Puntuacion.Edinburgh'] >= 24) & (df['Puntuacion.Edinburgh'] <= 37)).sum(),
        'Lefthanded_percentage':  f"{((df['Puntuacion.Edinburgh'] >= 38) & (df['Puntuacion.Edinburgh'] <= 50)).mean() * 100:.2f}%",
        'Mixedhanded_percentage': f"{((df['Puntuacion.Edinburgh'] >= 24) & (df['Puntuacion.Edinburgh'] <= 37)).mean() * 100:.2f}%",

        # Totals & missingness
        'Total_count':   len(df),
        'Missing_count': df.isnull().sum().to_dict(),
        'Ids_with_missing_musical_formation': (
            df[df['Formacion_musical.Tipo'].isnull()]['ID'].tolist()
            if 'Formacion_musical.Tipo' in df.columns else 'N/A'
        ),
        'Ids_with_missing_any_excl_musical_formation': (
            df[df.drop(columns=['Formacion_musical.Tipo']).isnull().any(axis=1)]['ID'].tolist()
            if 'Formacion_musical.Tipo' in df.columns
            else df[df.isnull().any(axis=1)]['ID'].tolist()
        ),
    }
    return metrics


# ── SAVE METRICS ──────────────────────────────────────────────────────────────
def save_metrics_txt(metrics, filename):
    with open(filename, 'w') as f:
        for key, value in metrics.items():
            try:
                f.write(f"{key}\t{float(value):.2f}\n")
            except (ValueError, TypeError):
                f.write(f"{key}\t{value}\n")


metrics_ezm2   = compute_metrics(df2)
metrics_ezm3   = compute_metrics(df3)
metrics_merged = compute_metrics(merged)

save_metrics_txt(metrics_ezm2,   results_path + 'EZM2_metrics.txt')
save_metrics_txt(metrics_ezm3,   results_path + 'EZM3_metrics.txt')
save_metrics_txt(metrics_merged, results_path + 'merged_metrics.txt')

print("Done — metrics saved to results folder.")