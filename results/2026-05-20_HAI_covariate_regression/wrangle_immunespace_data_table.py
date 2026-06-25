import sys
import pandas as pd 
import os
import numpy as np
import matplotlib.pyplot as plt
from patsy.contrasts import Treatment
import warnings
warnings.filterwarnings("ignore")

!{sys.executable} -m pip install ipykernel --upgrade --force-reinstall
!{sys.executable} -m pip install statsmodels
!{sys.executable} -m pip install matplotlib

import statsmodels.formula.api as smf

#directories 

SCRATCH = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/"

data = "/sc/arion/work/willij115/projects/HAI/data/2026-05-20_HAI_covariate_regression"
scratch = "/sc/arion/scratch/willij115/projects/HAI/2026-05-20_HAI_covariate_regression"
results = "/sc/arion/work/willij115/projects/HAI/results/2026-05-20_HAI_covariate_regression"

pd.set_option('display.max_columns', None)


#upload immunespace data tables
 
studies = pd.read_csv(os.path.join(SCRATCH, 'immunespaceHAI_studies_tables.csv'), names=['Study ID', 'Study Title', 'Study Data Release', 'PMID', 'Publication Title', 'Publication Date', 'Author Count', 'Authors'], header = 0)
arms = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_arms_tables.csv'),  header = 0)
participants = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_participants_tables.csv'), header = 0)
events = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_events_tables.csv'), header = 0)
assays = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_assays_tables.csv'), header = 0)



merged1 = pd.merge(
    studies,
    arms,
    how='outer',
    on='Study ID',
    suffixes=('_studies', '_arms')
)

merged2 = pd.merge(
    merged1,
    participants,
    how='outer',
    on= ['Study ID', 'Arm ID']
)

merged3 = pd.merge(
    merged2,
    events,
    how='outer',
    on= ['Study ID', 'Participant ID'],
    suffixes=('participants', '_events')         
)

merged4 = pd.merge(
    merged3,
    assays,
    how='outer',
    on= ['Study ID', 'Participant ID', 'Event ID'],
    suffixes=('_events', '_assays')
)


filtered = merged4.dropna(
    subset=['Target Entity Subtype', 'Study ID', 'Arm ID', 'Value' ]
)


##full covariate table (sex, age  + race for studies that have it )

full = filtered.dropna(
    subset=['Target Entity Subtype', 'Study ID', 'Arm ID', 'Value', 'Biological sex', 'Age', 'Race' ]
)

#sex and age only covariate table (for studies that don't have race 

sex_age = filtered = filtered.dropna(
    subset=['Target Entity Subtype', 'Study ID', 'Arm ID', 'Value', 'Biological sex', 'Age' ]
)


strain_summary = (
    filtered.groupby('Target Entity Subtype')['Value']
      .agg(['count', 'median', 'mean', 'std'])
      .sort_values('std', ascending=False)
)


analysis_df= filtered.copy()
analysis_df = filtered[['Value', 'Age', 'Biological sex', "Target Entity Subtype", "Study ID", "Arm ID", "Participant ID", "Event ID", "Count Participants"]].dropna().copy()
analysis_df['Value'] = pd.to_numeric(analysis_df['Value'])
analysis_df['Age'] = pd.to_numeric(analysis_df['Age'])

# Create a new column with log2-transformed HAI titers
analysis_df["Value_log2"] = np.log2(analysis_df["Value"].replace(0, np.nan))

analysis_df = analysis_df[
    analysis_df["Biological sex"].isin(["female", "male"])
]










analysis_df['Value_log2'].describe()

plt.hist(analysis_df['Value_log2'], bins=10, alpha=0.7)
plt.xlabel("Value (log2)")
plt.ylabel("Count")
plt.show()



#strain vs HI
analysis_df.boxplot(
    column="Value_log2",
    by="Target Entity Subtype",
    rot=90,
    figsize=(12, 6)
)
plt.suptitle("")
plt.xlabel("Influenza strain")
plt.ylabel("log2 HAI")
plt.show()



#age vs HI 
plt.figure(figsize=(6,4))
plt.scatter(
    analysis_df["Age"],
    analysis_df["Value_log2"],
    alpha=0.3
)
plt.xlabel("Age")
plt.ylabel("HI titer")
plt.show()

#sex vs HI
analysis_df.boxplot(
    column="Value_log2",
    by="Biological sex",
    figsize=(5,4)
)
plt.suptitle("")
plt.title("HI titer by Sex")
plt.ylabel("HI titer")
plt.show()










# %%
#!/usr/bin/env python3

"""
Minimal script to run HI titer screening analysis + plots.

Usage:
    python run_screening.py

Assumes: analysis_df is loaded in your Python environment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from scipy import stats
import warnings
warnings.filterwarnings("ignore")

print("="*80)
print("STAGE 1: STRAIN SELECTION")
print("="*80)

# Summary statistics by strain
strain_summary = (
    analysis_df.groupby("Target Entity Subtype")
    .agg({
        "Value_log2": ["median", "std", "mean", "count"],
        "Participant ID": "nunique",
        "Study ID": "nunique",
    })
    .round(2)
)

strain_summary.columns = [
    "Median_log2", "Std_log2", "Mean_log2", "N_obs",
    "N_participants", "N_studies"
]

# Filter: at least 20 unique participants
strain_summary = strain_summary[
    strain_summary["N_participants"] >= 20
].copy()

strain_summary = strain_summary.sort_values("Std_log2", ascending=False)

print(f"\n✓ Found {len(strain_summary)} strains with >=20 participants\n")
print("Top 10 strains by HI titer variability:\n")
print(strain_summary.head(10))

# ============ STAGE 2: REGRESS ON AGE & SEX ============
print("\n" + "="*80)
print("STAGE 2: Regression Analysis (HI ~ Age + Sex)")
print("="*80)

results = []

for strain in strain_summary.index:
    df = analysis_df[analysis_df["Target Entity Subtype"] == strain].copy()
    
    # Clean sex variable
    df["sex"] = df["Biological sex"].str.lower().str.strip()
    df = df[df["sex"].isin(["female", "male"])]
    df = df.dropna(subset=["Value_log2", "Age", "sex"])
    
    n_obs = len(df)
    n_part = df["Participant ID"].nunique()
    
    if n_obs < 10:
        print(f"\n{strain[:50]}")
        print(f"  ✗ SKIP: only {n_obs} obs")
        continue
    
    print(f"\n{strain[:60]}")
    print(f"  N = {n_obs} obs, {n_part} participants")
    
    # Regression
    import statsmodels.api as sm
    
    df["sex_numeric"] = (df["sex"] == "male").astype(int)
    X = df[["Age", "sex_numeric"]].copy()
    X = sm.add_constant(X)
    y = df["Value_log2"].copy()
    
    try:
        model = sm.OLS(y, X).fit()
        
        r2 = model.rsquared
        residuals = model.resid
        residual_var = np.var(residuals, ddof=2)
        total_var = np.var(y, ddof=0)
        
        age_coef = model.params["Age"]
        age_pval = model.pvalues["Age"]
        sex_coef = model.params["sex_numeric"]
        sex_pval = model.pvalues["sex_numeric"]
        
        var_explained_pct = 100 * r2
        var_residual_pct = 100 * (1 - r2)
        
        results.append({
            "Strain": strain,
            "N_Obs": n_obs,
            "N_Participants": n_part,
            "N_Studies": df["Study ID"].nunique(),
            "Median_HI_log2": df["Value_log2"].median(),
            "Std_HI_log2": df["Value_log2"].std(),
            "Total_Variance": total_var,
            "Age_Coefficient": age_coef,
            "Age_Pvalue": age_pval,
            "Sex_Coefficient": sex_coef,
            "Sex_Pvalue": sex_pval,
            "R_squared": r2,
            "Var_Explained_pct": var_explained_pct,
            "Residual_Variance_pct": var_residual_pct,
            "Residual_Variance_absolute": residual_var,
        })
        
        print(f"  ✓ Age: β={age_coef:.4f} (p={age_pval:.4f})")
        print(f"  ✓ Sex: β={sex_coef:.4f} (p={sex_pval:.4f})")
        print(f"  ✓ R² = {r2:.3f} → Explained: {var_explained_pct:.1f}% | Residual: {var_residual_pct:.1f}%")
        
    except Exception as e:
        print(f"  ✗ ERROR: {type(e).__name__}: {str(e)[:80]}")

# ============ RESULTS ============
print("\n" + "="*80)
print("RESULTS: STRAINS RANKED BY RESIDUAL VARIANCE")
print("="*80 + "\n")

results_df = pd.DataFrame(results)

if len(results_df) > 0:
    results_df = results_df.sort_values("Residual_Variance_pct", ascending=False)
    
    print(results_df[[
        "Strain", "N_Participants", "Std_HI_log2", 
        "Var_Explained_pct", "Residual_Variance_pct"
    ]].to_string(index=False))
    
    # Save
    results_df.to_csv("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hi_titer_screening_results.csv", index=False)
    print(f"\n✓ Results saved to 'hi_titer_screening_results.csv'")
    
    # ============ PLOTTING ============
    print("\n" + "="*80)
    print("GENERATING PLOTS")
    print("="*80 + "\n")
    
   
    
    # Plot 1: Main result - variance decomposition
    fig, ax = plt.subplots(figsize=(14, 8))
    
    strains_short = [s.replace("Influenza ", "").replace("virus (", "\n(")[:50] for s in results_df["Strain"]]
    x = np.arange(len(results_df))
    
    explained = results_df["Var_Explained_pct"]
    residual = results_df["Residual_Variance_pct"]
    
    ax.bar(x, explained, width=0.6, label="Explained by Age & Sex", color="#4472C4", alpha=0.9)
    ax.bar(x, residual, width=0.6, bottom=explained, label="Residual (potential genetic)", 
           color="#ED7D31", alpha=0.85)
    
    ax.set_xlabel("Strain", fontsize=12, fontweight="bold")
    ax.set_ylabel("Variance (%)", fontsize=12, fontweight="bold")
    ax.set_title("HI Titer Variance Decomposition\n(Orange = unexplained, potential genetic signal)", 
                 fontsize=13, fontweight="bold", pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(strains_short, rotation=45, ha="right", fontsize=9)
    ax.set_ylim([0, 105])
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/01_variance_decomposition.png", dpi=300, bbox_inches="tight")
    print("✓ Saved: 01_variance_decomposition.png")
    plt.close()
    
    # Plot 2: Residual variance ranking
    fig, ax = plt.subplots(figsize=(12, 7))
    
    top_n = min(15, len(results_df))
    top_strains = results_df.head(top_n).sort_values("Residual_Variance_pct")
    y_pos = np.arange(len(top_strains))
    
    colors = ["#d73027" if x > 60 else "#fee090" if x > 40 else "#91bfdb" 
              for x in top_strains["Residual_Variance_pct"]]
    
    bars = ax.barh(y_pos, top_strains["Residual_Variance_pct"], color=colors, 
                    edgecolor="black", linewidth=0.8)
    
    for i, (idx, row) in enumerate(top_strains.iterrows()):
        ax.text(row["Residual_Variance_pct"] + 1, i, f"{row['Residual_Variance_pct']:.1f}%", 
                va="center", fontsize=9, fontweight="bold")
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels([s.replace("Influenza ", "").replace("virus (", "\n(")[:45] for s in top_strains["Strain"]], 
                        fontsize=9)
    ax.set_xlabel("Residual Variance (%)", fontsize=11, fontweight="bold")
    ax.set_title("Top Strains by Residual Variance\nRed=high (>60%), Yellow=medium (40-60%), Blue=low (<40%)",
                 fontsize=12, fontweight="bold", pad=15)
    ax.set_xlim([0, 100])
    ax.axvline(x=60, color="red", linestyle="--", linewidth=1.5, alpha=0.5)
    ax.axvline(x=40, color="orange", linestyle="--", linewidth=1.5, alpha=0.5)
    ax.grid(axis="x", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/02_residual_variance_ranking.png", dpi=300, bbox_inches="tight")
    print("✓ Saved: 02_residual_variance_ranking.png")
    plt.close()
    
    # Plot 3: Scatter - residual variance vs. strain variability
    fig, ax = plt.subplots(figsize=(10, 6))
    
    scatter = ax.scatter(results_df["Std_HI_log2"], results_df["Residual_Variance_pct"], 
                         s=results_df["N_Participants"]*2, alpha=0.6, 
                         c=results_df["Residual_Variance_pct"], cmap="YlOrRd", 
                         edgecolors="black", linewidth=0.5)
    
    ax.set_xlabel("HI Titer Variability (Std of log₂ titers)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Residual Variance (%)", fontsize=11, fontweight="bold")
    ax.set_title("Residual Variance vs. HI Titer Variability\n(Bubble size = sample size)", 
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Residual Variance (%)", fontsize=10)
    
    plt.tight_layout()
    plt.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/03_residual_vs_variability.png", dpi=300, bbox_inches="tight")
    print("✓ Saved: 03_residual_vs_variability.png")
    plt.close()
    
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print(f"""
Generated 3 plots:
  1. 01_variance_decomposition.png - Main result
  2. 02_residual_variance_ranking.png - Top strains for follow-up
  3. 03_residual_vs_variability.png - Scatter plot relationship

Interpretation:
  • Strains with HIGH residual variance (>60%, red) are good candidates
    for genetic studies (more variation unexplained by age/sex)
  • Strains with LOW residual variance (<40%, blue) have most variation
    explained by demographics
    """)
else:
    print("✗ No strains analyzed successfully.")





















