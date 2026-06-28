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

demo = pd.read_csv(os.path.join(SCRATCH,'datatools_demographic_Table.csv'), header = 0)
hai = pd.read_csv(os.path.join(SCRATCH,'datatools_HAI_Table.csv'), header = 0)


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
analysis_df["log2_HAI"] = np.log2(analysis_df["Value"].replace(0, np.nan))

analysis_df = analysis_df[
    analysis_df["Biological sex"].isin(["female", "male"])
]










analysis_df['log2_HAI'].describe()

plt.hist(analysis_df['log2_HAI'], bins=10, alpha=0.7)
plt.xlabel("Value (log2)")
plt.ylabel("Count")
plt.show()



#strain vs HI
analysis_df.boxplot(
    column="log2_HAI",
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
    analysis_df["log2_HAI"],
    alpha=0.3
)
plt.xlabel("Age")
plt.ylabel("HI titer")
plt.show()

#sex vs HI
analysis_df.boxplot(
    column="log2_HAI",
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









##### NEW DATA 

merge11 = pd.merge(
    demo,
    hai,
    how='outer',
    on= ['Cohort', 'Participant ID'],
    suffixes=('_demo', '_hai')
       
)

merge1= merge11.copy()
merge1["log2_HAI"] = np.log2(merge1["Value Preferred"].replace(0, np.nan))

#merge1['Age Reported_demo'].isna().sum()
#merge1['Gender_demo'].isna().sum()
#merge1['Race_demo'].isna().sum()


#day 0, 3, 14, 28, 180

day0 =merge1[merge1['Study Time Collected'] == 0]
day3 =merge1[merge1['Study Time Collected'] == 3]
day14 =merge1[merge1['Study Time Collected'] == 14]
day28 =merge1[merge1['Study Time Collected'] == 28]
day180 =merge1[merge1['Study Time Collected'] == 180]




#start with day 0 and day 28 


analysis_df['log2_HAI'].describe()

plt.hist(analysis_df['Value_log2'], bins=10, alpha=0.7)
plt.xlabel("Value (log2)")
plt.ylabel("Count")
plt.show()



#strain vs HI
merge1.boxplot(
    column="log2_HAI",
    by="Virus",
    rot=90,
    figsize=(12, 6)
)
plt.title("Day 0 HAI titers by strain")
plt.xlabel("Influenza strain")
plt.ylabel("log2 HAI")
plt.show()



#cohort vs HI
merge1.boxplot(
    column="log2_HAI",
    by="Cohort",
    rot=90,
    figsize=(12, 6)
)
plt.title("Day 0 HAI titers by cohort")
plt.xlabel("Cohort")
plt.ylabel("log2 HAI")
plt.show()




import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import warnings

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------
# 1. LOAD DATA -- point this at your real file
# ---------------------------------------------------------------
df = day0.copy()

# ---------------------------------------------------------------
# 2. COLUMN MAP -- adjust right-hand values if your headers differ
# ---------------------------------------------------------------
COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_demo",
    "sex": "Gender_demo",
    "cohort": "Cohort",
    "virus": "Virus",
    "participant_id": "Participant ID",
}

clean = df.rename(columns={
    COLUMNS["outcome"]: "log2_HAI",
    COLUMNS["age"]: "age",
    COLUMNS["sex"]: "sex",
    COLUMNS["cohort"]: "cohort",
    COLUMNS["virus"]: "virus",
    COLUMNS["participant_id"]: "participant_id",
})

# ---------------------------------------------------------------
# 3. BASIC CLEANING
# ---------------------------------------------------------------
needed = ["log2_HAI", "age", "sex", "cohort", "virus"]
clean = clean.dropna(subset=needed).copy()

clean["sex"] = clean["sex"].astype("category")
clean["cohort"] = clean["cohort"].astype("category")
clean["age"] = pd.to_numeric(clean["age"], errors="coerce")
clean = clean.dropna(subset=["age"])

# ---------------------------------------------------------------
# 4. FIT MIXED MODEL SEPARATELY PER VIRUS STRAIN
#    log2_HAI ~ age + sex + (1 | cohort)
# ---------------------------------------------------------------
results = {}

for virus_name, sub in clean.groupby("virus"):
    sub = sub.copy()
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)

    print("=" * 70)
    print(f"Virus: {virus_name}   (n_obs={n_obs}, n_cohorts={n_cohorts})")
    print("=" * 70)

    if n_cohorts < 2:
        print("  -> Skipped: need >=2 cohorts to estimate a random effect.\n")
        continue
    if n_obs < 10:
        print("  -> Warning: very small sample size for this virus; results unstable.\n")

    model = smf.mixedlm(
        "log2_HAI ~ age + sex",
        data=sub,
        groups=sub["cohort"],
        re_formula="1",   # random intercept only
    )

    try:
        fit = model.fit(reml=True)
        results[virus_name] = fit
        print(fit.summary())
    except Exception as e:
        print(f"  -> Model failed to converge: {e}")

    print()

# ---------------------------------------------------------------
# 5. COMPACT SUMMARY TABLE ACROSS STRAINS
# ---------------------------------------------------------------
rows = []
for virus_name, fit in results.items():
    params = fit.params
    pvals = fit.pvalues
    ci = fit.conf_int()

    for term in params.index:
        if term in ("Group Var",):
            continue
        rows.append({
            "Virus": virus_name,
            "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4),
            "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Virus": virus_name,
        "Term": "Cohort (random intercept) variance",
        "Estimate": round(fit.cov_re.iloc[0, 0], 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })
    rows.append({
        "Virus": virus_name,
        "Term": "Residual variance",
        "Estimate": round(fit.scale, 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })

summary_df = pd.DataFrame(rows)
print("=" * 70)
print("SUMMARY TABLE (all strains)")
print("=" * 70)
print(summary_df.to_string(index=False))

#summary_df.to_csv("hai_mixedlm_summary.csv", index=False)






"""
Model diagnostics for the HAI mixed-effects regression
=======================================================
For each virus strain's fitted MixedLM model, produces three standard
diagnostic plots, arranged in a grid (rows = virus strains):

  1. Residuals vs Fitted        -> checks linearity / non-random patterns
  2. Scale-Location             -> checks homoscedasticity (constant variance);
                                    this is the "residuals vs variance" plot
  3. Normal Q-Q                 -> checks normality of residuals

HOW TO USE WITH YOUR REAL DATA
-------------------------------
This script re-fits the same models as hai_mixedlm.py, so update the same
two things here:
  - the file path in step 1
  - the COLUMNS dict in step 2 (use check_columns.py first if unsure)
"""

import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import scipy.stats as stats
import warnings

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------
# 1. LOAD DATA -- point this at your real file
# ---------------------------------------------------------------
df = day0.copy()

# ---------------------------------------------------------------
# 2. COLUMN MAP -- adjust right-hand values if your headers differ
# ---------------------------------------------------------------
df.columns = df.columns.str.strip()  # guard against stray whitespace in headers

COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_demo",
    "sex": "Gender_demo",
    "cohort": "Cohort",
    "virus": "Virus",
}

clean = df.rename(columns={
    COLUMNS["outcome"]: "log2_HAI",
    COLUMNS["age"]: "age",
    COLUMNS["sex"]: "sex",
    COLUMNS["cohort"]: "cohort",
    COLUMNS["virus"]: "virus",
})

missing = [c for c in ["log2_HAI", "age", "sex", "cohort", "virus"] if c not in clean.columns]
if missing:
    raise KeyError(f"Rename didn't produce expected columns: {missing}. "
                    f"Available columns are: {df.columns.tolist()}")

clean = clean.dropna(subset=["log2_HAI", "age", "sex", "cohort", "virus"]).copy()
clean["sex"] = clean["sex"].astype("category")
clean["cohort"] = clean["cohort"].astype("category")
clean["age"] = pd.to_numeric(clean["age"], errors="coerce")
clean = clean.dropna(subset=["age"])

# ---------------------------------------------------------------
# 3. FIT MODEL PER VIRUS STRAIN, COLLECT FITTED VALUES + RESIDUALS
# ---------------------------------------------------------------
fitted_models = {}

for virus_name, sub in clean.groupby("virus"):
    sub = sub.copy()
    if sub["cohort"].nunique() < 2:
        print(f"Skipping {virus_name}: needs >=2 cohorts for random effect.")
        continue

    model = smf.mixedlm("log2_HAI ~ age + sex", data=sub, groups=sub["cohort"], re_formula="1")
    try:
        fit = model.fit(reml=True)
    except Exception as e:
        print(f"Skipping {virus_name}: model failed to converge ({e}).")
        continue

    fitted_models[virus_name] = {
        "fitted": fit.fittedvalues,
        "resid": fit.resid,
        "fit_obj": fit,
    }

if not fitted_models:
    raise RuntimeError("No models converged — nothing to plot.")

# ---------------------------------------------------------------
# 4. DIAGNOSTIC PLOTS: one row per virus strain, 3 columns
# ---------------------------------------------------------------
n_models = len(fitted_models)
fig, axes = plt.subplots(n_models, 3, figsize=(15, 4.5 * n_models), squeeze=False)

for row_idx, (virus_name, m) in enumerate(fitted_models.items()):
    fitted_vals = m["fitted"]
    resid = m["resid"]

    # Standardized residuals (for scale-location and Q-Q)
    resid_std = (resid - resid.mean()) / resid.std()
    sqrt_abs_std_resid = np.sqrt(np.abs(resid_std))

    # --- (a) Residuals vs Fitted ---
    ax = axes[row_idx, 0]
    ax.scatter(fitted_vals, resid, alpha=0.7, edgecolor="k", linewidth=0.3)
    ax.axhline(0, color="red", linestyle="--", linewidth=1)
    ax.set_xlabel("Fitted values")
    ax.set_ylabel("Residuals")
    ax.set_title(f"{virus_name}\nResiduals vs Fitted")

    # --- (b) Scale-Location (sqrt|standardized residuals| vs fitted) ---
    ax = axes[row_idx, 1]
    ax.scatter(fitted_vals, sqrt_abs_std_resid, alpha=0.7, edgecolor="k", linewidth=0.3)
    ax.set_xlabel("Fitted values")
    ax.set_ylabel(r"$\sqrt{|\mathrm{Standardized\ residuals}|}$")
    ax.set_title("Scale-Location\n(checks variance is constant)")

    # --- (c) Normal Q-Q plot ---
    ax = axes[row_idx, 2]
    stats.probplot(resid_std, dist="norm", plot=ax)
    ax.set_title("Normal Q-Q")
    ax.get_lines()[0].set_markerfacecolor("steelblue")
    ax.get_lines()[0].set_markeredgecolor("k")
    ax.get_lines()[0].set_alpha(0.7)
    ax.get_lines()[1].set_color("red")

fig.tight_layout()
fig.savefig("/home/claude/hai_model_diagnostics.png", dpi=150, bbox_inches="tight")
print("Saved: /home/claude/hai_model_diagnostics.png")

# ---------------------------------------------------------------
# 5. QUICK NUMERIC CHECKS TO ACCOMPANY THE PLOTS
# ---------------------------------------------------------------
print("\n" + "=" * 70)
print("NUMERIC DIAGNOSTIC SUMMARY")
print("=" * 70)
for virus_name, m in fitted_models.items():
    resid = m["resid"]
    shapiro_stat, shapiro_p = stats.shapiro(resid)
    print(f"\n{virus_name}")
    print(f"  Residual mean:      {resid.mean():.4f}  (should be ~0)")
    print(f"  Residual std:       {resid.std():.4f}")
    print(f"  Shapiro-Wilk p:     {shapiro_p:.4f}  "
          f"({'residuals look normal' if shapiro_p > 0.05 else 'deviation from normality detected'})")
