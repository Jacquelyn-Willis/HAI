
!{sys.executable} -m pip install ipykernel --upgrade --force-reinstall
!{sys.executable} -m pip install statsmodels
!{sys.executable} -m pip install matplotlib
!{sys.executable} -m pip install patsy

import sys
import pandas as pd 
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from patsy.contrasts import Treatment
import warnings
warnings.filterwarnings("ignore")
import scipy.stats as stats
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#directories 

SCRATCH = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/"

data = "/sc/arion/work/willij115/projects/HAI/data/2026-05-20_HAI_covariate_regression"
scratch = "/sc/arion/scratch/willij115/projects/HAI/2026-05-20_HAI_covariate_regression"
results = "/sc/arion/work/willij115/projects/HAI/results/2026-05-20_HAI_covariate_regression"

pd.set_option('display.max_columns', None)
pd.set_option('display.max_r', None)



#upload immunespace data tables
 
studies = pd.read_csv(os.path.join(SCRATCH, 'immunespaceHAI_studies_tables.csv'), names=['Study ID', 'Study Title', 'Study Data Release', 'PMID', 'Publication Title', 'Publication Date', 'Author Count', 'Authors'], header = 0)
arms = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_arms_tables.csv'),  header = 0)
participants = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_participants_tables.csv'), header = 0)
events = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_events_tables.csv'), header = 0)
assays = pd.read_csv(os.path.join(SCRATCH,'immunespaceHAI_assays_tables.csv'), header = 0)

demo = pd.read_csv(os.path.join(SCRATCH,'datatools_demographic_Table.csv'), header = 0)
hai = pd.read_csv(os.path.join(SCRATCH,'datatools_HAI_Table.csv'), header = 0)



######## need to merge: merge1 , arms, and oarticipants to get study ID and cohort decriptions 

merge1 = pd.merge(
    demo,
    hai,
    how='inner',
    on= ['Cohort', 'Participant ID'],
    suffixes=('_demo', '_hai')
       
)

merge1_log2= merge1.copy()
merge1_log2["log2_HAI"] = np.log2(merge1["Value Preferred"].replace(0, np.nan))
merge1_log2["Study_ID"] = merge1_log2["Participant ID"].str.split(".").str[1]


participants["Study_ID"] = participants["Study ID"].str.split("Y").str[1]
arms["Study_ID"] = arms["Study ID"].str.split("Y").str[1]

merge2 = pd.merge(
    arms,
    participants[['Participant ID', 'Arm ID', 'Study ID', 'Study_ID']],
    how='right',
    on='Study_ID',
    suffixes=('_arms', '_par')
)

merge2["new_participant_id"] = (
     merge2["Participant ID"].astype(str).str.split(".").str[-1]
    + "."
    + merge2["Study_ID"].astype(str)
)



final_merge = pd.merge(
    merge2,
    merge1_log2,
    how='right',
    left_on=['new_participant_id', 'Study_ID'], 
    right_on=['Participant ID', 'Study_ID'],
    suffixes=('_merge2', '_merge1')    
)



strings_to_remove = [
    "placebo",
    "saline",
    "type 2 diabetes",
    "Pneunomax23",
    "young T2D",
    "old T2D"
]
pattern = "|".join(map(str, strings_to_remove))

mask = (
    final_merge["Cohort"].str.contains(pattern, case=False, na=False) |
    final_merge["Description_merge2"].str.contains(pattern, case=False, na=False)
)

final_merge = final_merge[~mask]

'''
study_cohort_desc = (


    final_merge[
        ["Study ID_par", "Cohort", "Description_merge2"]
    ]
    .drop_duplicates()
    .sort_values(["Study ID_par", "Cohort"])
)

study_cohort_desc2 = (


    final_merge2[
        ["Study ID_par", "Cohort", "Description_merge2"]
    ]
    .drop_duplicates()
    .sort_values(["Study ID_par", "Cohort"])
)


'''










STRAIN_COL = "Virus"
COHORT_COL = "Cohort"
MIN_N = 5

os.makedirs(os.path.join(SCRATCH, "hai_demo_plots_html"), exist_ok=True)

def plot_group_plotly(df, label, outdir= os.path.join(SCRATCH, "hai_demo_plots_html")):
    print(f"\n=== {label} (n={len(df)}) ===")
    print(df["log2_HAI"].describe())

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("HAI distribution", "Age vs HAI", "HAI by sex", ""),
        vertical_spacing=0.12,
        horizontal_spacing=0.10,
    )

    # Histogram
    fig.add_trace(
        go.Histogram(x=df["log2_HAI"], nbinsx=10, name="log2_HAI", opacity=0.7),
        row=1, col=1
    )

    # Scatter: age vs HAI
    fig.add_trace(
        go.Scatter(
            x=df["Age Reported_demo"],
            y=df["log2_HAI"],
            mode="markers",
            marker=dict(opacity=0.3),
            name="Age vs HAI",
        ),
        row=1, col=2
    )

    # Box: sex vs HAI
    if df["Gender_demo"].nunique() > 1:
        for gender, gdf in df.groupby("Gender_demo"):
            fig.add_trace(
                go.Box(
                    y=gdf["log2_HAI"],
                    name=str(gender),
                    boxmean=True,
                ),
                row=2, col=1
            )
    else:
        # leave empty if only one sex category
        pass

    fig.update_xaxes(title_text="Value (log2)", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)

    fig.update_xaxes(title_text="Age", row=1, col=2)
    fig.update_yaxes(title_text="log2 HAI", row=1, col=2)

    fig.update_xaxes(title_text="gender", row=2, col=1)
    fig.update_yaxes(title_text="log2 HAI", row=2, col=1)

    fig.update_layout(
        title=label,
        width=1200,
        height=900,
        showlegend=False,
        template="plotly_white",
    )

    filename = label.replace(" ", "_").replace("/", "_").replace("—", "_") + ".html"
    filepath = os.path.join(outdir, filename)
    fig.write_html(filepath, include_plotlyjs="cdn", full_html=True)
    return filepath


html_files = []

for (strain, cohort), sub_df in final_merge.groupby([STRAIN_COL, COHORT_COL]):
    if len(sub_df) < MIN_N:
        print(f"Skipping {strain} / {cohort} (n={len(sub_df)} < {MIN_N})")
        continue

    path = plot_group_plotly(sub_df, f"{strain} — {cohort}")
    html_files.append(path)

print("Saved HTML files:")
for f in html_files:
    print(f)







"""
HAI regression: auto LMM (random effect) or OLS, decided per strain
=====================================================================
For EACH virus strain, the script automatically decides which model to fit:

  - >=2 cohorts present for that strain:
        LMM:  log2_HAI ~ age + sex + (1 | cohort)
        Cohort is a random intercept -- appropriate because there's more
        than one cohort to estimate variance across.

  - Only 1 cohort present for that strain:
        OLS:  log2_HAI ~ age + sex
        Cohort can't be modeled (no variation across cohorts to estimate),
        so it's dropped entirely and only age + sex are covariates.

No model is fit comparing LMM vs OLS on the same strain -- each strain gets
exactly one model, chosen by its own cohort count.

HOW TO USE WITH YOUR REAL DATA
-------------------------------
Replace the line that reads "example_hai_data.csv" with your real file, e.g.:

    df = pd.read_csv("/mnt/user-data/uploads/your_file.csv")
    # or: df = pd.read_excel("/mnt/user-data/uploads/your_file.xlsx")

If your real headers differ from the ones below, run check_columns.py first
and update the COLUMNS dict accordingly.
"""





# ---------------------------------------------------------------
# 1. LOAD DATA -- swap this line for your real file
# ---------------------------------------------------------------
df = final_merge.copy()  # or day3, day14, day28, day180, or merge1 for all days
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

missing = [c for c in ["log2_HAI", "age", "sex", "cohort", "virus"] if c not in clean.columns]
if missing:
    raise KeyError(f"Rename didn't produce expected columns: {missing}. "
                    f"Available columns are: {df.columns.tolist()}")

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
# 3b. MARGINAL / CONDITIONAL R^2 FOR MIXED MODELS
#     (Nakagawa & Schielzeth 2013 method, theoretical/delta approach
#     for a linear mixed model with a single random intercept)
#
#     Marginal R^2   = Var(fixed)              / [Var(fixed) + Var(random) + Var(resid)]
#     Conditional R^2 = [Var(fixed) + Var(random)] / [Var(fixed) + Var(random) + Var(resid)]
#
#     Var(fixed)  = variance of fitted values using ONLY the fixed-effect
#                   part of the model (X @ beta), i.e. excluding random
#                   intercepts entirely
#     Var(random) = the random intercept variance estimated by the model
#                   (fit.cov_re), i.e. "Group Var" / "Cohort variance"
#     Var(resid)  = residual variance (fit.scale)
# ---------------------------------------------------------------
def mixedlm_r2(fit):
    """Returns (marginal_r2, conditional_r2) for a fitted MixedLM (random intercept only)."""
    fixed_effects_fitted = fit.model.exog @ fit.fe_params  # X @ beta, fixed part only
    var_fixed = np.var(fixed_effects_fitted, ddof=0)
    var_random = fit.cov_re.iloc[0, 0]
    var_resid = fit.scale

    total = var_fixed + var_random + var_resid
    marginal_r2 = var_fixed / total
    conditional_r2 = (var_fixed + var_random) / total
    return marginal_r2, conditional_r2

# ---------------------------------------------------------------
# 4. PER-STRAIN MODEL SELECTION + FIT
#    Decision rule:
#      n_cohorts(strain) >= 2  -> LMM with random intercept for cohort
#      n_cohorts(strain) == 1  -> OLS, age + sex only (no cohort term)
#      n_cohorts(strain) == 0 / no data -> skipped
# ---------------------------------------------------------------
lmm_results = {}
ols_results = {}
model_choice = {}  # virus -> "LMM" or "OLS", for the summary table
r2_values = {}     # virus -> (marginal_r2, conditional_r2)

for virus_name, sub in clean.groupby("virus"):
    sub = sub.copy()
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)

    print("=" * 70)
    print(f"Virus: {virus_name}   (n_obs={n_obs}, n_cohorts={n_cohorts})")
    print("=" * 70)

    if n_obs < 5:
        print("  -> Skipped: too few observations to fit anything meaningful.\n")
        continue

    if n_cohorts >= 2:
        # ---- LMM: random intercept for cohort ----
        model_choice[virus_name] = "LMM"
        print("  -> Decision: >=2 cohorts present -> fitting LMM (1 | cohort)")
        model = smf.mixedlm(
            "log2_HAI ~ age + sex",
            data=sub,
            groups=sub["cohort"],
            re_formula="1",
        )
        try:
            fit = model.fit(reml=True)
            lmm_results[virus_name] = fit
            print("\n--- LMM: log2_HAI ~ age + sex + (1 | cohort) ---")
            print(fit.summary())

            marg_r2, cond_r2 = mixedlm_r2(fit)
            r2_values[virus_name] = (marg_r2, cond_r2)
            print(f"\n  Marginal R^2 (fixed effects only):        {marg_r2:.4f}")
            print(f"  Conditional R^2 (fixed + random effects): {cond_r2:.4f}")
        except Exception as e:
            print(f"  -> LMM failed to converge: {e}")

    else:
        # ---- OLS: age + sex only, no cohort term (only 1 cohort present) ----
        model_choice[virus_name] = "OLS"
        print("  -> Decision: only 1 cohort present -> fitting OLS (age + sex only)")
        try:
            ols_fit = smf.ols("log2_HAI ~ age + sex", data=sub).fit()
            ols_results[virus_name] = ols_fit
            print("\n--- OLS: log2_HAI ~ age + sex ---")
            print(ols_fit.summary())

            # No random effect exists for OLS, so marginal R^2 == conditional R^2
            # == standard R^2 (fixed effects explain all the explained variance)
            r2_values[virus_name] = (ols_fit.rsquared, ols_fit.rsquared)
            print(f"\n  R^2 (= marginal = conditional, no random effect modeled): {ols_fit.rsquared:.4f}")
        except Exception as e:
            print(f"  -> OLS failed: {e}")

    print()

# ---------------------------------------------------------------
# 5. UNIFIED SUMMARY TABLE ACROSS ALL STRAINS
#    (one row per term, regardless of which model type was used)
# ---------------------------------------------------------------
rows = []

for virus_name, fit in lmm_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()
    for term in params.index:
        if term == "Group Var":
            continue
        rows.append({
            "Virus": virus_name, "Model": "LMM", "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4),
            "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Virus": virus_name, "Model": "LMM", "Term": "Cohort (random intercept) variance",
        "Estimate": round(fit.cov_re.iloc[0, 0], 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })
    rows.append({
        "Virus": virus_name, "Model": "LMM", "Term": "Residual variance",
        "Estimate": round(fit.scale, 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })

for virus_name, fit in ols_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()
    for term in params.index:
        rows.append({
            "Virus": virus_name, "Model": "OLS", "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4),
            "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Virus": virus_name, "Model": "OLS", "Term": "Residual variance",
        "Estimate": round(fit.mse_resid, 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })
    # No "Cohort variance" row for OLS strains -- cohort wasn't modeled at all

summary_df = pd.DataFrame(rows)
print("=" * 70)
print("SUMMARY TABLE (all strains, model chosen automatically per strain)")
print("=" * 70)
print(summary_df.to_string(index=False))

#summary_df.to_csv("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_regression_summary.csv", index=False)
print("\nSaved: /Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_regression_summary.csv")

# ---------------------------------------------------------------
# 6. MODEL CHOICE LOG -- which strains got which model, and why
# ---------------------------------------------------------------
choice_rows = []
for virus_name, sub in clean.groupby("virus"):
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)
    marg_r2, cond_r2 = r2_values.get(virus_name, (np.nan, np.nan))
    choice_rows.append({
        "Virus": virus_name,
        "N_Obs": n_obs,
        "N_Cohorts": n_cohorts,
        "Model_Used": model_choice.get(virus_name, "Skipped (n_obs<5)"),
        "Reason": ("Cohort modeled as random intercept (>=2 cohorts)" if n_cohorts >= 2
                   else "Cohort not modeled -- only 1 cohort present" if n_cohorts == 1
                   else "No cohort data"),
        "Marginal_R2": round(marg_r2, 4) if pd.notna(marg_r2) else np.nan,
        "Conditional_R2": round(cond_r2, 4) if pd.notna(cond_r2) else np.nan,
    })

choice_df = pd.DataFrame(choice_rows)
#choice_df.to_csv("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_model_choice_log.csv", index=False)
print("\n" + "=" * 70)
print("MODEL CHOICE LOG")
print("=" * 70)
print(choice_df.to_string(index=False))
print("\nSaved: /Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_model_choice_log.csv")






# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
df = pd.read_csv("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_model_choice_log.csv")
df = df.sort_values("Conditional_R2", ascending=False).reset_index(drop=True)

def short_label(name, maxlen=30):
    return str(name).replace("Influenza ", "")[:maxlen]

df["Strain_short"] = df["Virus"].apply(short_label)

# ---------------------------------------------------------------
# PLOT 1: Grouped bar chart -- Marginal vs Conditional R^2 per strain
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(max(8, 1.8 * len(df)), 7))

x = np.arange(len(df))
width = 0.35

bars1 = ax.bar(x - width/2, df["Marginal_R2"], width, label="Marginal R² (fixed effects only)",
               color="#4472C4", edgecolor="black", linewidth=0.6)
bars2 = ax.bar(x + width/2, df["Conditional_R2"], width, label="Conditional R² (fixed + random effects)",
               color="#ED7D31", edgecolor="black", linewidth=0.6)

for b in bars1:
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.02, f"{b.get_height():.3f}",
            ha="center", va="bottom", fontsize=9)
for b in bars2:
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.02, f"{b.get_height():.3f}",
            ha="center", va="bottom", fontsize=9)

# annotate model type used per strain
for i, model_used in enumerate(df["Model_Used"]):
    ax.text(i, -0.06, f"[{model_used}]", ha="center", va="top", fontsize=8,
            color="gray", style="italic")

ax.set_xlabel("Strain", fontsize=12, fontweight="bold")
ax.set_ylabel("R²", fontsize=12, fontweight="bold")
ax.set_title("Marginal vs Conditional R² by Strain", fontsize=13, fontweight="bold", pad=15)
ax.set_xticks(x)
ax.set_xticklabels(df["Strain_short"], rotation=45, ha="right", fontsize=9)
ax.set_ylim([0, 1.0])
ax.legend(loc="upper right", fontsize=10)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
fig.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/r2_marginal_vs_conditional.png", dpi=300, bbox_inches="tight")
print("Saved: r2_marginal_vs_conditional.png")
plt.close(fig)

# ---------------------------------------------------------------
# PLOT 2: Gap plot -- how much cohort (random effect) adds, per strain
# ---------------------------------------------------------------
df["R2_gap"] = df["Conditional_R2"] - df["Marginal_R2"]
df_gap_sorted = df.sort_values("R2_gap", ascending=True)

fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(df))))

y_pos = np.arange(len(df_gap_sorted))
colors = ["#91bfdb" if g == 0 else "#4472C4" for g in df_gap_sorted["R2_gap"]]

bars = ax.barh(y_pos, df_gap_sorted["R2_gap"], color=colors, edgecolor="black", linewidth=0.7)

for i, (_, row) in enumerate(df_gap_sorted.iterrows()):
    label = f"{row['R2_gap']:.3f}" if row["R2_gap"] > 0 else "0 (OLS, no random effect)"
    ax.text(row["R2_gap"] + 0.01, i, label, va="center", fontsize=9)

ax.set_yticks(y_pos)
ax.set_yticklabels(df_gap_sorted["Strain_short"], fontsize=9)
ax.set_xlabel("Conditional R² − Marginal R²\n(variance explained by cohort random effect)",
              fontsize=11, fontweight="bold")
ax.set_title("How Much Does Cohort (Random Effect) Add to R²?", fontsize=12, fontweight="bold", pad=15)
ax.set_xlim([0, max(0.1, df["R2_gap"].max() * 1.3)])
ax.grid(axis="x", alpha=0.3)

plt.tight_layout()
fig.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/r2_cohort_contribution_gap.png", dpi=300, bbox_inches="tight")
print("Saved: r2_cohort_contribution_gap.png")
plt.close(fig)

print("\nDone. Strains with R2_gap = 0 used OLS (single cohort, no random effect to add).")














"""
HAI variance decomposition plots
=================================
Builds the variance-decomposition table directly from `summary_df`
(long format: one row per Virus x Term) and makes three plots:

  1. Stacked bar: Cohort vs Residual variance share, per strain (%)
  2. Ranked horizontal bar: strains sorted by Residual Variance %
  3. Scatter: Residual Variance % vs raw titer variability (Std of log2_HAI),
     bubble size = N participants

NOTE ON LABELING
-----------------
The "residual" share is labeled "Unexplained / individual-level variance",
NOT "potential genetic signal". The model only adjusts for age, sex, and
cohort -- whatever variance is left over could reflect genetics, but could
equally reflect assay noise, batch effects, prior exposure history, timing
of sample collection, etc. The model has no way to distinguish these, so
the plot does not claim a genetic interpretation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

STRAIN_COL = "Virus"   # column in analysis_df matching summary_df's "Virus"

# ---------------------------------------------------------------
# 1. RESHAPE summary_df (long) -> one row per strain with variance %s
# ---------------------------------------------------------------
var_rows = summary_df[summary_df["Term"].isin(
    ["Cohort (random intercept) variance", "Residual variance"]
)].copy()

wide_var = var_rows.pivot(index="Virus", columns="Term", values="Estimate").reset_index()
wide_var = wide_var.rename(columns={
    "Cohort (random intercept) variance": "Cohort_Variance",
    "Residual variance": "Residual_Variance",
    "Virus": "Strain"
})

wide_var["Total_Variance"] = wide_var["Cohort_Variance"] + wide_var["Residual_Variance"]
wide_var["Cohort_Variance_pct"] = 100 * wide_var["Cohort_Variance"] / wide_var["Total_Variance"]
wide_var["Residual_Variance_pct"] = 100 * wide_var["Residual_Variance"] / wide_var["Total_Variance"]

# ---------------------------------------------------------------
# 2. PULL Std_HAI_log2 and N_Participants from analysis_df, per strain
# ---------------------------------------------------------------
stats_by_strain = (
    day0.groupby(STRAIN_COL)["log2_HAI"]
    .agg(Std_HAI_log2="std", N_Participants="count")
    .reset_index()
    .rename(columns={STRAIN_COL: "Strain"})
)

df = wide_var.merge(stats_by_strain, on="Strain", how="left")

# Sort by residual variance % for consistent ordering across plots
df = df.sort_values("Residual_Variance_pct", ascending=False).reset_index(drop=True)

# Shorten strain names for axis labels
def short_label(name, maxlen=40):
    name = str(name).replace("Influenza ", "")
    return name[:maxlen]

df["Strain_short"] = df["Strain"].apply(short_label)

EXPLAINED_LABEL = "Explained by Age, Sex & Cohort"
RESIDUAL_LABEL = "Unexplained / individual-level variance"

# ---------------------------------------------------------------
# PLOT 1: Variance decomposition -- stacked bar, all strains
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(df)), 7))

x = np.arange(len(df))
explained_pct = df["Cohort_Variance_pct"]
residual_pct = df["Residual_Variance_pct"]

ax.bar(x, explained_pct, width=0.6, label="Cohort (random intercept) variance",
       color="#4472C4", alpha=0.9)
ax.bar(x, residual_pct, width=0.6, bottom=explained_pct, label=RESIDUAL_LABEL,
       color="#ED7D31", alpha=0.85)

ax.set_xlabel("Strain", fontsize=12, fontweight="bold")
ax.set_ylabel("Share of unexplained variance (%)", fontsize=12, fontweight="bold")
ax.set_title("HAI Titer Variance Decomposition\n(Orange = variance not attributable to cohort, after age/sex are modeled)",
             fontsize=13, fontweight="bold", pad=20)
ax.set_xticks(x)
ax.set_xticklabels(df["Strain_short"], rotation=45, ha="right", fontsize=9)
ax.set_ylim([0, 105])
ax.legend(loc="upper right", fontsize=10)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
fig.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/01_variance_decomposition.png", dpi=300, bbox_inches="tight")
print("Saved: 01_variance_decomposition.png")
plt.close(fig)

# ---------------------------------------------------------------
# PLOT 2: Ranked horizontal bar -- Residual Variance % by strain
# ---------------------------------------------------------------
top_n = min(15, len(df))
top_strains = df.head(top_n).sort_values("Residual_Variance_pct")
y_pos = np.arange(len(top_strains))

colors = ["#d73027" if v > 60 else "#fee090" if v > 40 else "#91bfdb"
          for v in top_strains["Residual_Variance_pct"]]

fig, ax = plt.subplots(figsize=(12, max(5, 0.45 * len(top_strains))))
bars = ax.barh(y_pos, top_strains["Residual_Variance_pct"], color=colors,
                edgecolor="black", linewidth=0.8)

for i, (_, row) in enumerate(top_strains.iterrows()):
    ax.text(row["Residual_Variance_pct"] + 1, i, f"{row['Residual_Variance_pct']:.1f}%",
            va="center", fontsize=9, fontweight="bold")

ax.set_yticks(y_pos)
ax.set_yticklabels(top_strains["Strain_short"], fontsize=9)
ax.set_xlabel("Unexplained / individual-level variance (%)", fontsize=11, fontweight="bold")
ax.set_title("Strains Ranked by Unexplained Variance\nRed >60%, Yellow 40-60%, Blue <40%",
             fontsize=12, fontweight="bold", pad=15)
ax.set_xlim([0, 100])
ax.axvline(x=60, color="red", linestyle="--", linewidth=1.5, alpha=0.5)
ax.axvline(x=40, color="orange", linestyle="--", linewidth=1.5, alpha=0.5)
ax.grid(axis="x", alpha=0.3)

plt.tight_layout()
fig.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/02_residual_variance_ranking.png", dpi=300, bbox_inches="tight")
print("Saved: 02_residual_variance_ranking.png")
plt.close(fig)

# ---------------------------------------------------------------
# PLOT 3: Scatter -- Residual Variance % vs raw titer variability
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

scatter = ax.scatter(df["Std_HAI_log2"], df["Residual_Variance_pct"],
                      s=df["N_Participants"] * 4, alpha=0.65,
                      c=df["Residual_Variance_pct"], cmap="YlOrRd",
                      edgecolors="black", linewidth=0.5)

ax.set_xlabel("HAI Titer Variability (Std of log2 titers)", fontsize=11, fontweight="bold")
ax.set_ylabel("Unexplained / individual-level variance (%)", fontsize=11, fontweight="bold")
ax.set_title("Unexplained Variance vs. Raw HAI Titer Variability\n(Bubble size = N participants)",
             fontsize=12, fontweight="bold")
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Unexplained variance (%)", fontsize=10)

plt.tight_layout()
fig.savefig("/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/03_residual_vs_variability.png", dpi=300, bbox_inches="tight")
print("Saved: 03_residual_vs_variability.png")
plt.close(fig)

print("\n" + "=" * 70)
print("DONE -- 3 plots generated")
print("=" * 70)





#post results analysis 

STRAIN_COL = "Virus"
COHORT_COL = "Cohort"
OUTPUT_DIR = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_analysis_by_group"

def safe_name(s):
    """Make a string filesystem-safe for use as a folder/file name."""
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s

def plot_group(df, label, out_dir, folder_name=None):
    print(f"\n=== {label} (n={len(df)}) ===")
    print(df['log2_HAI'].describe())

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # log2 HAI distribution
    axes[0,0].hist(df['log2_HAI'], bins=10, alpha=0.7)
    axes[0,0].set_xlabel("Value (log2)")
    axes[0,0].set_ylabel("Count")
    axes[0,0].set_title("HAI distribution")

    # age vs HI
    axes[0,1].scatter(df["Age Reported_demo"], df["log2_HAI"], alpha=0.3)
    axes[0,1].set_xlabel("Age")
    axes[0,1].set_ylabel("log2 HAI")
    axes[0,1].set_title("Age vs HAI")

    # sex vs HI
    if df["Gender_demo"].nunique() > 1:
        df.boxplot(column="log2_HAI", by="Gender_demo", ax=axes[1,0])
        axes[1,0].set_title("HAI by sex")
        axes[1,0].set_xlabel("gender")
        axes[1,0].set_ylabel("log2 HAI")
    else:
        axes[1,0].axis("off")

    axes[1,1].axis("off")

    plt.suptitle(label)
    plt.tight_layout()

    # --- save to folder named by group ---
    if folder_name is None:
        folder_name = safe_name(label)
    group_folder = os.path.join(out_dir, folder_name)
    os.makedirs(group_folder, exist_ok=True)

    fig_path = os.path.join(group_folder, "summary_plots.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    stats_path = os.path.join(group_folder, "summary_stats.csv")
    df['log2_HAI'].describe().to_csv(stats_path)

    print(f"Saved -> {fig_path}")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- per strain ---
for strain, sub_df in day0.groupby(STRAIN_COL):
    plot_group(
        sub_df,
        f"Strain: {strain}",
        OUTPUT_DIR,
        folder_name=safe_name(strain)
    )