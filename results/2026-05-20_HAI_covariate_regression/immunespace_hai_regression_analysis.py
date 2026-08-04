#!{sys.executable} -m pip install ipykernel --upgrade --force-reinstall
#!{sys.executable} -m pip install statsmodels
#!{sys.executable} -m pip install matplotlib
#!{sys.executable} -m pip install patsy
#!{sys.executable} -m pip install upsetplot

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

import io
import base64
import statsmodels.formula.api as smf

#directories 

SCRATCH = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression"

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
    left_on=['new_participant_id', 'Study_ID', "Name"], 
    right_on=['Participant ID', 'Study_ID', "Cohort"],
    suffixes=('_merge2', '_merge1')    
)


#remove non influenza and unhealthy , and unvaccinated cohorts
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


final_merge["Race_hai"].value_counts(dropna=False) 

#add subtype mapping to final_merge dataframe
subtype_map = {
    "A/South Dakota/06/2007": "H1N1",
    "A/Uruguay/716/2007": "H3N2",
    "B/Florida/4/2006": "Yamagata",
    "A/Brisbane/59/2007": "H1N1",
    "B/Brisbane/3/2007": "Yamagata",
    "A/Solomon Islands/3/2006": "H1N1",
    "A/Wisconsin/67/2005": "H3N2",
    "B/Malaysia/2506/2004": "Victoria",
    "A/California/7/2009": "H1N1",
    "A/Perth/16/2009": "H3N2",
    "B/Brisbane/60/2008": "Victoria",
    "A/Indonesia/5/2005": "H5N1",
    "A/Brisbane/10/2007": "H3N2",
    "A/Victoria/361/2011": "H3N2",
    "B/Wisconsin/01/2010": "Yamagata",
    "B/Massachusetts/02/2012": "Yamagata",
    "A/Puerto Rico/8/1934": "H1N1",
    "A/Victoria/3/1975": "H3N2",
    "B/Lee/1940": "B Pre-lineage",
    "A/Texas/50/2012": "H3N2",
    "A/Perth/19/2009": "H3N2",
}



strain_fix = {
    "B/Wisonsin/01/2010": "B/Wisconsin/01/2010",
    "B/Massachusetts/2/2012": "B/Massachusetts/02/2012",
}

final_merge["Virus"] = final_merge["Virus"].replace(strain_fix)
final_merge["subtype"] = final_merge["Virus"].map(subtype_map)

final_merge.loc[final_merge["subtype"].isna(), "Virus"].unique()



#check that all participant IDs in final_merge are present in the hai table after filtering
set(final_merge['Participant ID_merge1']).issubset(set(hai['Participant ID']))

#additional filters
final_merge["Cohort for regression"] = final_merge["Description_merge2"].fillna(final_merge["Cohort"])

final_merge = final_merge[final_merge["Age Reported_demo"] >= 18]

cohorts_to_remove = ["18-30 year old monozygotic twins trivalent influenza vaccine","40-64 year old monozygotic twins trivalent influenza vaccine", "40-59 year old monozygotic twins trivalent influenza vaccine", "18-30 year old dizygotic twins trivalent influenza vaccine", "70-100 year old monozygotic or dizygotic twin pairs given  IIV3", "40-59 year old dizygotic twins trivalent influenza vaccine", "40-64 year old dizygotic twins trivalent influenza vaccine", "18-30 year old monozygotic twins pairs given  LAIV3", 
                     "18-30 year old dizygotic twin pairs given  LAIV3", "40-49 year old dizygotic twin pairs given  LAIV3", "40-49 year old monozygotic twin pairs given  LAIV3"]

final_merge = final_merge[~final_merge["Cohort for regression"].isin(cohorts_to_remove)]


mask_keep = (
    final_merge["Phenotype"].str.contains("Non-twin", case=False, na=False) |
    final_merge["Phenotype"].str.contains("Non-Twin", case=False, na=False) |
    ~final_merge["Phenotype"].str.contains("twin", case=False, na=False)
)

final_merge = final_merge[mask_keep]





OUTPUT_ROOT = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/immunespace_hai_regression_results"
os.makedirs(OUTPUT_ROOT, exist_ok=True)



### pre HAI demographics and distribution plots
####plots needed for lab notebook
 
"""
UpSet plots: Participants x Days, computed separately per strain.

Sets     = Study Time Collected (Day 0, Day 28, ...)
Elements = participants (new_participant_id)
Split by = strain (subtype, or Virus for finer granularity)

For each strain, this shows how many participants were sampled at each day,
and which combinations of days each participant was sampled at -- so you can
compare e.g. "H1N1 day-overlap" vs "H3N2 day-overlap".


"""

import pandas as pd
from upsetplot import from_memberships, UpSet
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Load your data here. Replace this with your actual source (csv, etc.)


# df = pd.read_csv("your_file.csv")
# ------------------------------------------------------------------
df = final_merge.copy()  # assumes `df` is already the dataframe shown in your message

# use `subtype` (H1N1/H3N2/Yamagata/...) for the strain grouping.
# swap to "Virus" instead if you want individual reference-strain granularity
# (e.g. A/South Dakota/06/2007 vs A/Uruguay/716/2007) rather than subtype.
strain_col = "Virus"

# one row per (participant, day, strain), regardless of everything else
part_day_strain = (
    df[["Participant ID_merge1", "Study Time Collected", strain_col]]
    .dropna(subset=["Study Time Collected", strain_col])
    .drop_duplicates()
)

strains = sorted(part_day_strain[strain_col].unique())

for strain in strains:
    subset = part_day_strain[part_day_strain[strain_col] == strain]

    # for each participant, which set of days were they sampled at for this strain?
    memberships = (
        subset.groupby("Participant ID_merge1")["Study Time Collected"]
        .apply(lambda s: [f"Day {int(d)}" for d in sorted(s.unique())])
        .tolist()
    )

    upset_data = from_memberships(memberships)

    fig = plt.figure(figsize=(10, 6))
    UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
    ).plot(fig=fig)
    plt.suptitle(f"Participant overlap across collection days -- {strain}")

    safe_name = str(strain).replace("/", "-").replace(" ", "_")
    plt.savefig(OUTPUT_ROOT + f"/upset_participants_by_day_{safe_name}.png",
        dpi=150,
        bbox_inches="tight",
    )

plt.show()




### pre HAI demographics and distribution plots
### One PNG per strain: 3 facets (distribution / age vs HAI / HAI by sex),
### with Day plotted as different colors (shared legend) within each facet.
### PNGs are organized into subfolders by subtype for easy cross-strain comparison.

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd

# ---- config ----
STRAIN_COL = "Virus"
COHORT_COL = "Cohort"
DAY_COL = "Study Time Collected"
DAY_UNIT_COL = "Study Time Collected Unit"
SUBTYPE_COL = "subtype"
AGE_COL = "Age Reported_hai"
GENDER_COL = "Gender_hai"
HAI_COL = "log2_HAI"
MIN_N = 15  # per strain/day slice; days below this are dropped before plotting


outdir = os.path.join(OUTPUT_ROOT, "final_hai_strain_faceted_demo_data_png")
os.makedirs(outdir, exist_ok=True)

sns.set_style("whitegrid")


def slugify(s):
    return (
        str(s)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("—", "_")
        .replace(":", "")
    )


def make_strain_facet_png(strain_df, strain, subtype, day_col, out_path):
    """Build a 1x4 facet PNG for one strain: distribution, age vs HAI, HAI by sex, HAI by phenotype.
    Days are colored consistently across all four panels with one shared legend.
    """
    days = sorted(strain_df[day_col].dropna().unique())
    palette = dict(zip(days, sns.color_palette("tab10", n_colors=len(days))))

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # --- Panel 1: HAI distribution, colored by day ---
    # Outlined step histograms (no fill) so overlapping day-colors stay
    # readable as distinct lines instead of blending into a solid mass.
    sns.histplot(
        data=strain_df,
        x=HAI_COL,
        hue=day_col,
        palette=palette,
        multiple="layer",
        fill=False,
        element="step",
        linewidth=2.2,
        bins=10,
        ax=axes[0],
        legend=False,
    )
    axes[0].set_title("HAI distribution")
    axes[0].set_xlabel("log2 HAI")
    axes[0].set_ylabel("Count")

    # --- Panel 2: Age vs HAI, colored by day ---
    sns.scatterplot(
        data=strain_df,
        x=AGE_COL,
        y=HAI_COL,
        hue=day_col,
        palette=palette,
        alpha=0.5,
        ax=axes[1],
        legend=False,
    )
    axes[1].set_title("Age vs HAI")
    axes[1].set_xlabel("Age")
    axes[1].set_ylabel("log2 HAI")

    # --- Panel 3: HAI by sex, colored by day (grouped boxplot) ---
    sns.boxplot(
        data=strain_df,
        x=GENDER_COL,
        y=HAI_COL,
        hue=day_col,
        palette=palette,
        ax=axes[2],
    )
    axes[2].set_title("HAI by sex")
    axes[2].set_xlabel("Gender")
    axes[2].set_ylabel("log2 HAI")
    axes[2].legend_.remove()  # remove per-axis legend, add one shared legend below
    


    # --- shared legend for Day, placed once for the whole figure ---
    handles = [mpatches.Patch(color=palette[d], label=f"Day {d}") for d in days]
    fig.legend(
        handles=handles,
        title=day_col,
        loc="lower center",
        ncol=min(len(days), 8),
        bbox_to_anchor=(0.5, -0.05),
        frameon=False,
    )

    n_total = len(strain_df)
    fig.suptitle(f"{strain}  (Subtype: {subtype}, N={n_total})", fontsize=14, y=1.03)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---- build one PNG per strain, grouped into subtype subfolders ----
final_merge["_subtype_filled"] = (
    final_merge[SUBTYPE_COL].fillna("NA") if SUBTYPE_COL in final_merge.columns else "NA"
)

subtypes = sorted(final_merge["_subtype_filled"].dropna().unique())
saved_files = []

for subtype in subtypes:
    subtype_df = final_merge[final_merge["_subtype_filled"] == subtype]
    subtype_dir = os.path.join(outdir, f"subtype_{slugify(subtype)}")
    os.makedirs(subtype_dir, exist_ok=True)

    strains = sorted(subtype_df[STRAIN_COL].dropna().unique())

    for strain in strains:
        strain_df_full = subtype_df[subtype_df[STRAIN_COL] == strain].copy()

        # drop days with too few samples, same as before, but keep the rest together
        day_counts = strain_df_full[DAY_COL].value_counts()
        valid_days = day_counts[day_counts >= MIN_N].index
        strain_df = strain_df_full[strain_df_full[DAY_COL].isin(valid_days)]

        if strain_df.empty:
            print(f"Skipping {strain} (subtype {subtype}): no day group with N >= {MIN_N}")
            continue

        dropped = set(strain_df_full[DAY_COL].dropna().unique()) - set(valid_days)
        if dropped:
            print(f"{strain}: dropping days {sorted(dropped)} (N < {MIN_N})")

        fname = f"{slugify(strain)}.png"
        out_path = os.path.join(subtype_dir, fname)

        make_strain_facet_png(strain_df, strain, subtype, DAY_COL, out_path)
        saved_files.append(out_path)
        print(f"Saved: {out_path}")

print(f"\nDone. {len(saved_files)} PNGs saved under {outdir}")



#table of  phenotype value counts

final_merge["Phenotype"].value_counts()
















"""
HAI regression + variance analysis, organized by Subtype -> Strain -> Day
===========================================================================
Builds on the original per-strain-per-day LMM/OLS logic, but:

  1. Everything is grouped/output by SUBTYPE -> STRAIN -> DAY, not just
     strain-day, so results for the same subtype sit together.
  2. Every strain gets an annotation block: total N samples, N distinct
     studies, N distinct cohorts (computed once, across ALL days/visits
     for that strain).
  3. New "Std Dev vs Residual Variance" bubble plot, with the N for each
     group written directly next to its bubble (not just encoded in size).
  4. All descriptive + model-fit + variance plots for a subtype/strain/day
     are collected into ONE self-contained HTML report you can open from
     a single link, with a jump-to table of contents.

Run this after `final_merge` is already in memory (same as your original
script).
"""





# ---------------------------------------------------------------
# 1. LOAD + CLEAN (same as original)
# ---------------------------------------------------------------
df = final_merge.copy()
df.columns = df.columns.str.strip()

COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_hai",
    "sex": "Gender_demo",
    "cohort": "Cohort for regression",
    "virus": "Virus",
    "subtype": "subtype",
    "study_id": "Study_ID",
    "participant_id": "Participant ID",
    "day": "Study Time Collected",
}

clean = df.rename(columns={
    COLUMNS["outcome"]: "log2_HAI",
    COLUMNS["age"]: "age",
    COLUMNS["sex"]: "sex",
    COLUMNS["cohort"]: "Cohort for regression",
    COLUMNS["virus"]: "virus",
    COLUMNS["subtype"]: "subtype",
    COLUMNS["study_id"]: "study_id",
    COLUMNS["participant_id"]: "participant_id",
    COLUMNS["day"]: "day",
})

needed = ["log2_HAI", "age", "sex", "Cohort for regression", "virus", "subtype", "study_id", "day"]
missing = [c for c in needed if c not in clean.columns]
if missing:
    raise KeyError(f"Rename didn't produce expected columns: {missing}. "
                    f"Available columns are: {df.columns.tolist()}")

clean = clean.dropna(subset=needed).copy()
clean["sex"] = clean["sex"].astype("category")
clean["cohort"] = clean["Cohort for regression"].astype("category")
clean["subtype"] = clean["subtype"].astype("category")
clean["age"] = pd.to_numeric(clean["age"], errors="coerce")
clean["day"] = pd.to_numeric(clean["day"], errors="coerce")
clean = clean.dropna(subset=["age", "day"])

if "Study Time Collected Unit" in df.columns:
    bad_units = df.loc[clean.index, "Study Time Collected Unit"].dropna().unique()
    bad_units = [u for u in bad_units if str(u).strip().lower() != "days"]
    if bad_units:
        raise ValueError(
            f"'Study Time Collected' isn't uniformly in Days -- found units: {bad_units}."
        )

clean["vaccinated"] = np.where(clean["day"] >= 14, "Vaccinated", "Not Vaccinated")
clean["vaccinated"] = clean["vaccinated"].astype("category")


def mixedlm_r2(fit):
    fixed_effects_fitted = fit.model.exog @ fit.fe_params
    var_fixed = np.var(fixed_effects_fitted, ddof=0)
    var_random = fit.cov_re.iloc[0, 0]
    var_resid = fit.scale
    total = var_fixed + var_random + var_resid
    return var_fixed / total, (var_fixed + var_random) / total


def safe_name(s):
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s


# ---------------------------------------------------------------
# 2. STRAIN-LEVEL ANNOTATION COUNTS
#    (computed ACROSS ALL DAYS for that strain -- this is the
#    "how many samples / studies / cohorts" annotation)
# ---------------------------------------------------------------
strain_counts = (
    clean.groupby(["subtype", "virus"])
    .agg(
        N_Samples=("log2_HAI", "count"),
        N_Studies=("study_id", "nunique"),
        N_Cohorts=("Cohort for regression", "nunique"),
        N_Days=("day", "nunique"),
    )
    .reset_index()
)
strain_counts.to_csv(os.path.join(OUTPUT_ROOT, "strain_level_counts.csv"), index=False)
print("Strain-level counts (samples / studies / cohorts):")
print(strain_counts.to_string(index=False))


def strain_annotation_text(virus_name):
    row = strain_counts[strain_counts["virus"] == virus_name]
    if row.empty:
        return "N_Samples=?, N_Studies=?, N_Cohorts=?"
    r = row.iloc[0]
    return f"N_Samples={r['N_Samples']} | N_Studies={r['N_Studies']} | N_Cohorts={r['N_Cohorts']} | N_Days={r['N_Days']}"


# ---------------------------------------------------------------
# 3. PER (subtype, virus, day) MODEL FIT -- same decision rule as before
# ---------------------------------------------------------------
lmm_results, ols_results, model_choice, r2_values, vacc_label = {}, {}, {}, {}, {}

for (subtype_name, virus_name, day_val), sub in clean.groupby(["subtype", "virus", "day"]):
    sub = sub.copy()
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)
    vacc_status = sub["vaccinated"].iloc[0]
    key = (subtype_name, virus_name, day_val)
    vacc_label[key] = vacc_status

    if n_obs < 20: #If the number of observations is less than 20, skip this iteration. 
        continue

    if n_cohorts >= 2:
        model_choice[key] = "LMM"
        model = smf.mixedlm("log2_HAI ~ age + sex", data=sub, groups=sub["cohort"], re_formula="1")
        try:
            fit = model.fit(reml=True)
            lmm_results[key] = fit
            r2_values[key] = mixedlm_r2(fit)
        except Exception as e:
            print(f"  -> LMM failed for {key}: {e}")
    else:
        model_choice[key] = "OLS"
        try:
            ols_fit = smf.ols("log2_HAI ~ age + sex", data=sub).fit()
            ols_results[key] = ols_fit
            r2_values[key] = (ols_fit.rsquared, ols_fit.rsquared)
        except Exception as e:
            print(f"  -> OLS failed for {key}: {e}")

# ---------------------------------------------------------------
# 4. UNIFIED SUMMARY TABLE (now carries subtype explicitly)
# ---------------------------------------------------------------
rows = []
for (subtype_name, virus_name, day_val), fit in lmm_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()
    for term in params.index:
        if term == "Group Var":
            continue
        rows.append({
            "Subtype": subtype_name, "Virus": virus_name, "Day": day_val,
            "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
            "Model": "LMM", "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4), "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Subtype": subtype_name, "Virus": virus_name, "Day": day_val,
        "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
        "Model": "LMM", "Term": "Cohort (random intercept) variance",
        "Estimate": round(fit.cov_re.iloc[0, 0], 4), "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })
    rows.append({
        "Subtype": subtype_name, "Virus": virus_name, "Day": day_val,
        "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
        "Model": "LMM", "Term": "Residual variance",
        "Estimate": round(fit.scale, 4), "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })

for (subtype_name, virus_name, day_val), fit in ols_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()
    for term in params.index:
        rows.append({
            "Subtype": subtype_name, "Virus": virus_name, "Day": day_val,
            "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
            "Model": "OLS", "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4), "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Subtype": subtype_name, "Virus": virus_name, "Day": day_val,
        "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
        "Model": "OLS", "Term": "Residual variance",
        "Estimate": round(fit.mse_resid, 4), "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })

summary_df = pd.DataFrame(rows)
summary_df.to_csv(os.path.join(OUTPUT_ROOT, "hai_regression_summary.csv"), index=False)

# ---------------------------------------------------------------
# 5. MODEL CHOICE LOG
# ---------------------------------------------------------------
choice_rows = []
for (subtype_name, virus_name, day_val), sub in clean.groupby(["subtype", "virus", "day"]):
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)
    key = (subtype_name, virus_name, day_val)
    marg_r2, cond_r2 = r2_values.get(key, (np.nan, np.nan))
    choice_rows.append({
        "Subtype": subtype_name, "Virus": virus_name, "Day": day_val,
        "Vaccinated": vacc_label.get(key, sub["vaccinated"].iloc[0]),
        "N_Obs": n_obs, "N_Cohorts": n_cohorts,
        "Model_Used": model_choice.get(key, "Skipped (n_obs<20)"),
        "Marginal_R2": round(marg_r2, 4) if pd.notna(marg_r2) else np.nan,
        "Conditional_R2": round(cond_r2, 4) if pd.notna(cond_r2) else np.nan,
    })
choice_df = pd.DataFrame(choice_rows)
choice_df.to_csv(os.path.join(OUTPUT_ROOT, "hai_model_choice_log.csv"), index=False)

# ---------------------------------------------------------------
# 6. VARIANCE DECOMPOSITION TABLE (adds subtype + strain-level N annotation)
# ---------------------------------------------------------------
var_rows = summary_df[summary_df["Term"].isin(
    ["Cohort (random intercept) variance", "Residual variance"]
)].copy()
wide_var = var_rows.pivot_table(
    index=["Subtype", "Virus", "Day", "Vaccinated"], columns="Term", values="Estimate"
).reset_index()
wide_var = wide_var.rename(columns={
    "Cohort (random intercept) variance": "Cohort_Variance",
    "Residual variance": "Residual_Variance",
})
if "Cohort_Variance" not in wide_var.columns:
    wide_var["Cohort_Variance"] = 0.0
wide_var["Cohort_Variance"] = wide_var["Cohort_Variance"].fillna(0.0)
wide_var["Total_Variance"] = wide_var["Cohort_Variance"] + wide_var["Residual_Variance"]
wide_var["Cohort_Variance_pct"] = 100 * wide_var["Cohort_Variance"] / wide_var["Total_Variance"]
wide_var["Residual_Variance_pct"] = 100 * wide_var["Residual_Variance"] / wide_var["Total_Variance"]

stats_by_group = (
    clean.groupby(["subtype", "virus", "day"])["log2_HAI"]
    .agg(Std_HAI_log2="std", N_Participants="count")
    .reset_index()
    .rename(columns={"subtype": "Subtype", "virus": "Virus", "day": "Day"})
)
vdf = wide_var.merge(stats_by_group, on=["Subtype", "Virus", "Day"], how="left")
vdf = vdf.merge(strain_counts.rename(columns={"virus": "Virus", "subtype": "Subtype"}),
                 on=["Subtype", "Virus"], how="left")


def short_label(virus, day, maxlen=28):
    v = str(virus).replace("Influenza ", "")[:maxlen]
    return f"{v} (D{int(day)})"


vdf["Group_short"] = vdf.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)
# label used on every plot below: strain + day + subtype, so a bubble/bar is identifiable at a glance
vdf["Group_label_full"] = vdf.apply(
    lambda r: f"{short_label(r['Virus'], r['Day'])} [{r['Subtype']}]", axis=1
)
vdf.to_csv(os.path.join(OUTPUT_ROOT, "variance_decomposition_table.csv"), index=False)

overview_plot_paths = []  # collected so they can be embedded in the HTML report



# =================================================================
# 7a. PLOT: Grouped bar chart -- Marginal vs Conditional R² per group
# Output: one PNG per subtype
# =================================================================

df_plot = choice_df.dropna(subset=["Model_Used"]).copy()
df_plot = df_plot[df_plot["Model_Used"] != "Skipped (n_obs<20)"]

df_plot["Group_short"] = df_plot.apply(
    lambda r: short_label(r["Virus"], r["Day"]), axis=1
)
df_plot["Group_label_full"] = df_plot.apply(
    lambda r: f"{short_label(r['Virus'], r['Day'])} [{r['Subtype']}]",
    axis=1,
)

for subtype in sorted(df_plot["Subtype"].unique()):

    sub = (
        df_plot[df_plot["Subtype"] == subtype]
        .sort_values("Conditional_R2", ascending=False)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(max(8, 1.8 * len(sub)), 7))

    x = np.arange(len(sub))
    width = 0.35

    bars1 = ax.bar(
        x - width / 2,
        sub["Marginal_R2"],
        width,
        label="Marginal R²",
        color="#4472C4",
        edgecolor="black",
        linewidth=0.6,
    )

    bars2 = ax.bar(
        x + width / 2,
        sub["Conditional_R2"],
        width,
        label="Conditional R²",
        color="#ED7D31",
        edgecolor="black",
        linewidth=0.6,
    )

    for b in bars1:
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.02,
            f"{b.get_height():.3f}",
            ha="center",
            fontsize=9,
        )

    for b in bars2:
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.02,
            f"{b.get_height():.3f}",
            ha="center",
            fontsize=9,
        )

    for i, (model_used, vacc) in enumerate(zip(sub["Model_Used"], sub["Vaccinated"])):
        ax.text(
            i,
            -0.06,
            f"[{model_used} | {vacc}]",
            ha="center",
            va="top",
            fontsize=7.5,
            color="gray",
            style="italic",
        )

    ax.set_xlabel("Strain (Day)", fontsize=12, fontweight="bold")
    ax.set_ylabel("R²", fontsize=12, fontweight="bold")
    ax.set_title(
        f"Marginal vs Conditional R²\nSubtype: {subtype}",
        fontsize=13,
        fontweight="bold",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(sub["Group_short"], rotation=45, ha="right", fontsize=8.5)
    ax.set_ylim([0, 1.0])
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT,
        f"01_r2_marginal_vs_conditional_{subtype.replace('/', '_').replace(' ', '_')}.png",
    )

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    overview_plot_paths.append(
        (f"Marginal vs Conditional R² ({subtype})", outfile)
    )

    print(f"Saved: {outfile}")


# =================================================================
# 7b. PLOT: Gap plot -- how much cohort (random effect) adds, by Subtype
# =================================================================
df_plot = choice_df.dropna(subset=["Model_Used", "Subtype"]).copy()
df_plot = df_plot[df_plot["Model_Used"] != "Skipped (n_obs<20)"].copy()

df_plot["Group_short"] = df_plot.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)
df_plot["Group_label_full"] = df_plot.apply(
    lambda r: f"{short_label(r['Virus'], r['Day'])} [{r['Subtype']}]",
    axis=1
)
df_plot["R2_gap"] = df_plot["Conditional_R2"] - df_plot["Marginal_R2"]

for subtype in sorted(df_plot["Subtype"].dropna().unique()):

    sub = (
        df_plot[df_plot["Subtype"] == subtype]
        .sort_values("R2_gap", ascending=True)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(sub))))

    y_pos = np.arange(len(sub))
    colors = ["#91bfdb" if g == 0 else "#4472C4" for g in sub["R2_gap"]]

    ax.barh(y_pos, sub["R2_gap"], color=colors, edgecolor="black", linewidth=0.7)

    for i, (_, row) in enumerate(sub.iterrows()):
        label = f"{row['R2_gap']:.3f}" if row["R2_gap"] > 0 else "0 (OLS, no random effect)"
        ax.text(row["R2_gap"] + 0.01, i, label, va="center", fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([short_label(v, d) for v, d in zip(sub["Virus"], sub["Day"])], fontsize=8.5)
    ax.set_xlabel("Conditional R² − Marginal R²", fontsize=11, fontweight="bold")
    ax.set_title(f"Cohort Contribution to R²\nSubtype: {subtype}", fontsize=12, fontweight="bold")
    ax.set_xlim([0, max(0.1, df_plot["R2_gap"].max() * 1.3)])
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()

    safe_subtype = subtype.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")
    outfile = os.path.join(OUTPUT_ROOT, f"02_r2_cohort_contribution_gap_{safe_subtype}.png")
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {outfile}")
    
    

# =================================================================
# 7c. PLOT: Variance decomposition -- stacked bar, faceted by Subtype
# =================================================================
for subtype in sorted(vdf["Subtype"].unique()):

    sub = (
        vdf[vdf["Subtype"] == subtype]
        .sort_values("Residual_Variance_pct", ascending=False)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(sub)), 7))

    x = np.arange(len(sub))

    ax.bar(
        x,
        sub["Cohort_Variance_pct"],
        width=0.6,
        color="#4472C4",
        label="Cohort variance"
    )

    ax.bar(
        x,
        sub["Residual_Variance_pct"],
        width=0.6,
        bottom=sub["Cohort_Variance_pct"],
        color="#ED7D31",
        label="Residual variance"
    )

    for i, vacc in enumerate(sub["Vaccinated"]):
        ax.text(i, -6, vacc,
                ha="center", fontsize=7,
                color="gray", style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(sub["Group_label_full"],
                       rotation=45, ha="right")

    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    ax.set_title(f"Variance Decomposition\n{subtype}")

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT,
        f"03_variance_decomposition_{subtype}.png"
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {outfile}")



# =================================================================
# 7d. PLOT: Ranked horizontal bar -- Residual Variance % by group, faceted by Subtype
# =================================================================
for subtype in sorted(vdf["Subtype"].unique()):

    sub = vdf[vdf["Subtype"] == subtype]

    top_groups = (
        sub.sort_values("Residual_Variance_pct", ascending=False)
        .head(15)
        .sort_values("Residual_Variance_pct")
    )

    fig, ax = plt.subplots(
        figsize=(12, max(5, 0.45 * len(top_groups)))
    )

    y_pos = np.arange(len(top_groups))

    colors = [
        "#d73027" if v > 60
        else "#fee090" if v > 40
        else "#91bfdb"
        for v in top_groups["Residual_Variance_pct"]
    ]

    ax.barh(
        y_pos,
        top_groups["Residual_Variance_pct"],
        color=colors,
        edgecolor="black"
    )

    for i, (_, row) in enumerate(top_groups.iterrows()):
        ax.text(
            row["Residual_Variance_pct"] + 1,
            i,
            f"{row['Residual_Variance_pct']:.1f}%",
            va="center",
            fontsize=9,
            fontweight="bold"
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_groups["Group_label_full"])
    ax.set_xlim(0, 100)

    ax.axvline(40, ls="--", color="orange")
    ax.axvline(60, ls="--", color="red")

    ax.grid(axis="x", alpha=0.3)
    ax.set_title(f"Residual Variance Ranking\n{subtype}")

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT,
        f"04_residual_variance_ranking_{subtype}.png"
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {outfile}")



# =================================================================
# 7e. Std Dev vs Residual Variance bubble plot -- colored by strain,
#     shaped by day, sized by N, labeled directly on every bubble
# =================================================================
import itertools

MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]

def plot_std_vs_residual(sub_vdf, title, out_path):
    fig, ax = plt.subplots(figsize=(14, 8))

    strains = sorted(sub_vdf["Virus"].dropna().unique())
    days = sorted(sub_vdf["Day"].dropna().unique())

    strain_colors = plt.cm.tab20(np.linspace(0, 1, max(len(strains), 1)))
    strain_color_map = dict(zip(strains, strain_colors))

    marker_cycle = itertools.cycle(MARKER_CYCLE)
    day_marker_map = {day: next(marker_cycle) for day in days}

    # plot one small scatter per (strain, day) group so marker shape can vary
    for (virus, day), grp in sub_vdf.groupby(["Virus", "Day"]):
        ax.scatter(
            grp["Std_HAI_log2"], grp["Residual_Variance"],
            s=grp["N_Participants"] * 6, alpha=0.75,
            color=strain_color_map[virus],
            marker=day_marker_map[day],
            edgecolors="black", linewidth=0.6,
        )

    # annotate each bubble with its strain, day, and N
    for _, row in sub_vdf.iterrows():
        ax.annotate(
            f"{row['Group_short']}\nN={int(row['N_Participants'])}",
            (row["Std_HAI_log2"], row["Residual_Variance"]),
            textcoords="offset points", xytext=(7, 5), fontsize=7,
            color="black", linespacing=1.3,
        )

    ax.set_xlabel("HAI Titer Std Dev (log2 scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Residual Variance (unexplained, after age/sex/cohort adjustment)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)

    # --- Legend 1: Strain (color) ---
    strain_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=strain_color_map[s], markeredgecolor="black",
                   markersize=8, label=s)
        for s in strains
    ]
    legend1 = ax.legend(
        handles=strain_handles, title="Strain",
        loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=True, fontsize=8,
    )
    ax.add_artist(legend1)

    # --- Legend 2: Day (marker shape) ---
    day_handles = [
        plt.Line2D([0], [0], marker=day_marker_map[d], color="gray",
                   linestyle="None", markeredgecolor="black",
                   markersize=8, label=f"Day {d}")
        for d in days
    ]
    legend2 = ax.legend(
        handles=day_handles, title="Day",
        loc="upper left", bbox_to_anchor=(1.02, 0.55),
        frameon=True, fontsize=8,
    )
    ax.add_artist(legend2)

    # --- Legend 3: Sample size (bubble size) ---
    size_values = [10, 25, 50, 100]
    size_values = [n for n in size_values if n <= sub_vdf["N_Participants"].max()]
    size_handles = [
        plt.scatter([], [], s=n * 6, color="gray", alpha=0.75,
                    edgecolors="black", linewidth=0.6, label=f"N={n}")
        for n in size_values
    ]
    ax.legend(
        handles=size_handles, title="Sample size",
        loc="upper left", bbox_to_anchor=(1.02, 0.15),
        frameon=True, fontsize=8,
    )

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


p5 = os.path.join(OUTPUT_ROOT, "05_stddev_vs_residual_all.png")
plot_std_vs_residual(vdf, "Std Dev vs Residual Variance (color=strain, shape=day, size=N)", p5)
overview_plot_paths.append(("Std Dev vs Residual Variance — All Subtypes", p5))

for st in vdf["Subtype"].unique():
    sub = vdf[vdf["Subtype"] == st]
    p_st = os.path.join(OUTPUT_ROOT, f"05_stddev_vs_residual_{safe_name(st)}.png")
    plot_std_vs_residual(sub, f"Std Dev vs Residual Variance — Subtype {st}", p_st)
    overview_plot_paths.append((f"Std Dev vs Residual Variance — Subtype {st}", p_st))




# =================================================================
# 8. DESCRIPTIVE PLOTS per (subtype, virus, day), saved into a
#    subtype/strain/day folder tree
# =================================================================
BASE_DESC_DIR = os.path.join(OUTPUT_ROOT, "descriptive_by_subtype_strain_day")
os.makedirs(BASE_DESC_DIR, exist_ok=True)


def plot_group_matplotlib(sub_df, label, out_dir):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].hist(sub_df["log2_HAI"], bins=10, alpha=0.7, color="#4472C4", edgecolor="black")
    axes[0, 0].set_xlabel("Value (log2)")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("HAI distribution")

    axes[0, 1].scatter(sub_df["age"], sub_df["log2_HAI"], alpha=0.3, color="#ED7D31")
    axes[0, 1].set_xlabel("Age")
    axes[0, 1].set_ylabel("log2 HAI")
    axes[0, 1].set_title("Age vs HAI")

    if sub_df["sex"].nunique() > 1:
        sub_df.boxplot(column="log2_HAI", by="sex", ax=axes[1, 0])
        axes[1, 0].set_title("HAI by sex")
        axes[1, 0].set_xlabel("sex")
        axes[1, 0].set_ylabel("log2 HAI")
        fig.suptitle("")
    else:
        axes[1, 0].axis("off")

    axes[1, 1].axis("off")
    fig.suptitle(label, fontsize=11)
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, "summary_plots.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    sub_df["log2_HAI"].describe().to_csv(os.path.join(out_dir, "summary_stats.csv"))
    return fig_path


desc_plot_paths = {}  # key -> png path, for embedding in the HTML report

for (subtype_name, virus_name, day_val), sub_df in clean.groupby(["subtype", "virus", "day"]):
    vacc_status = sub_df["vaccinated"].iloc[0]
    n_obs = len(sub_df)
    label = f"{virus_name} | Subtype {subtype_name} | Day {day_val} | {vacc_status} | N={n_obs}"

    out_dir = os.path.join(
        BASE_DESC_DIR, safe_name(subtype_name), safe_name(virus_name),
        f"Day{int(day_val)}__{safe_name(vacc_status)}",
    )
    fig_path = plot_group_matplotlib(sub_df, label, out_dir)
    desc_plot_paths[(subtype_name, virus_name, day_val)] = fig_path

print(f"\nDescriptive plots saved under: {BASE_DESC_DIR}")


# ---------------------------------------------------------------
# PLOT: Scatter -- Residual Variance % vs raw titer variability
# One PNG per subtype, colored by strain, shaped by day,
# sized by N, with a sample-size legend
# ---------------------------------------------------------------
import itertools

def safe_name(s):
    return (
        str(s)
        .replace("/", "_")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
    )

# Marker shapes to cycle through for "Day"
MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]

for subtype in sorted(vdf["Subtype"].dropna().unique()):

    sub = vdf[vdf["Subtype"] == subtype].copy()

    fig, ax = plt.subplots(figsize=(10, 6))

    # --- build strain -> color and day -> marker mappings ---
    strains = sorted(sub["Virus"].dropna().unique())
    days = sorted(sub["Day"].dropna().unique())

    strain_colors = plt.cm.tab20(np.linspace(0, 1, max(len(strains), 1)))
    strain_color_map = dict(zip(strains, strain_colors))

    marker_cycle = itertools.cycle(MARKER_CYCLE)
    day_marker_map = {day: next(marker_cycle) for day in days}

    # --- plot one small scatter per (strain, day) group ---
    for (virus, day), grp in sub.groupby(["Virus", "Day"]):
        ax.scatter(
            grp["Std_HAI_log2"],
            grp["Residual_Variance_pct"],
            s=grp["N_Participants"] * 4,
            alpha=0.75,
            color=strain_color_map[virus],
            marker=day_marker_map[day],
            edgecolors="black",
            linewidth=0.5,
        )

    # Point labels: strain + day, plus sample size
    for _, row in sub.iterrows():
        ax.annotate(
            f"{short_label(row['Virus'], row['Day'])}\nN={row['N_Participants']}",
            (row["Std_HAI_log2"], row["Residual_Variance_pct"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7.5,
        )

    ax.set_xlabel("HAI Titer Std Dev (log2 scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Unexplained variance (%)", fontsize=11, fontweight="bold")
    ax.set_title(
        f"Raw HAI Titer Variability vs. Residual Variance (%) After Covariate Adjustment\nSubtype: {subtype}",
        fontsize=12,
        fontweight="bold",
    )
    ax.grid(True, alpha=0.3)

    # --- Legend 1: Strain (color) ---
    strain_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=strain_color_map[s], markeredgecolor="black",
                   markersize=8, label=s)
        for s in strains
    ]
    legend1 = ax.legend(
        handles=strain_handles,
        title="Strain",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        fontsize=8,
    )
    ax.add_artist(legend1)

    # --- Legend 2: Day (marker shape) ---
    day_handles = [
        plt.Line2D([0], [0], marker=day_marker_map[d], color="gray",
                   linestyle="None", markeredgecolor="black",
                   markersize=8, label=f"Day {d}")
        for d in days
    ]
    legend2 = ax.legend(
        handles=day_handles,
        title="Day",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.55),
        frameon=True,
        fontsize=8,
    )
    ax.add_artist(legend2)

    # --- Legend 3: Sample size (bubble size) ---
    size_values = [10, 25, 50, 100]
    size_values = [n for n in size_values if n <= sub["N_Participants"].max()]
    size_handles = [
        plt.scatter([], [], s=n * 4, color="gray", alpha=0.75,
                    edgecolors="black", linewidth=0.5, label=f"N={n}")
        for n in size_values
    ]
    ax.legend(
        handles=size_handles,
        title="Sample size",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.15),
        frameon=True,
        fontsize=8,
    )

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT,
        f"03_residual(relative)_vs_std_dev_{safe_name(subtype)}.png"
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

print(f"Saved: {outfile}")




"""
HAI model diagnostics: QQ plots, residual checks, and performance metrics
===========================================================================
Subtype-aware version: saves diagnostics into subtype-specific folders and
adds Subtype to the performance summary.
"""
# ---------------------------------------------------------------
# Diagnostics: use the same subtype source as the rest of the script
# ---------------------------------------------------------------
DIAG_DIR = os.path.join(OUTPUT_ROOT, "final_model_diagnostics_by_subtype")
os.makedirs(DIAG_DIR, exist_ok=True)

def safe_name(s):
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s

def get_resid_fitted(fit):
    resid = np.asarray(fit.resid)
    fitted = np.asarray(fit.fittedvalues)
    return resid, fitted

def diagnostic_plots(subtype, virus_name, day_val, vacc_status, fit, model_type, out_dir):
    resid, fitted = get_resid_fitted(fit)
    std_resid = (resid - resid.mean()) / resid.std(ddof=1)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    ax = axes[0, 0]
    ax.scatter(fitted, resid, alpha=0.5, edgecolor="black", linewidth=0.3)
    ax.axhline(0, color="red", linestyle="--", linewidth=1)
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess as _lowess
        sm_line = _lowess(resid, fitted, frac=0.6)
        ax.plot(sm_line[:, 0], sm_line[:, 1], color="blue", linewidth=1.5)
    except Exception:
        pass
    ax.set_xlabel("Fitted values")
    ax.set_ylabel("Residuals")
    ax.set_title("Residuals vs Fitted")
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    stats.probplot(std_resid, dist="norm", plot=ax)
    ax.set_title("Normal Q-Q")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    sqrt_abs_std_resid = np.sqrt(np.abs(std_resid))
    ax.scatter(fitted, sqrt_abs_std_resid, alpha=0.5, edgecolor="black", linewidth=0.3)
    try:
        sm_line2 = _lowess(sqrt_abs_std_resid, fitted, frac=0.6)
        ax.plot(sm_line2[:, 0], sm_line2[:, 1], color="blue", linewidth=1.5)
    except Exception:
        pass
    ax.set_xlabel("Fitted values")
    ax.set_ylabel("sqrt(|Standardized Residuals|)")
    ax.set_title("Scale-Location")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.hist(resid, bins=15, density=True, alpha=0.7, color="#4472C4", edgecolor="black")
    xs = np.linspace(resid.min(), resid.max(), 200)
    ax.plot(xs, stats.norm.pdf(xs, resid.mean(), resid.std(ddof=1)), color="red", linewidth=2)
    ax.set_xlabel("Residuals")
    ax.set_ylabel("Density")
    ax.set_title("Residual Distribution")
    ax.grid(alpha=0.3)

    plt.suptitle(
        f"{subtype} | {virus_name} | Day {day_val} | {vacc_status} [{model_type}]",
        fontsize=12,
        fontweight="bold",
    )
    plt.tight_layout()

    subtype_dir = os.path.join(out_dir, safe_name(subtype))
    group_dir = os.path.join(
        subtype_dir,
        f"{safe_name(virus_name)}__Day{int(day_val)}__{safe_name(vacc_status)}"
    )
    os.makedirs(group_dir, exist_ok=True)

    fig_path = os.path.join(group_dir, "diagnostics.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return fig_path

perf_rows = []

all_fits = [(key, f, "LMM") for key, f in lmm_results.items()] + \
           [(key, f, "OLS") for key, f in ols_results.items()]

for (subtype_name, virus_name, day_val), fit, model_type in all_fits:
    vacc_status = vacc_label.get((subtype_name, virus_name, day_val), "Unknown")

    resid, fitted = get_resid_fitted(fit)
    n_obs = int(fit.nobs) if hasattr(fit, "nobs") else len(resid)
    n_params = len(fit.params)

    rmse = float(np.sqrt(np.mean(resid ** 2)))
    mae = float(np.mean(np.abs(resid)))
    aic = getattr(fit, "aic", np.nan)
    bic = getattr(fit, "bic", np.nan)

    if 3 <= n_obs <= 5000:
        try:
            sw_stat, sw_p = stats.shapiro(resid)
        except Exception:
            sw_stat, sw_p = np.nan, np.nan
    else:
        sw_stat, sw_p = np.nan, np.nan

    marg_r2, cond_r2 = r2_values.get((subtype_name, virus_name, day_val), (np.nan, np.nan))

    fig_path = diagnostic_plots(
        subtype=subtype_name,
        virus_name=virus_name,
        day_val=day_val,
        vacc_status=vacc_status,
        fit=fit,
        model_type=model_type,
        out_dir=DIAG_DIR,
    )

    perf_rows.append({
        "Subtype": subtype_name,
        "Virus": virus_name,
        "Day": day_val,
        "Vaccinated": vacc_status,
        "Model": model_type,
        "N_Obs": n_obs,
        "N_Params": n_params,
        "AIC": round(aic, 2) if pd.notna(aic) else np.nan,
        "BIC": round(bic, 2) if pd.notna(bic) else np.nan,
        "RMSE": round(rmse, 4),
        "MAE": round(mae, 4),
        "Marginal_R2": round(marg_r2, 4) if pd.notna(marg_r2) else np.nan,
        "Conditional_R2": round(cond_r2, 4) if pd.notna(cond_r2) else np.nan,
        "Shapiro_W": round(sw_stat, 4) if pd.notna(sw_stat) else np.nan,
        "Shapiro_p": round(sw_p, 4) if pd.notna(sw_p) else np.nan,
        "Residuals_Normal_at_0.05": (sw_p > 0.05) if pd.notna(sw_p) else np.nan,
        "Diagnostic_Plot": fig_path,
    })

perf_df = pd.DataFrame(perf_rows).sort_values(["Subtype", "Model", "AIC"]).reset_index(drop=True)
perf_df.to_csv(os.path.join(OUTPUT_ROOT, "model_performance_summary_by_subtype.csv"), index=False)



"""
build_html_report.py
=====================
Assembles ONE self-contained HTML page from everything the pipeline script
produces, in this order:

  1. Participant / Day Overlap (UpSet plots)
  2. Demographic / HAI Distribution Plots (faceted, per subtype -> strain)
  3. Phenotype Value Counts
  4. Strain-Level Counts (N samples / studies / cohorts / days)
  5. Regression Overview Plots
       - Marginal vs Conditional R^2 (01_r2_marginal_vs_conditional_*)
       - Cohort contribution gap (02_r2_cohort_contribution_gap_*)
       - Variance decomposition stacked bars (03_variance_decomposition_*)
       - Residual variance ranking (04_residual_variance_ranking_*)
       - Std Dev vs Residual Variance bubble plots (05_stddev_vs_residual_*)
       - Residual (%) vs Std Dev scatter (03_residual(relative)_vs_std_dev_*)
  6. Regression Results (Coefficients) -- estimate / SE / p / CI per term
  7. Model Choice Log (LMM vs OLS vs skipped, per subtype/strain/day)
  8. Variance Decomposition Table (cohort vs residual variance, % split)
  9. Descriptive Plots by Subtype -> Strain -> Day (summary_plots.png +
     summary_stats.csv from descriptive_by_subtype_strain_day/)
  10. Model Performance & Diagnostics (QQ / residual plots + fit metrics)

Images are embedded as base64 so the resulting .html file is fully
portable -- no external image files needed.

USAGE
-----
Run this after the full pipeline has executed in the same session, so
these are already in memory: OUTPUT_ROOT, DIAG_DIR, final_merge, perf_df,
lmm_results, ols_results, vacc_label, choice_df, vdf, strain_counts.
Then call:

    generate_master_html_report()

It writes: OUTPUT_ROOT/master_report.html

Anything not in memory is re-loaded from the CSVs the pipeline already
saves to OUTPUT_ROOT (strain_level_counts.csv, hai_model_choice_log.csv,
variance_decomposition_table.csv, model_performance_summary_by_subtype.csv),
so this also works run standalone in a later session -- just set
OUTPUT_ROOT and DIAG_DIR first.
"""

import os
import re
import glob
import base64
import html as html_lib
import pandas as pd


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------

def _slugify(s):
    s = str(s).strip()
    s = re.sub(r"[^\w\-]+", "_", s)
    return s


def _anchor(text):
    return _slugify(text).lower()


def _img_to_base64_tag(path, max_width_px=900):
    """Read a PNG from disk and return an <img> tag with base64-embedded data."""
    if not path or not os.path.exists(path):
        return "<p><em>[image not found]</em></p>"
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return (
        f'<img src="data:image/png;base64,{data}" '
        f'style="max-width:{max_width_px}px; width:100%; height:auto; '
        f'border:1px solid #ddd; border-radius:6px; margin:8px 0;" />'
    )


def _load_csv_if_exists(path):
    return pd.read_csv(path) if os.path.exists(path) else None


CSS = """
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       max-width: 1150px; margin: 0 auto; padding: 24px; color: #1a1a1a; line-height: 1.5; }
h1 { border-bottom: 3px solid #4472C4; padding-bottom: 8px; }
h2 { margin-top: 52px; border-bottom: 2px solid #ccc; padding-bottom: 6px; color: #2b4a8b; }
h3 { margin-top: 30px; color: #333; }
h4 { margin-top: 18px; color: #555; }
h5 { margin-top: 10px; color: #666; }
table { border-collapse: collapse; width: 100%; margin: 12px 0; font-size: 13px; }
th, td { border: 1px solid #ddd; padding: 5px 9px; text-align: left; }
th { background: #4472C4; color: white; position: sticky; top: 0; }
tr:nth-child(even) { background: #f7f7f9; }
.toc { background: #f7f7f9; border: 1px solid #ddd; border-radius: 8px; padding: 16px 24px; margin-bottom: 32px; }
.toc ul { margin: 4px 0; }
.section-block { margin-bottom: 40px; }
.day-block { margin: 10px 0 10px 12px; padding: 8px 0 8px 12px; border-left: 3px solid #eee; }
.plot-block { margin: 18px 0; padding: 10px 0 10px 12px; border-left: 3px solid #e5e5e5; }
a { color: #2b4a8b; text-decoration: none; }
a:hover { text-decoration: underline; }
.back-to-top { font-size: 12px; }
.note { color: #777; font-size: 13px; font-style: italic; }
"""


# ------------------------------------------------------------------
# 1. UpSet plots
# ------------------------------------------------------------------

def _build_upset_section(output_root):
    files = sorted(glob.glob(os.path.join(output_root, "upset_participants_by_day_*.png")))
    if not files:
        return "", []

    toc_items = []
    parts = ["<h2 id='upset-plots'>1. Participant / Day Overlap (UpSet plots)</h2>"]
    for f in files:
        strain = os.path.basename(f).replace("upset_participants_by_day_", "").replace(".png", "")
        anchor = _anchor(f"upset-{strain}")
        toc_items.append(f"<li><a href='#{anchor}'>{html_lib.escape(strain)}</a></li>")
        parts.append(f"<h3 id='{anchor}'>{html_lib.escape(strain)}</h3>")
        parts.append(_img_to_base64_tag(f, max_width_px=800))
        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 2. Faceted demographic plots
# ------------------------------------------------------------------

def _build_faceted_demo_section(output_root):
    base_dir = os.path.join(output_root, "final_hai_strain_faceted_demo_data_png")
    if not os.path.isdir(base_dir):
        return "", []

    subtype_dirs = sorted(
        d for d in glob.glob(os.path.join(base_dir, "subtype_*")) if os.path.isdir(d)
    )
    if not subtype_dirs:
        return "", []

    toc_items = []
    parts = ["<h2 id='faceted-demo'>2. Demographic / HAI Distribution Plots</h2>"]
    for sdir in subtype_dirs:
        subtype = os.path.basename(sdir).replace("subtype_", "")
        sub_anchor = _anchor(f"demo-{subtype}")
        strain_files = sorted(glob.glob(os.path.join(sdir, "*.png")))

        toc_items.append(
            f"<li><a href='#{sub_anchor}'>{html_lib.escape(subtype)}</a></li>"
        )

        parts.append(f"<h3 id='{sub_anchor}'>Subtype: {html_lib.escape(subtype)}</h3>")
        for sf in strain_files:
            strain_name = os.path.basename(sf).replace(".png", "")
            parts.append(f"<h4>{html_lib.escape(strain_name)}</h4>")
            parts.append(_img_to_base64_tag(sf, max_width_px=1050))
        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 3. Phenotype value counts
# ------------------------------------------------------------------

def _build_phenotype_section(final_merge):
    if final_merge is None or "Phenotype" not in getattr(final_merge, "columns", []):
        return "", []

    counts = final_merge["Phenotype"].value_counts().reset_index()
    counts.columns = ["Phenotype", "Count"]
    table_html = counts.to_html(index=False, border=0)

    toc_items = ["<li><a href='#phenotype-counts'>Phenotype value counts</a></li>"]
    parts = [
        "<h2 id='phenotype-counts'>3. Phenotype Value Counts</h2>",
        table_html,
        "<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>",
    ]
    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 4. Strain-level counts
# ------------------------------------------------------------------

def _build_strain_counts_section(output_root, strain_counts_df=None):
    df = strain_counts_df
    if df is None:
        df = _load_csv_if_exists(os.path.join(output_root, "strain_level_counts.csv"))
    if df is None or len(df) == 0:
        return "", []

    toc_items = ["<li><a href='#strain-counts'>Strain-level counts</a></li>"]
    parts = [
        "<h2 id='strain-counts'>4. Strain-Level Counts</h2>",
        "<p class='note'>N samples / studies / cohorts / days, computed across all "
        "days/visits for each (subtype, virus).</p>",
        df.to_html(index=False, border=0),
        "<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>",
    ]
    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 5. Regression overview plots
# ------------------------------------------------------------------

_OVERVIEW_PLOT_GROUPS = [
    ("01_r2_marginal_vs_conditional_", "Marginal vs Conditional R\u00b2"),
    ("02_r2_cohort_contribution_gap_", "Cohort Contribution to R\u00b2 (gap plot)"),
    ("03_variance_decomposition_", "Variance Decomposition (stacked bars)"),
    ("04_residual_variance_ranking_", "Residual Variance Ranking"),
    ("05_stddev_vs_residual_", "Std Dev vs Residual Variance (bubble plot)"),
    ("03_residual(relative)_vs_std_dev_", "Residual Variance (%) vs Std Dev"),
]


def _build_overview_plots_section(output_root):
    found_any = False
    toc_items = []
    parts = ["<h2 id='overview-plots'>5. Regression Overview Plots</h2>"]

    for prefix, title in _OVERVIEW_PLOT_GROUPS:
        files = sorted(glob.glob(os.path.join(output_root, f"{prefix}*.png")))
        if not files:
            continue
        found_any = True
        group_anchor = _anchor(f"overview-{prefix}")
        toc_items.append(f"<li><a href='#{group_anchor}'>{html_lib.escape(title)}</a></li>")
        parts.append(f"<h3 id='{group_anchor}'>{html_lib.escape(title)}</h3>")

        for f in files:
            label = os.path.basename(f).replace(prefix, "").replace(".png", "")
            label = label.strip("_")
            parts.append(f"<h4>{html_lib.escape(label) if label else 'All groups'}</h4>")
            parts.append(_img_to_base64_tag(f, max_width_px=1050))

        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    if not found_any:
        return "", []
    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 6. Regression coefficients
# ------------------------------------------------------------------

def _fit_coef_table(fit):
    """Pull term/estimate/SE/p-value/CI out of a statsmodels fit object."""
    try:
        params = fit.params
        bse = fit.bse
        pvals = fit.pvalues
        terms = list(params.index)

        rows = {
            "Term": terms,
            "Estimate": [round(float(v), 4) for v in params.values],
            "Std Error": [round(float(v), 4) for v in bse.values],
            "p-value": [round(float(v), 4) for v in pvals.values],
        }

        try:
            conf = fit.conf_int()
            conf = conf.values if hasattr(conf, "values") else conf
            rows["CI Lower"] = [round(float(v), 4) for v in conf[:, 0]]
            rows["CI Upper"] = [round(float(v), 4) for v in conf[:, 1]]
        except Exception:
            pass

        df = pd.DataFrame(rows)
        df["Sig."] = df["p-value"].apply(
            lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        )
        return df
    except Exception:
        return None


def _build_regression_results_section_from_fits(lmm_results, ols_results, vacc_label):
    """Preferred path: build coefficient tables directly from the fit objects
    (includes cohort/residual variance rows for LMMs)."""
    if not lmm_results and not ols_results:
        return None

    all_fits = [(key, f, "LMM") for key, f in (lmm_results or {}).items()] + \
               [(key, f, "OLS") for key, f in (ols_results or {}).items()]
    if not all_fits:
        return None

    by_subtype = {}
    for (subtype_name, virus_name, day_val), fit, model_type in all_fits:
        by_subtype.setdefault(subtype_name, []).append(
            (virus_name, day_val, fit, model_type)
        )

    toc_items = []
    parts = ["<h2 id='regression-results'>6. Regression Results (Coefficients)</h2>"]

    for subtype_name in sorted(by_subtype.keys(), key=str):
        sub_anchor = _anchor(f"coef-{subtype_name}")
        toc_items.append(f"<li><a href='#{sub_anchor}'>{html_lib.escape(str(subtype_name))}</a></li>")
        parts.append(f"<h3 id='{sub_anchor}'>Subtype: {html_lib.escape(str(subtype_name))}</h3>")

        entries = sorted(by_subtype[subtype_name], key=lambda e: (str(e[0]), e[1], e[3]))
        for virus_name, day_val, fit, model_type in entries:
            vacc_status = (vacc_label or {}).get((subtype_name, virus_name, day_val), "Unknown")
            coef_df = _fit_coef_table(fit)

            parts.append(
                f"<div class='day-block'><b>{html_lib.escape(str(virus_name))} &mdash; "
                f"Day {day_val} &mdash; {html_lib.escape(str(vacc_status))} &mdash; {model_type}</b>"
            )
            if coef_df is not None:
                parts.append(coef_df.to_html(index=False, border=0))
            else:
                parts.append("<p><em>[could not extract coefficient table for this fit]</em></p>")
            parts.append("</div>")

        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    return "\n".join(parts), toc_items


def _build_regression_results_section_from_csv(output_root):
    """Fallback path: read hai_regression_summary.csv (Subtype/Virus/Day/
    Vaccinated/Model/Term/Estimate/CI_low/CI_high/p_value) if fit objects
    aren't available in memory."""
    df = _load_csv_if_exists(os.path.join(output_root, "hai_regression_summary.csv"))
    if df is None or len(df) == 0:
        return "", []

    toc_items = []
    parts = ["<h2 id='regression-results'>6. Regression Results (Coefficients)</h2>"]

    for subtype_name, sdf in df.groupby("Subtype"):
        sub_anchor = _anchor(f"coef-{subtype_name}")
        toc_items.append(f"<li><a href='#{sub_anchor}'>{html_lib.escape(str(subtype_name))}</a></li>")
        parts.append(f"<h3 id='{sub_anchor}'>Subtype: {html_lib.escape(str(subtype_name))}</h3>")

        for (virus_name, day_val, vacc_status, model_type), gdf in sdf.groupby(
            ["Virus", "Day", "Vaccinated", "Model"]
        ):
            parts.append(
                f"<div class='day-block'><b>{html_lib.escape(str(virus_name))} &mdash; "
                f"Day {day_val} &mdash; {html_lib.escape(str(vacc_status))} &mdash; {model_type}</b>"
            )
            cols = [c for c in ["Term", "Estimate", "CI_low", "CI_high", "p_value"] if c in gdf.columns]
            parts.append(gdf[cols].to_html(index=False, border=0))
            parts.append("</div>")

        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    return "\n".join(parts), toc_items


def _build_regression_results_section(output_root, lmm_results, ols_results, vacc_label):
    result = _build_regression_results_section_from_fits(lmm_results, ols_results, vacc_label)
    if result is not None:
        return result
    return _build_regression_results_section_from_csv(output_root)


# ------------------------------------------------------------------
# 7. Model choice log
# ------------------------------------------------------------------

def _build_model_choice_section(output_root, choice_df_in=None):
    df = choice_df_in
    if df is None:
        df = _load_csv_if_exists(os.path.join(output_root, "hai_model_choice_log.csv"))
    if df is None or len(df) == 0:
        return "", []

    toc_items = ["<li><a href='#model-choice'>Model choice log</a></li>"]
    parts = [
        "<h2 id='model-choice'>7. Model Choice Log</h2>",
        "<p class='note'>Which model (LMM vs OLS) was fit for each subtype/strain/day, "
        "or why a group was skipped (N &lt; 20).</p>",
        df.to_html(index=False, border=0),
        "<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>",
    ]
    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 8. Variance decomposition table
# ------------------------------------------------------------------

def _build_variance_decomposition_table_section(output_root, vdf_in=None):
    df = vdf_in
    if df is None:
        df = _load_csv_if_exists(os.path.join(output_root, "variance_decomposition_table.csv"))
    if df is None or len(df) == 0:
        return "", []

    preferred_cols = [
        "Subtype", "Virus", "Day", "Vaccinated",
        "Cohort_Variance", "Residual_Variance", "Total_Variance",
        "Cohort_Variance_pct", "Residual_Variance_pct",
        "Std_HAI_log2", "N_Participants", "N_Samples", "N_Studies", "N_Cohorts", "N_Days",
    ]
    cols = [c for c in preferred_cols if c in df.columns] + \
           [c for c in df.columns if c not in preferred_cols and c not in ("Group_short", "Group_label_full")]

    toc_items = ["<li><a href='#variance-table'>Variance decomposition table</a></li>"]
    parts = [
        "<h2 id='variance-table'>8. Variance Decomposition Table</h2>",
        df[cols].to_html(index=False, border=0),
        "<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>",
    ]
    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 9. Descriptive plots by subtype -> strain -> day
# ------------------------------------------------------------------

def _build_descriptive_plots_section(output_root):
    base_dir = os.path.join(output_root, "descriptive_by_subtype_strain_day")
    if not os.path.isdir(base_dir):
        return "", []

    plot_files = sorted(glob.glob(os.path.join(base_dir, "*", "*", "*", "summary_plots.png")))
    if not plot_files:
        return "", []

    toc_items = []
    parts = ["<h2 id='descriptive-plots'>9. Descriptive Plots by Subtype &rarr; Strain &rarr; Day</h2>"]

    by_subtype = {}
    for f in plot_files:
        day_dir = os.path.dirname(f)
        virus_dir = os.path.dirname(day_dir)
        subtype_dir = os.path.dirname(virus_dir)
        subtype = os.path.basename(subtype_dir)
        virus = os.path.basename(virus_dir)
        day_label = os.path.basename(day_dir)
        by_subtype.setdefault(subtype, {}).setdefault(virus, []).append((day_label, f))

    for subtype in sorted(by_subtype.keys()):
        sub_anchor = _anchor(f"desc-{subtype}")
        toc_items.append(f"<li><a href='#{sub_anchor}'>{html_lib.escape(subtype)}</a></li>")
        parts.append(f"<h3 id='{sub_anchor}'>Subtype: {html_lib.escape(subtype)}</h3>")

        for virus in sorted(by_subtype[subtype].keys()):
            parts.append(f"<h4>{html_lib.escape(virus)}</h4>")
            for day_label, fpath in sorted(by_subtype[subtype][virus]):
                stats_path = os.path.join(os.path.dirname(fpath), "summary_stats.csv")
                parts.append(f"<div class='plot-block'><h5>{html_lib.escape(day_label)}</h5>")
                parts.append(_img_to_base64_tag(fpath, max_width_px=900))
                if os.path.exists(stats_path):
                    stats_df = pd.read_csv(stats_path)
                    stats_df.columns = ["Statistic", "log2_HAI"] if len(stats_df.columns) == 2 else stats_df.columns
                    parts.append(stats_df.to_html(index=False, border=0))
                parts.append("</div>")

        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# 10. Model performance & diagnostics
# ------------------------------------------------------------------

def _build_diagnostics_section(perf_df, diag_dir):
    if perf_df is None or len(perf_df) == 0:
        return "", []

    toc_items = []
    parts = ["<h2 id='diagnostics'>10. Model Performance &amp; Diagnostics</h2>"]

    summary_cols = [c for c in perf_df.columns if c != "Diagnostic_Plot"]
    parts.append("<h3>Full performance summary table</h3>")
    parts.append(perf_df[summary_cols].to_html(index=False, border=0))
    parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    for subtype, sdf in perf_df.groupby("Subtype"):
        sub_anchor = _anchor(f"diag-{subtype}")
        toc_items.append(f"<li><a href='#{sub_anchor}'>{html_lib.escape(str(subtype))}</a></li>")
        parts.append(f"<h3 id='{sub_anchor}'>Subtype: {html_lib.escape(str(subtype))}</h3>")

        for virus, vdf_g in sdf.groupby("Virus"):
            parts.append(f"<h4>{html_lib.escape(str(virus))}</h4>")

            for _, row in vdf_g.sort_values(["Day", "Vaccinated", "Model"]).iterrows():
                parts.append(
                    f"<div class='day-block'><b>Day {row['Day']} &mdash; "
                    f"{html_lib.escape(str(row['Vaccinated']))} &mdash; {row['Model']}</b><br>"
                    f"N={row['N_Obs']}, AIC={row['AIC']}, BIC={row['BIC']}, "
                    f"RMSE={row['RMSE']}, Marginal R&sup2;={row['Marginal_R2']}, "
                    f"Conditional R&sup2;={row['Conditional_R2']}, "
                    f"Shapiro p={row['Shapiro_p']}"
                    f"</div>"
                )
                parts.append(_img_to_base64_tag(row.get("Diagnostic_Plot"), max_width_px=950))

        parts.append("<p class='back-to-top'><a href='#toc'>&uarr; back to top</a></p>")

    return "\n".join(parts), toc_items


# ------------------------------------------------------------------
# main entry point
# ------------------------------------------------------------------

def generate_master_html_report(
    output_root=None,
    diag_dir=None,
    final_merge_df=None,
    perf_df_in=None,
    lmm_results_in=None,
    ols_results_in=None,
    vacc_label_in=None,
    choice_df_in=None,
    vdf_in=None,
    strain_counts_in=None,
    out_filename="master_report.html",
):
    """
    Build one self-contained HTML page combining every plot/table the
    pipeline produced. Anything not passed explicitly is pulled from the
    current global namespace (OUTPUT_ROOT, DIAG_DIR, final_merge, perf_df,
    lmm_results, ols_results, vacc_label, choice_df, vdf, strain_counts),
    or re-loaded from the CSVs the pipeline already writes to OUTPUT_ROOT.
    """
    g = globals()
    output_root = output_root or g.get("OUTPUT_ROOT")
    diag_dir = diag_dir or g.get("DIAG_DIR")
    final_merge_df = final_merge_df if final_merge_df is not None else g.get("final_merge")
    perf_df_in = perf_df_in if perf_df_in is not None else g.get("perf_df")
    lmm_results_in = lmm_results_in if lmm_results_in is not None else g.get("lmm_results")
    ols_results_in = ols_results_in if ols_results_in is not None else g.get("ols_results")
    vacc_label_in = vacc_label_in if vacc_label_in is not None else g.get("vacc_label")
    choice_df_in = choice_df_in if choice_df_in is not None else g.get("choice_df")
    vdf_in = vdf_in if vdf_in is not None else g.get("vdf")
    strain_counts_in = strain_counts_in if strain_counts_in is not None else g.get("strain_counts")

    if output_root is None:
        raise ValueError("output_root not found -- pass it explicitly (e.g. OUTPUT_ROOT).")

    if perf_df_in is None:
        csv_path = os.path.join(output_root, "model_performance_summary_by_subtype.csv")
        if os.path.exists(csv_path):
            perf_df_in = pd.read_csv(csv_path)

    upset_html, upset_toc = _build_upset_section(output_root)
    demo_html, demo_toc = _build_faceted_demo_section(output_root)
    pheno_html, pheno_toc = _build_phenotype_section(final_merge_df)
    strain_counts_html, strain_counts_toc = _build_strain_counts_section(output_root, strain_counts_in)
    overview_html, overview_toc = _build_overview_plots_section(output_root)
    coef_html, coef_toc = _build_regression_results_section(
        output_root, lmm_results_in, ols_results_in, vacc_label_in
    )
    choice_html, choice_toc = _build_model_choice_section(output_root, choice_df_in)
    vartable_html, vartable_toc = _build_variance_decomposition_table_section(output_root, vdf_in)
    desc_html, desc_toc = _build_descriptive_plots_section(output_root)
    diag_html, diag_toc = _build_diagnostics_section(perf_df_in, diag_dir)

    toc_html = f"""
    <div class="toc" id="toc">
      <h2 style="margin-top:0;">Table of Contents</h2>
      <ol>
        <li><a href="#upset-plots">Participant / Day Overlap (UpSet plots)</a><ul>{''.join(upset_toc)}</ul></li>
        <li><a href="#faceted-demo">Demographic / HAI Distribution Plots</a><ul>{''.join(demo_toc)}</ul></li>
        <li><a href="#phenotype-counts">Phenotype Value Counts</a></li>
        <li><a href="#strain-counts">Strain-Level Counts</a></li>
        <li><a href="#overview-plots">Regression Overview Plots</a><ul>{''.join(overview_toc)}</ul></li>
        <li><a href="#regression-results">Regression Results (Coefficients)</a><ul>{''.join(coef_toc)}</ul></li>
        <li><a href="#model-choice">Model Choice Log</a></li>
        <li><a href="#variance-table">Variance Decomposition Table</a></li>
        <li><a href="#descriptive-plots">Descriptive Plots by Subtype &rarr; Strain &rarr; Day</a><ul>{''.join(desc_toc)}</ul></li>
        <li><a href="#diagnostics">Model Performance &amp; Diagnostics</a><ul>{''.join(diag_toc)}</ul></li>
      </ol>
    </div>
    """

    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>HAI Pipeline -- Master Report</title>
<style>{CSS}</style>
</head>
<body>
<h1>HAI Pipeline &mdash; Master Report</h1>
<p><em>Auto-generated report combining all pipeline outputs.</em></p>
{toc_html}
{upset_html}
{demo_html}
{pheno_html}
{strain_counts_html}
{overview_html}
{coef_html}
{choice_html}
{vartable_html}
{desc_html}
{diag_html}
</body>
</html>"""

    out_path = os.path.join(output_root, out_filename)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(full_html)

    print(f"Master HTML report written to: {out_path}")
    return out_path


if __name__ == "__main__":
    generate_master_html_report()























##post confirmation regression analysis 

#1 filter final_merge for final strain list and only d0 and d28 for a single cohort 
strain_list = ['A/Solomon Islands/3/2006', 'A/California/7/2009', 'A/Perth/16/2009',
               'A/Victoria/361/2011', 'B/Malaysia/2506/2004', 'B/Brisbane/3/2007', 'B/Wisconsin/01/2010']

post_df = final_merge[(final_merge['Virus'].isin(strain_list)) & 
                      (final_merge['Study Time Collected'].isin([0, 28]))]


#h1n1
h1n1 = post_df[post_df['Virus'].isin(['A/California/7/2009', 'A/Solomon Islands/3/2006'])]
solomon = h1n1[(h1n1["Virus"] == 'A/Solomon Islands/3/2006') & (h1n1['Cohort for regression'] == 'Older participants aged 60 to 89 years, vaccinated with Fluzone')]
cali = h1n1[(h1n1["Virus"] == 'A/California/7/2009') &   (h1n1['Cohort for regression'] == '150 healthy adults, 50-74 yo')]


#h3n2
h3n2 = post_df[(post_df['Virus'].isin(['A/Perth/16/2009', 'A/Victoria/361/2011']))]
perth = h3n2[(h3n2['Virus'] == 'A/Perth/16/2009') & (h3n2['Cohort for regression'] == '150 healthy adults, 50-74 yo')]
vic = h3n2[(h3n2['Virus'] == 'A/Victoria/361/2011') & (h3n2['Cohort for regression'] == 'Healthy Adults 2012 - 2013')]


#victoria
victoria = post_df[post_df['Virus'] == "B/Malaysia/2506/2004"] 
malaysia = victoria[victoria['Cohort for regression'] == 'Older participants aged 60 to 89 years, vaccinated with Fluzone']



#yamagata
yamagata = post_df[post_df['Virus'].isin(['B/Brisbane/3/2007', 'B/Wisconsin/01/2010'])]
brisbane = yamagata[(yamagata['Virus'] == 'B/Brisbane/3/2007')] #single cohort for regression already #add data from publication for confirmation
wisconsin = yamagata[(yamagata['Virus'] == 'B/Wisconsin/01/2010') & (yamagata['Cohort for regression'] == 'Healthy Adults 2012 - 2013')]


 





# Put all 6 DataFrames in a list
replication_df_list = [solomon, cali, perth, vic, malaysia, brisbane, wisconsin]

# Stacking them into 1 single DataFrame
rep_df = pd.concat(replication_df_list, ignore_index=True)





OUTPUT_ROOT2 = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/immunespace_replication_hai_regression_results"
os.makedirs(OUTPUT_ROOT2, exist_ok=True)

#
"""
UpSet plots: Participants x Days, computed separately per strain.

Sets     = Study Time Collected (Day 0, Day 28, ...)
Elements = participants (new_participant_id)
Split by = strain (subtype, or Virus for finer granularity)

For each strain, this shows how many participants were sampled at each day,
and which combinations of days each participant was sampled at -- so you can
compare e.g. "H1N1 day-overlap" vs "H3N2 day-overlap".


"""

import pandas as pd
from upsetplot import from_memberships, UpSet
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Load your data here. Replace this with your actual source (csv, etc.)


# df = pd.read_csv("your_file.csv")
# ------------------------------------------------------------------
df2 = rep_df.copy()  # assumes `df` is already the dataframe shown in your message

# use `subtype` (H1N1/H3N2/Yamagata/...) for the strain grouping.
# swap to "Virus" instead if you want individual reference-strain granularity
# (e.g. A/South Dakota/06/2007 vs A/Uruguay/716/2007) rather than subtype.
strain_col = "Virus"

# one row per (participant, day, strain), regardless of everything else
part_day_strain = (
    df2[["Participant ID_merge1", "Study Time Collected", strain_col]]
    .dropna(subset=["Study Time Collected", strain_col])
    .drop_duplicates()
)

strains = sorted(part_day_strain[strain_col].unique())

for strain in strains:
    subset = part_day_strain[part_day_strain[strain_col] == strain]

    # for each participant, which set of days were they sampled at for this strain?
    memberships = (
        subset.groupby("Participant ID_merge1")["Study Time Collected"]
        .apply(lambda s: [f"Day {int(d)}" for d in sorted(s.unique())])
        .tolist()
    )

    upset_data = from_memberships(memberships)

    fig = plt.figure(figsize=(10, 6))
    UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
    ).plot(fig=fig)
    plt.suptitle(f"Participant overlap across collection days -- {strain}")

    safe_name = str(strain).replace("/", "-").replace(" ", "_")
    plt.savefig(OUTPUT_ROOT2 + f"/replication_upset_participants_by_day_{safe_name}.png",
        dpi=150,
        bbox_inches="tight",
    )

plt.show()



"""
HAI regression + variance analysis, organized by Subtype -> Strain -> Day
===========================================================================
Builds on the original per-strain-per-day LMM/OLS logic, but:

  1. Everything is grouped/output by SUBTYPE -> STRAIN -> DAY, not just
     strain-day, so results for the same subtype sit together.
  2. Every strain gets an annotation block: total N samples, N distinct
     studies, N distinct cohorts (computed once, across ALL days/visits
     for that strain).
  3. New "Std Dev vs Residual Variance" bubble plot, with the N for each
     group written directly next to its bubble (not just encoded in size).
  4. All descriptive + model-fit + variance plots for a subtype/strain/day
     are collected into ONE self-contained HTML report you can open from
     a single link, with a jump-to table of contents.

Run this after `final_merge` is already in memory (same as your original
script).
"""






# ---------------------------------------------------------------
# 1. LOAD + CLEAN (same as original)
# ---------------------------------------------------------------
df2 = rep_df.copy()
df2.columns = df2.columns.str.strip()

COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_hai",
    "sex": "Gender_hai",
    "cohort": "Cohort for regression",
    "virus": "Virus",
    "subtype": "subtype",
    "study_id": "Study_ID",
    "participant_id": "Participant ID",
    "day": "Study Time Collected",
}

clean = df2.rename(columns={
    COLUMNS["outcome"]: "log2_HAI",
    COLUMNS["age"]: "age",
    COLUMNS["sex"]: "sex",
    COLUMNS["cohort"]: "Cohort for regression",
    COLUMNS["virus"]: "virus",
    COLUMNS["subtype"]: "subtype",
    COLUMNS["study_id"]: "study_id",
    COLUMNS["participant_id"]: "participant_id",
    COLUMNS["day"]: "day",
})

needed = ["log2_HAI", "age", "sex", "Cohort for regression", "virus", "subtype", "study_id", "day"]
missing = [c for c in needed if c not in clean.columns]
if missing:
    raise KeyError(f"Rename didn't produce expected columns: {missing}. "
                    f"Available columns are: {df.columns.tolist()}")

clean = clean.dropna(subset=needed).copy()
clean["sex"] = clean["sex"].astype("category")
clean["cohort"] = clean["Cohort for regression"].astype("category")
clean["subtype"] = clean["subtype"].astype("category")
clean["age"] = pd.to_numeric(clean["age"], errors="coerce")
clean["day"] = pd.to_numeric(clean["day"], errors="coerce")
clean = clean.dropna(subset=["age", "day"])

if "Study Time Collected Unit" in df2.columns:
    bad_units = df2.loc[clean.index, "Study Time Collected Unit"].dropna().unique()
    bad_units = [u for u in bad_units if str(u).strip().lower() != "days"]
    if bad_units:
        raise ValueError(
            f"'Study Time Collected' isn't uniformly in Days -- found units: {bad_units}."
        )

clean["vaccinated"] = np.where(clean["day"] >= 14, "Vaccinated", "Not Vaccinated")
clean["vaccinated"] = clean["vaccinated"].astype("category")


def mixedlm_r2(fit):
    fixed_effects_fitted = fit.model.exog @ fit.fe_params
    var_fixed = np.var(fixed_effects_fitted, ddof=0)
    var_random = fit.cov_re.iloc[0, 0]
    var_resid = fit.scale
    total = var_fixed + var_random + var_resid
    return var_fixed / total, (var_fixed + var_random) / total


def safe_name(s):
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s


# ---------------------------------------------------------------
# 2. STRAIN-LEVEL ANNOTATION COUNTS
#    (computed ACROSS ALL DAYS for that strain -- this is the
#    "how many samples / studies / cohorts" annotation)
# ---------------------------------------------------------------
strain_counts = (
    clean.groupby(["subtype", "virus"])
    .agg(
        N_Samples=("log2_HAI", "count"),
        N_Studies=("study_id", "nunique"),
        N_Cohorts=("Cohort for regression", "nunique"),
        N_Days=("day", "nunique"),
    )
    .reset_index()
)
strain_counts.to_csv(os.path.join(OUTPUT_ROOT2, "replication_strain_level_counts.csv"), index=False)
print("Strain-level counts (samples / studies / cohorts):")
print(strain_counts.to_string(index=False))


def strain_annotation_text(virus_name):
    row = strain_counts[strain_counts["virus"] == virus_name]
    if row.empty:
        return "N_Samples=?, N_Studies=?, N_Cohorts=?"
    r = row.iloc[0]
    return f"N_Samples={r['N_Samples']} | N_Studies={r['N_Studies']} | N_Cohorts={r['N_Cohorts']} | N_Days={r['N_Days']}"

# ---------------------------------------------------------------
# 3. PER (subtype, virus, day) MODEL FIT -- OLS only
# ---------------------------------------------------------------
ols_results, r2_values, vacc_label = {}, {}, {}

for (subtype_name, virus_name, day_val), sub in clean.groupby(["subtype", "virus", "day"]):
    sub = sub.copy()
    n_obs = len(sub)
    vacc_status = sub["vaccinated"].iloc[0]   # label only, not modeled
    key = (subtype_name, virus_name, day_val)
    vacc_label[key] = vacc_status

    if n_obs < 19:
        continue

    try:
        ols_fit = smf.ols("log2_HAI ~ age + sex", data=sub).fit()
        ols_results[key] = ols_fit
        r2_values[key] = ols_fit.rsquared
    except Exception as e:
        print(f"  -> OLS failed for {key}: {e}")

# ---------------------------------------------------------------
# 4. UNIFIED SUMMARY TABLE (OLS only)
# ---------------------------------------------------------------
rows = []
for (subtype_name, virus_name, day_val), fit in ols_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()

    for term in params.index:
        rows.append({
            "Subtype": subtype_name,
            "Virus": virus_name,
            "Day": day_val,
            "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
            "Model": "OLS",
            "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4),
            "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })

    rows.append({
        "Subtype": subtype_name,
        "Virus": virus_name,
        "Day": day_val,
        "Vaccinated": vacc_label[(subtype_name, virus_name, day_val)],
        "Model": "OLS",
        "Term": "Residual variance",
        "Estimate": round(fit.mse_resid, 4),
        "CI_low": np.nan,
        "CI_high": np.nan,
        "p_value": np.nan,
    })

summary_df = pd.DataFrame(rows)
summary_df.to_csv(os.path.join(OUTPUT_ROOT2, "replication_hai_regression_summary.csv"), index=False)

# ---------------------------------------------------------------
# 5. MODEL CHOICE LOG
# ---------------------------------------------------------------
choice_rows = []
for (subtype_name, virus_name, day_val), sub in clean.groupby(["subtype", "virus", "day"]):
    n_obs = len(sub)
    key = (subtype_name, virus_name, day_val)
    r2 = r2_values.get(key, np.nan)

    choice_rows.append({
        "Subtype": subtype_name,
        "Virus": virus_name,
        "Day": day_val,
        "Vaccinated": vacc_label.get(key, sub["vaccinated"].iloc[0]),
        "N_Obs": n_obs,
        "N_Cohorts": 1,
        "Model_Used": "OLS" if key in ols_results else "Skipped (n_obs<19)",
        "R2": round(r2, 4) if pd.notna(r2) else np.nan,
        "Residual_Unexplained_pct": round(100 * (1 - r2), 4) if pd.notna(r2) else np.nan,
    })

choice_df = pd.DataFrame(choice_rows)
choice_df.to_csv(os.path.join(OUTPUT_ROOT2, "replication_hai_model_choice_log.csv"), index=False)

# ---------------------------------------------------------------
# 6. VARIANCE / EXPLAINED-VARIANCE TABLE (OLS only)
# ---------------------------------------------------------------
resid_rows = summary_df[summary_df["Term"] == "Residual variance"].copy()
resid_rows = resid_rows.rename(columns={"Estimate": "Residual_Variance"})

vdf = choice_df.merge(
    resid_rows[["Subtype", "Virus", "Day", "Residual_Variance"]],
    on=["Subtype", "Virus", "Day"],
    how="left",
)

stats_by_group = (
    clean.groupby(["subtype", "virus", "day"])["log2_HAI"]
    .agg(Std_HAI_log2="std", N_Participants="count")
    .reset_index()
    .rename(columns={"subtype": "Subtype", "virus": "Virus", "day": "Day"})
)

vdf = vdf.merge(stats_by_group, on=["Subtype", "Virus", "Day"], how="left")
vdf = vdf.merge(
    strain_counts.rename(columns={"virus": "Virus", "subtype": "Subtype"}),
    on=["Subtype", "Virus"],
    how="left",
)

# OLS interpretation:
# explained variance = R²
# unexplained variance = 1 - R²
vdf["Explained_Variance_pct"] = 100 * vdf["R2"]
vdf["Residual_Unexplained_pct"] = 100 * (1 - vdf["R2"])

def short_label(virus, day, maxlen=28):
    v = str(virus).replace("Influenza ", "")[:maxlen]
    return f"{v} (D{int(day)})"

vdf["Group_short"] = vdf.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)
vdf["Group_label_full"] = vdf.apply(
    lambda r: f"{short_label(r['Virus'], r['Day'])} [{r['Subtype']}]",
    axis=1
)

vdf.to_csv(os.path.join(OUTPUT_ROOT2, "replication_variance_decomposition_table.csv"), index=False)

overview_plot_paths = []

# =================================================================
# 7a. PLOT: OLS R² per group -- one bar per strain/day
# =================================================================
df_plot = choice_df.dropna(subset=["Model_Used"]).copy()
df_plot = df_plot[df_plot["Model_Used"] != "Skipped (n_obs<20)"].copy()

df_plot["Group_short"] = df_plot.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)
df_plot["Group_label_full"] = df_plot.apply(
    lambda r: f"{short_label(r['Virus'], r['Day'])} [{r['Subtype']}]",
    axis=1
)

for subtype in sorted(df_plot["Subtype"].dropna().unique()):
    sub = (
        df_plot[df_plot["Subtype"] == subtype]
        .sort_values("R2", ascending=False)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(max(8, 1.8 * len(sub)), 7))
    x = np.arange(len(sub))

    bars = ax.bar(
        x,
        sub["R2"],
        width=0.6,
        color="#4472C4",
        edgecolor="black",
        linewidth=0.6,
    )

    for b in bars:
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.02,
            f"{b.get_height():.3f}",
            ha="center",
            fontsize=9,
        )

    for i, vacc in enumerate(sub["Vaccinated"]):
        ax.text(
            i,
            -0.06,
            vacc,
            ha="center",
            va="top",
            fontsize=7.5,
            color="gray",
            style="italic",
        )

    ax.set_xlabel("Strain (Day)", fontsize=12, fontweight="bold")
    ax.set_ylabel("R²", fontsize=12, fontweight="bold")
    ax.set_title(
        f"OLS R² by Strain/Day\nSubtype: {subtype}",
        fontsize=13,
        fontweight="bold",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(sub["Group_short"], rotation=45, ha="right", fontsize=8.5)
    ax.set_ylim([0, 1.0])
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT2,
        f"01_replication_r2_{subtype.replace('/', '_').replace(' ', '_')}.png",
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    overview_plot_paths.append((f"OLS R² ({subtype})", outfile))
    print(f"Saved: {outfile}")

# =================================================================
# 7b. PLOT: Residual unexplained variance (%) by group
# =================================================================
for subtype in sorted(vdf["Subtype"].dropna().unique()):
    sub = (
        vdf[vdf["Subtype"] == subtype]
        .sort_values("Residual_Unexplained_pct", ascending=False)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(sub)), 7))
    x = np.arange(len(sub))

    ax.bar(
        x,
        sub["Residual_Unexplained_pct"],
        width=0.6,
        color="#ED7D31",
        edgecolor="black",
        linewidth=0.6,
    )

    for i, row in sub.iterrows():
        ax.text(
            i,
            row["Residual_Unexplained_pct"] + 1,
            f"{row['Residual_Unexplained_pct']:.1f}%",
            ha="center",
            fontsize=9,
            fontweight="bold",
        )

    for i, vacc in enumerate(sub["Vaccinated"]):
        ax.text(
            i,
            -6,
            vacc,
            ha="center",
            fontsize=7,
            color="gray",
            style="italic",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(sub["Group_label_full"], rotation=45, ha="right")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title(f"Residual Unexplained Variance (%)\n{subtype}")

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT2,
        f"02_replication_residual_unexplained_{safe_name(subtype)}.png",
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {outfile}")

# =================================================================
# 7c. PLOT: Explained vs Unexplained variance -- stacked bar
# =================================================================
for subtype in sorted(vdf["Subtype"].unique()):
    sub = (
        vdf[vdf["Subtype"] == subtype]
        .sort_values("Residual_Unexplained_pct", ascending=False)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(sub)), 7))
    x = np.arange(len(sub))

    ax.bar(
        x,
        sub["Explained_Variance_pct"],
        width=0.6,
        color="#4472C4",
        label="Explained variance (R²)",
    )

    ax.bar(
        x,
        sub["Residual_Unexplained_pct"],
        width=0.6,
        bottom=sub["Explained_Variance_pct"],
        color="#ED7D31",
        label="Unexplained variance (1 - R²)",
    )

    for i, vacc in enumerate(sub["Vaccinated"]):
        ax.text(
            i,
            -6,
            vacc,
            ha="center",
            fontsize=7,
            color="gray",
            style="italic",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(sub["Group_label_full"], rotation=45, ha="right")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3)
    ax.legend()
    ax.set_title(f"OLS Variance Split (Explained vs Unexplained)\n{subtype}")

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT2,
        f"03_replication_variance_split_{safe_name(subtype)}.png",
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {outfile}")

# =================================================================
# 7d. PLOT: Ranked horizontal bar -- Residual unexplained variance %
# =================================================================
for subtype in sorted(vdf["Subtype"].unique()):
    sub = vdf[vdf["Subtype"] == subtype]

    top_groups = (
        sub.sort_values("Residual_Unexplained_pct", ascending=False)
        .head(15)
        .sort_values("Residual_Unexplained_pct")
    )

    fig, ax = plt.subplots(figsize=(12, max(5, 0.45 * len(top_groups))))
    y_pos = np.arange(len(top_groups))

    colors = [
        "#d73027" if v > 60
        else "#fee090" if v > 40
        else "#91bfdb"
        for v in top_groups["Residual_Unexplained_pct"]
    ]

    ax.barh(
        y_pos,
        top_groups["Residual_Unexplained_pct"],
        color=colors,
        edgecolor="black",
    )

    for i, (_, row) in enumerate(top_groups.iterrows()):
        ax.text(
            row["Residual_Unexplained_pct"] + 1,
            i,
            f"{row['Residual_Unexplained_pct']:.1f}%",
            va="center",
            fontsize=9,
            fontweight="bold",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_groups["Group_label_full"])
    ax.set_xlim(0, 100)
    ax.axvline(40, ls="--", color="orange")
    ax.axvline(60, ls="--", color="red")
    ax.grid(axis="x", alpha=0.3)
    ax.set_title(f"Residual Unexplained Variance Ranking\n{subtype}")

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT2,
        f"04_replication_residual_unexplained_ranking_{safe_name(subtype)}.png",
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {outfile}")

# =================================================================
# 7e. Std Dev vs Unexplained Variance bubble plot
# =================================================================
import itertools

MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]

def plot_std_vs_unexplained(sub_vdf, title, out_path):
    fig, ax = plt.subplots(figsize=(14, 8))

    strains = sorted(sub_vdf["Virus"].dropna().unique())
    days = sorted(sub_vdf["Day"].dropna().unique())

    strain_colors = plt.cm.tab20(np.linspace(0, 1, max(len(strains), 1)))
    strain_color_map = dict(zip(strains, strain_colors))

    marker_cycle = itertools.cycle(MARKER_CYCLE)
    day_marker_map = {day: next(marker_cycle) for day in days}

    for (virus, day), grp in sub_vdf.groupby(["Virus", "Day"]):
        ax.scatter(
            grp["Std_HAI_log2"],
            grp["Residual_Unexplained_pct"],
            s=grp["N_Participants"] * 6,
            alpha=0.75,
            color=strain_color_map[virus],
            marker=day_marker_map[day],
            edgecolors="black",
            linewidth=0.6,
        )

    for _, row in sub_vdf.iterrows():
        ax.annotate(
            f"{row['Group_short']}\nN={int(row['N_Participants'])}",
            (row["Std_HAI_log2"], row["Residual_Unexplained_pct"]),
            textcoords="offset points",
            xytext=(7, 5),
            fontsize=7,
            color="black",
            linespacing=1.3,
        )

    ax.set_xlabel("HAI Titer Std Dev (log2 scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Unexplained variance (%)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)

    strain_handles = [
        plt.Line2D(
            [0], [0],
            marker="o",
            color="w",
            markerfacecolor=strain_color_map[s],
            markeredgecolor="black",
            markersize=8,
            label=s,
        )
        for s in strains
    ]
    legend1 = ax.legend(
        handles=strain_handles,
        title="Strain",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        fontsize=8,
    )
    ax.add_artist(legend1)

    day_handles = [
        plt.Line2D(
            [0], [0],
            marker=day_marker_map[d],
            color="gray",
            linestyle="None",
            markeredgecolor="black",
            markersize=8,
            label=f"Day {d}",
        )
        for d in days
    ]
    legend2 = ax.legend(
        handles=day_handles,
        title="Day",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.55),
        frameon=True,
        fontsize=8,
    )
    ax.add_artist(legend2)

    size_values = [10, 25, 50, 100]
    size_values = [n for n in size_values if n <= sub_vdf["N_Participants"].max()]
    size_handles = [
        plt.scatter(
            [],
            [],
            s=n * 6,
            color="gray",
            alpha=0.75,
            edgecolors="black",
            linewidth=0.6,
            label=f"N={n}",
        )
        for n in size_values
    ]
    ax.legend(
        handles=size_handles,
        title="Sample size",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.15),
        frameon=True,
        fontsize=8,
    )

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


p5 = os.path.join(OUTPUT_ROOT2, "05_replication_stddev_vs_unexplained_all.png")
plot_std_vs_unexplained(vdf, "Std Dev vs Unexplained Variance (color=strain, shape=day, size=N)", p5)
overview_plot_paths.append(("Std Dev vs Unexplained Variance — All Subtypes", p5))

for st in vdf["Subtype"].unique():
    sub = vdf[vdf["Subtype"] == st]
    p_st = os.path.join(OUTPUT_ROOT2, f"05_replication_stddev_vs_unexplained_{safe_name(st)}.png")
    plot_std_vs_unexplained(sub, f"Std Dev vs Unexplained Variance — Subtype {st}", p_st)
    overview_plot_paths.append((f"Std Dev vs Unexplained Variance — Subtype {st}", p_st))




# =================================================================
# 8. DESCRIPTIVE PLOTS per (subtype, virus, day), saved into a
#    subtype/strain/day folder tree
# =================================================================
BASE_DESC_DIR = os.path.join(OUTPUT_ROOT2, "replication_descriptive_by_subtype_strain_day")
os.makedirs(BASE_DESC_DIR, exist_ok=True)


def plot_group_matplotlib(sub_df, label, out_dir):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].hist(sub_df["log2_HAI"], bins=10, alpha=0.7, color="#4472C4", edgecolor="black")
    axes[0, 0].set_xlabel("Value (log2)")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("HAI distribution")

    axes[0, 1].scatter(sub_df["age"], sub_df["log2_HAI"], alpha=0.3, color="#ED7D31")
    axes[0, 1].set_xlabel("Age")
    axes[0, 1].set_ylabel("log2 HAI")
    axes[0, 1].set_title("Age vs HAI")

    if sub_df["sex"].nunique() > 1:
        sub_df.boxplot(column="log2_HAI", by="sex", ax=axes[1, 0])
        axes[1, 0].set_title("HAI by sex")
        axes[1, 0].set_xlabel("sex")
        axes[1, 0].set_ylabel("log2 HAI")
        fig.suptitle("")
    else:
        axes[1, 0].axis("off")

    axes[1, 1].axis("off")
    fig.suptitle(label, fontsize=11)
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, "replication_summary_plots.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    sub_df["log2_HAI"].describe().to_csv(os.path.join(out_dir, "replication_summary_stats.csv"))
    return fig_path


desc_plot_paths = {}  # key -> png path, for embedding in the HTML report

for (subtype_name, virus_name, day_val), sub_df in clean.groupby(["subtype", "virus", "day"]):
    vacc_status = sub_df["vaccinated"].iloc[0]
    n_obs = len(sub_df)
    label = f"{virus_name} | Subtype {subtype_name} | Day {day_val} | {vacc_status} | N={n_obs}"

    out_dir = os.path.join(
        BASE_DESC_DIR, safe_name(subtype_name), safe_name(virus_name),
        f"Day{int(day_val)}__{safe_name(vacc_status)}",
    )
    fig_path = plot_group_matplotlib(sub_df, label, out_dir)
    desc_plot_paths[(subtype_name, virus_name, day_val)] = fig_path

print(f"\nDescriptive plots saved under: {BASE_DESC_DIR}")





# ---------------------------------------------------------------
# PLOT: Scatter -- Residual Unexplained Variance % vs raw titer variability
# One PNG for All data, then one PNG per subtype,
# colored by strain, shaped by day,
# sized by N, with a sample-size legend
# ---------------------------------------------------------------
import os
import itertools
import numpy as np
import matplotlib.pyplot as plt

def safe_name(s):
    return (
        str(s)
        .replace("/", "_")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
    )

# Marker shapes to cycle through for "Day"
MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]

def short_label(virus, day):
    return f"{virus}, Day {day}"

# Add "All" first, then each subtype
plot_groups = ["All"] + sorted(vdf["Subtype"].dropna().unique())

for subtype in plot_groups:

    if subtype == "All":
        sub = vdf.copy()
        subtype_label = "All Subtypes"
        subtype_file = "All"
    else:
        sub = vdf[vdf["Subtype"] == subtype].copy()
        subtype_label = f"Subtype: {subtype}"
        subtype_file = subtype

    fig, ax = plt.subplots(figsize=(10, 6))

    # --- build strain -> color and day -> marker mappings ---
    strains = sorted(sub["Virus"].dropna().unique())
    days = sorted(sub["Day"].dropna().unique())

    strain_colors = plt.cm.tab20(np.linspace(0, 1, max(len(strains), 1)))
    strain_color_map = dict(zip(strains, strain_colors))

    marker_cycle = itertools.cycle(MARKER_CYCLE)
    day_marker_map = {day: next(marker_cycle) for day in days}

    # --- plot one small scatter per (strain, day) group ---
    for (virus, day), grp in sub.groupby(["Virus", "Day"]):
        ax.scatter(
            grp["Std_HAI_log2"],
            grp["Residual_Unexplained_pct"],
            s=grp["N_Participants"] * 4,
            alpha=0.75,
            color=strain_color_map[virus],
            marker=day_marker_map[day],
            edgecolors="black",
            linewidth=0.5,
        )

    # Point labels: strain + day, plus sample size
    for _, row in sub.iterrows():
        ax.annotate(
            f"{short_label(row['Virus'], row['Day'])}\nN={row['N_Participants']}",
            (row["Std_HAI_log2"], row["Residual_Unexplained_pct"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7.5,
        )

    ax.set_xlabel("HAI Titer Std Dev (log2 scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Unexplained variance (%)", fontsize=11, fontweight="bold")
    ax.set_title(
        f"Raw HAI Titer Variability vs. Residual Variance (%) After Covariate Adjustment\n{subtype_label}",
        fontsize=12,
        fontweight="bold",
    )
    ax.grid(True, alpha=0.3)

    # --- Legend 1: Strain (color) ---
    strain_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=strain_color_map[s], markeredgecolor="black",
                   markersize=8, label=s)
        for s in strains
    ]
    legend1 = ax.legend(
        handles=strain_handles,
        title="Strain",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        fontsize=8,
    )
    ax.add_artist(legend1)

    # --- Legend 2: Day (marker shape) ---
    day_handles = [
        plt.Line2D([0], [0], marker=day_marker_map[d], color="gray",
                   linestyle="None", markeredgecolor="black",
                   markersize=8, label=f"Day {d}")
        for d in days
    ]
    legend2 = ax.legend(
        handles=day_handles,
        title="Day",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.55),
        frameon=True,
        fontsize=8,
    )
    ax.add_artist(legend2)

    # --- Legend 3: Sample size (bubble size) ---
    size_values = [10, 25, 50, 100]
    size_values = [n for n in size_values if n <= sub["N_Participants"].max()]
    size_handles = [
        plt.scatter([], [], s=n * 4, color="gray", alpha=0.75,
                    edgecolors="black", linewidth=0.5, label=f"N={n}")
        for n in size_values
    ]
    ax.legend(
        handles=size_handles,
        title="Sample size",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.15),
        frameon=True,
        fontsize=8,
    )

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_ROOT2,
        f"03_replication_residual(relative)_vs_std_dev_{safe_name(subtype_file)}.png"
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

print("Saved all plots.")







# ---------------------------------------------------------------
# 3b. REGRESSION DIAGNOSTICS per (subtype, virus, day) -- QQ plot,
#     residuals vs fitted, scale-location, residuals vs leverage
# ---------------------------------------------------------------
from scipy import stats
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.stattools import durbin_watson
from statsmodels.stats.outliers_influence import OLSInfluence

BASE_DIAG_DIR = os.path.join(OUTPUT_ROOT2, "replication_diagnostics_by_subtype_strain_day")
os.makedirs(BASE_DIAG_DIR, exist_ok=True)

diag_rows = []


def plot_ols_diagnostics(fit, sub, label, out_dir):
    resid = np.asarray(fit.resid)
    fitted = np.asarray(fit.fittedvalues)
    std_resid = (resid - resid.mean()) / resid.std(ddof=1)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # 1. Residuals vs Fitted
    ax = axes[0, 0]
    ax.scatter(fitted, resid, alpha=0.5, edgecolor="black", linewidth=0.3)
    ax.axhline(0, color="red", linestyle="--", linewidth=1)
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess as _lowess
        sm_line = _lowess(resid, fitted, frac=0.6)
        ax.plot(sm_line[:, 0], sm_line[:, 1], color="blue", linewidth=1.5)
    except Exception:
        pass
    ax.set_xlabel("Fitted values")
    ax.set_ylabel("Residuals")
    ax.set_title("Residuals vs Fitted")
    ax.grid(alpha=0.3)

    # 2. Normal Q-Q
    ax = axes[0, 1]
    stats.probplot(std_resid, dist="norm", plot=ax)
    ax.set_title("Normal Q-Q")
    ax.grid(alpha=0.3)

    # 3. Scale-Location
    ax = axes[1, 0]
    sqrt_abs_std_resid = np.sqrt(np.abs(std_resid))
    ax.scatter(fitted, sqrt_abs_std_resid, alpha=0.5, edgecolor="black", linewidth=0.3)
    try:
        sm_line2 = _lowess(sqrt_abs_std_resid, fitted, frac=0.6)
        ax.plot(sm_line2[:, 0], sm_line2[:, 1], color="blue", linewidth=1.5)
    except Exception:
        pass
    ax.set_xlabel("Fitted values")
    ax.set_ylabel("sqrt(|Standardized Residuals|)")
    ax.set_title("Scale-Location")
    ax.grid(alpha=0.3)

    # 4. Residual Distribution (density)
    ax = axes[1, 1]
    ax.hist(resid, bins=15, density=True, alpha=0.7, color="#4472C4", edgecolor="black")
    xs = np.linspace(resid.min(), resid.max(), 200)
    ax.plot(xs, stats.norm.pdf(xs, resid.mean(), resid.std(ddof=1)), color="red", linewidth=2)
    ax.set_xlabel("Residuals")
    ax.set_ylabel("Density")
    ax.set_title("Residual Distribution")
    ax.grid(alpha=0.3)

    fig.suptitle(label, fontsize=11, fontweight="bold")
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, "replication_diagnostic_plots.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig_path

for (subtype_name, virus_name, day_val), fit in ols_results.items():
    sub = clean[
        (clean["subtype"] == subtype_name)
        & (clean["virus"] == virus_name)
        & (clean["day"] == day_val)
    ]
    n_obs = len(sub)
    vacc_status = vacc_label[(subtype_name, virus_name, day_val)]
    label = f"{virus_name} | Subtype {subtype_name} | Day {day_val} | {vacc_status} | N={n_obs}"

    out_dir = os.path.join(
        BASE_DIAG_DIR, safe_name(subtype_name), safe_name(virus_name),
        f"Day{int(day_val)}__{safe_name(vacc_status)}",
    )
    fig_path = plot_ols_diagnostics(fit, sub, label, out_dir)  # unchanged call

    # --- normality of residuals ---
    shapiro_stat, shapiro_p = stats.shapiro(fit.resid) if n_obs <= 5000 else (np.nan, np.nan)

    # --- heteroscedasticity (Breusch-Pagan) ---
    try:
        bp_stat, bp_p, _, _ = het_breuschpagan(fit.resid, fit.model.exog)
    except Exception:
        bp_stat, bp_p = np.nan, np.nan

    # --- autocorrelation of residuals ---
    dw_stat = durbin_watson(fit.resid)

    # --- influential points ---
    infl = OLSInfluence(fit)
    cooks_d = infl.cooks_distance[0]
    n_influential = int(np.sum(cooks_d > 4 / n_obs)) if n_obs > 0 else 0

    diag_rows.append({
        "Subtype": subtype_name,
        "Virus": virus_name,
        "Day": day_val,
        "Vaccinated": vacc_status,
        "N_Obs": n_obs,
        "Shapiro_W": round(shapiro_stat, 4) if pd.notna(shapiro_stat) else np.nan,
        "Shapiro_p": round(shapiro_p, 4) if pd.notna(shapiro_p) else np.nan,
        "Normal_Residuals": (shapiro_p > 0.05) if pd.notna(shapiro_p) else np.nan,
        "BreuschPagan_stat": round(bp_stat, 4) if pd.notna(bp_stat) else np.nan,
        "BreuschPagan_p": round(bp_p, 4) if pd.notna(bp_p) else np.nan,
        "Homoscedastic": (bp_p > 0.05) if pd.notna(bp_p) else np.nan,
        "DurbinWatson": round(dw_stat, 4),
        "N_Influential_CooksD": n_influential,
        "Diagnostic_Plot_Path": fig_path,
    })

    print(f"Saved diagnostics: {fig_path}")

diag_df = pd.DataFrame(diag_rows)
diag_df.to_csv(os.path.join(OUTPUT_ROOT2, "replication_regression_diagnostics_summary.csv"), index=False)
print("\nDiagnostics summary:")
print(diag_df.drop(columns=["Diagnostic_Plot_Path"]).to_string(index=False))










"""
================================================================================
SECTION 9. BUILD ONE SELF-CONTAINED HTML REPORT
================================================================================
Appends to the script above. Assumes the following already exist in memory
from the code you pasted:

    OUTPUT_ROOT2        - output directory
    strain_counts       - per-strain annotation counts (df)
    choice_df           - model choice log (df)
    summary_df          - regression summary (df)
    vdf                 - variance decomposition table (df)
    diag_df             - regression diagnostics summary (df)
    overview_plot_paths - list of (title, path) tuples, built during section 7
    desc_plot_paths     - dict {(subtype, virus, day): path}, built in section 8
    safe_name(s)        - helper already defined above

This section only ADDS to that state -- it does not redefine or recompute
anything already produced above. Images are embedded as base64 so the report
is a single file you can open/share without the folder of PNGs alongside it.

Run this after the rest of the script has completed.
"""

import os
import re
import base64
import html as html_lib
from datetime import datetime

REPORT_PATH2 = os.path.join(OUTPUT_ROOT2, "replication_hai_report.html")


# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
def safe_name2(s):
    """Local copy, independent of any global `safe_name` that may have been
    overwritten (e.g. accidentally reassigned to a string) earlier in the
    notebook session."""
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s


def img_to_b642(path):
    """Read a PNG from disk and return a data: URI, or None if missing."""
    if not path or not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        data = f.read()
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def esc2(x):
    return html_lib.escape(str(x))


def df_to_html_table2(df, table_id=None, max_rows=None):
    """Render a dataframe as a styled HTML table (no external deps)."""
    if df is None or len(df) == 0:
        return "<p><em>No data.</em></p>"
    d = df.copy()
    if max_rows is not None and len(d) > max_rows:
        d = d.head(max_rows)
        truncated_note = (
            f"<p class='note'>Showing first {max_rows} of {len(df)} rows.</p>"
        )
    else:
        truncated_note = ""

    cols = list(d.columns)
    id_attr = f" id='{esc2(table_id)}'" if table_id else ""
    parts = [f"<table class='data-table'{id_attr}>", "<thead><tr>"]
    for c in cols:
        parts.append(f"<th>{esc2(c)}</th>")
    parts.append("</tr></thead><tbody>")
    for _, row in d.iterrows():
        parts.append("<tr>")
        for c in cols:
            val = row[c]
            if isinstance(val, float):
                val = "" if val != val else round(val, 4)  # NaN check
            parts.append(f"<td>{esc2(val)}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table>")
    return truncated_note + "\n".join(parts)


def img_block2(title, path, anchor_id=None):
    b64 = img_to_b642(path)
    id_attr = f" id='{esc2(anchor_id)}'" if anchor_id else ""
    if b64 is None:
        return (
            f"<div class='plot-block'{id_attr}>"
            f"<h3>{esc2(title)}</h3>"
            f"<p class='note'>Image not found: {esc2(path)}</p></div>"
        )
    return (
        f"<div class='plot-block'{id_attr}>"
        f"<h3>{esc2(title)}</h3>"
        f"<img src='{b64}' alt='{esc2(title)}' loading='lazy'/>"
        f"</div>"
    )


def slug2(*parts):
    return safe_name2("_".join(str(p) for p in parts))


# ---------------------------------------------------------------
# Build TOC + sections
# ---------------------------------------------------------------
toc_items2 = []
body_sections2 = []

# --- Section: Overview plots (7a-7e) ---
toc_items2.append("<li><a href='#sec-overview'>Overview Plots</a><ul>")
overview_html2 = ["<h2 id='sec-overview'>Overview Plots</h2>"]
for i, (title, path) in enumerate(overview_plot_paths):
    anchor2 = f"ov-{i}-{slug2(title)}"
    toc_items2.append(f"<li><a href='#{anchor2}'>{esc2(title)}</a></li>")
    overview_html2.append(img_block2(title, path, anchor_id=anchor2))
toc_items2.append("</ul></li>")
body_sections2.append("\n".join(overview_html2))

# --- Section: Strain-level annotation counts ---
toc_items2.append("<li><a href='#sec-strain-counts'>Strain-Level Counts</a></li>")
body_sections2.append(
    "<h2 id='sec-strain-counts'>Strain-Level Annotation Counts</h2>"
    "<p class='note'>Computed once, across all days/visits, per strain.</p>"
    + df_to_html_table2(strain_counts)
)

# --- Section: Model choice log ---
toc_items2.append("<li><a href='#sec-model-choice'>Model Choice Log</a></li>")
body_sections2.append(
    "<h2 id='sec-model-choice'>Model Choice Log</h2>"
    + df_to_html_table2(choice_df, max_rows=500)
)

# --- Section: Variance decomposition table ---
toc_items2.append("<li><a href='#sec-variance'>Variance Decomposition Table</a></li>")
body_sections2.append(
    "<h2 id='sec-variance'>Variance Decomposition Table</h2>"
    + df_to_html_table2(vdf, max_rows=500)
)

# --- Section: Regression summary (coefficients) ---
toc_items2.append("<li><a href='#sec-summary'>Regression Summary (OLS Coefficients)</a></li>")
body_sections2.append(
    "<h2 id='sec-summary'>Regression Summary (OLS Coefficients)</h2>"
    + df_to_html_table2(summary_df, max_rows=1000)
)

# --- Section: Diagnostics summary ---
if "diag_df" in globals() and len(diag_df) > 0:
    toc_items2.append("<li><a href='#sec-diagnostics'>Regression Diagnostics Summary</a></li>")
    diag_table_df2 = diag_df.drop(columns=["Diagnostic_Plot_Path"], errors="ignore")
    body_sections2.append(
        "<h2 id='sec-diagnostics'>Regression Diagnostics Summary</h2>"
        + df_to_html_table2(diag_table_df2, max_rows=500)
    )

# --- Section: Per Subtype -> Strain -> Day (descriptive + diagnostic plots) ---
toc_items2.append("<li><a href='#sec-bystrain'>By Subtype &rarr; Strain &rarr; Day</a><ul>")
bystrain_html2 = ["<h2 id='sec-bystrain'>By Subtype &rarr; Strain &rarr; Day</h2>"]

# Build lookup from diag_df for the diagnostic plot paths
diag_lookup2 = {}
if "diag_df" in globals():
    for _, r in diag_df.iterrows():
        diag_lookup2[(r["Subtype"], r["Virus"], r["Day"])] = r.get("Diagnostic_Plot_Path")

# organize keys by subtype for nested TOC
from collections import defaultdict
by_subtype2 = defaultdict(list)
for key in desc_plot_paths.keys():
    subtype_name2, virus_name2, day_val2 = key
    by_subtype2[subtype_name2].append(key)

for subtype_name2 in sorted(by_subtype2.keys(), key=str):
    subtype_anchor2 = f"subtype-{slug2(subtype_name2)}"
    toc_items2.append(
        f"<li><a href='#{subtype_anchor2}'>{esc2(subtype_name2)}</a><ul>"
    )
    bystrain_html2.append(f"<h3 id='{subtype_anchor2}'>Subtype: {esc2(subtype_name2)}</h3>")

    # strain annotation block (once per strain within this subtype)
    strains_in_subtype2 = sorted(
        {k[1] for k in by_subtype2[subtype_name2]}, key=str
    )
    for virus_name2 in strains_in_subtype2:
        ann_row2 = strain_counts[
            (strain_counts["subtype"] == subtype_name2)
            & (strain_counts["virus"] == virus_name2)
        ]
        if len(ann_row2):
            r = ann_row2.iloc[0]
            ann_text2 = (
                f"N_Samples={r['N_Samples']} | N_Studies={r['N_Studies']} | "
                f"N_Cohorts={r['N_Cohorts']} | N_Days={r['N_Days']}"
            )
        else:
            ann_text2 = "N_Samples=? | N_Studies=? | N_Cohorts=? | N_Days=?"

        strain_anchor2 = f"strain-{slug2(subtype_name2, virus_name2)}"
        toc_items2.append(f"<li><a href='#{strain_anchor2}'>{esc2(virus_name2)}</a></li>")
        bystrain_html2.append(
            f"<div class='strain-block' id='{strain_anchor2}'>"
            f"<h4>{esc2(virus_name2)}</h4>"
            f"<p class='annotation'>{esc2(ann_text2)}</p>"
        )

        days_for_strain2 = sorted(
            d for (s, v, d) in by_subtype2[subtype_name2] if v == virus_name2
        )
        for day_val2 in days_for_strain2:
            key = (subtype_name2, virus_name2, day_val2)
            desc_path2 = desc_plot_paths.get(key)
            diag_path2 = diag_lookup2.get(key)

            bystrain_html2.append(f"<div class='day-block'><h5>Day {esc2(day_val2)}</h5>")
            bystrain_html2.append(
                img_block2(f"Descriptive plots — Day {day_val2}", desc_path2)
            )
            if diag_path2:
                bystrain_html2.append(
                    img_block2(f"Regression diagnostics — Day {day_val2}", diag_path2)
                )
            bystrain_html2.append("</div>")  # day-block

        bystrain_html2.append("</div>")  # strain-block

    toc_items2.append("</ul></li>")

toc_items2.append("</ul></li>")
body_sections2.append("\n".join(bystrain_html2))


# ---------------------------------------------------------------
# Assemble final HTML
# ---------------------------------------------------------------
generated_at2 = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

CSS2 = """
:root {
  --bg: #f7f8fa;
  --panel: #ffffff;
  --border: #e2e5ea;
  --text: #1f2430;
  --muted: #667085;
  --accent: #4472C4;
  --accent2: #ED7D31;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  display: flex;
}
#toc {
  width: 300px;
  min-width: 300px;
  height: 100vh;
  overflow-y: auto;
  position: sticky;
  top: 0;
  background: var(--panel);
  border-right: 1px solid var(--border);
  padding: 20px 16px;
  font-size: 13.5px;
}
#toc h1 {
  font-size: 15px;
  margin: 0 0 4px 0;
}
#toc .meta {
  color: var(--muted);
  font-size: 11.5px;
  margin-bottom: 16px;
}
#toc ul {
  list-style: none;
  padding-left: 14px;
  margin: 4px 0;
}
#toc > ul { padding-left: 0; }
#toc li { margin: 3px 0; }
#toc a {
  color: var(--text);
  text-decoration: none;
}
#toc a:hover { color: var(--accent); text-decoration: underline; }
#content {
  flex: 1;
  padding: 28px 40px 80px 40px;
  max-width: 1200px;
}
h2 {
  border-bottom: 2px solid var(--accent);
  padding-bottom: 6px;
  margin-top: 48px;
}
h3 {
  color: var(--accent);
  margin-top: 36px;
}
h4 {
  margin-top: 24px;
  background: #eef1f8;
  padding: 6px 10px;
  border-left: 3px solid var(--accent);
}
h5 {
  margin-top: 16px;
  color: var(--muted);
}
.annotation {
  font-size: 13px;
  color: var(--muted);
  font-family: monospace;
  background: #f0f1f4;
  display: inline-block;
  padding: 4px 8px;
  border-radius: 4px;
}
.strain-block {
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 12px 18px 18px 18px;
  margin: 16px 0;
  background: var(--panel);
}
.day-block {
  border-top: 1px dashed var(--border);
  padding-top: 10px;
  margin-top: 14px;
}
.plot-block {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 14px;
  margin: 16px 0;
}
.plot-block img {
  max-width: 100%;
  height: auto;
  display: block;
  margin-top: 8px;
  border-radius: 4px;
}
.note {
  color: var(--muted);
  font-size: 12.5px;
  font-style: italic;
}
table.data-table {
  border-collapse: collapse;
  width: 100%;
  font-size: 12.5px;
  margin: 10px 0 24px 0;
  background: var(--panel);
}
table.data-table th, table.data-table td {
  border: 1px solid var(--border);
  padding: 5px 9px;
  text-align: left;
  white-space: nowrap;
}
table.data-table th {
  background: #eef1f8;
  position: sticky;
  top: 0;
}
table.data-table tr:nth-child(even) { background: #fafbfc; }
.table-wrap {
  overflow-x: auto;
}
"""

html_doc2 = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HAI Replication Regression Report</title>
<style>{CSS2}</style>
</head>
<body>
<nav id="toc">
  <h1>HAI Replication Report</h1>
  <div class="meta">Generated {esc2(generated_at2)}</div>
  <ul>
    {"".join(toc_items2)}
  </ul>
</nav>
<main id="content">
  <h1>HAI Regression &amp; Variance Analysis</h1>
  <p class="note">Organized by Subtype &rarr; Strain &rarr; Day. Generated {esc2(generated_at2)}.</p>
  {"".join(f"<div class='table-wrap'>{s}</div>" if "<table" in s else s for s in body_sections2)}
</main>
</body>
</html>
"""

with open(REPORT_PATH2, "w", encoding="utf-8") as f:
    f.write(html_doc2)

print(f"\nHTML report written to: {REPORT_PATH2}")

