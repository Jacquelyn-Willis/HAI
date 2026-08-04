
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


final_merge["Race_demo"].value_counts(dropna=False) 

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


final_merge["Cohort for regression"] = final_merge["Description_merge2"].fillna(final_merge["Cohort"])

### pre HAI demographics and distribution plots


# ---- config ----
STRAIN_COL = "Virus"
COHORT_COL = "Cohort for regression"
DAY_COL = "Study Time Collected"
DAY_UNIT_COL = "Study Time Collected Unit"
SUBTYPE_COL = "subtype"
MIN_N = 10

outdir = os.path.join(SCRATCH, "final_hai_demo_plots_html")
os.makedirs(outdir, exist_ok=True)


def make_group_figure(df, title):
    """Build the 2x2 subplot figure for a single strain/day slice."""
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("HAI distribution", "Age vs HAI", "HAI by sex", ""),
        vertical_spacing=0.12,
        horizontal_spacing=0.10,
    )

    # Histogram
    fig.add_trace(
        go.Histogram(
            x=df["log2_HAI"], nbinsx=10, name="log2_HAI", opacity=0.7
        ),
        row=1,
        col=1,
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
        row=1,
        col=2,
    )

    # Box: sex vs HAI
    if df["Gender_demo"].nunique() > 1:
        for gender, gdf in df.groupby("Gender_demo"):
            fig.add_trace(
                go.Box(y=gdf["log2_HAI"], name=str(gender), boxmean=True),
                row=2,
                col=1,
            )

    fig.update_xaxes(title_text="Value (log2)", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_xaxes(title_text="Age", row=1, col=2)
    fig.update_yaxes(title_text="log2 HAI", row=1, col=2)
    fig.update_xaxes(title_text="Gender", row=2, col=1)
    fig.update_yaxes(title_text="log2 HAI", row=2, col=1)

    fig.update_layout(
        title=title,
        width=1100,
        height=750,
        showlegend=False,
        template="plotly_white",
    )
    return fig


def slugify(s):
    return (
        str(s)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("—", "_")
        .replace(":", "")
    )


# ---- build one combined HTML report ----
strains = sorted(final_merge[STRAIN_COL].dropna().unique())

toc_entries = []  # (anchor_id, display_label)
body_sections = []  # rendered HTML chunks
total_plots = 0

for strain in strains:
    strain_df = final_merge[final_merge[STRAIN_COL] == strain]
    days = sorted(strain_df[DAY_COL].dropna().unique())

    strain_anchor = f"strain_{slugify(strain)}"
    strain_section_plots = []

    for day in days:
        day_df = strain_df[strain_df[DAY_COL] == day]
        n = len(day_df)

        if n < MIN_N:
            print(f"Skipping {strain} / Day {day} (n={n} < {MIN_N})")
            continue

        subtype = (
            day_df[SUBTYPE_COL].dropna().iloc[0]
            if SUBTYPE_COL in day_df.columns
            and day_df[SUBTYPE_COL].notna().any()
            else "NA"
        )
        day_unit = (
            day_df[DAY_UNIT_COL].dropna().iloc[0]
            if DAY_UNIT_COL in day_df.columns
            and day_df[DAY_UNIT_COL].notna().any()
            else ""
        )
        cohorts = ", ".join(
            sorted(day_df[COHORT_COL].dropna().unique().astype(str))
        )

        label = f"{strain} — Day {day} {day_unit} — Subtype: {subtype} — N={n}"
        print(f"\n=== {label} ===")
        print(day_df["log2_HAI"].describe())

        fig = make_group_figure(day_df, label)

        day_anchor = f"{strain_anchor}_day_{slugify(day)}"
        toc_entries.append((day_anchor, label))

        # include_plotlyjs='cdn' ensures a clean single-script dependency without inline bundle bloating
        section_html = f"""
        <div id="{day_anchor}" class="plot-section">
            <h3>{label}</h3>
            <p class="meta">Cohort(s): {cohorts}</p>
            {fig.to_html(full_html=False, include_plotlyjs="cdn")}
        </div>
        """
        strain_section_plots.append(section_html)
        total_plots += 1

    if strain_section_plots:
        body_sections.append(f'<h2 id="{strain_anchor}">{strain}</h2>')
        body_sections.extend(strain_section_plots)

# ---- assemble final HTML ----
toc_html = (
    "<ul>"
    + "".join(
        f'<li><a href="#{anchor}">{label}</a></li>'
        for anchor, label in toc_entries
    )
    + "</ul>"
)

html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>HAI Distribution Report by Strain and Day</title>
<style>
    body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.5; }}
    h1 {{ border-bottom: 2px solid #333; padding-bottom: 8px; }}
    h2 {{ margin-top: 50px; color: #1a5276; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
    h3 {{ margin-top: 30px; color: #444; }}
    .meta {{ color: #666; font-size: 0.9em; }}
    .plot-section {{ margin-bottom: 40px; }}
    #toc {{ background: #f7f7f7; padding: 15px 25px; border-radius: 8px; margin-bottom: 30px; }}
    #toc a {{ text-decoration: none; color: #1a5276; }}
    #toc a:hover {{ text-decoration: underline; }}
    #toc ul {{ margin-bottom: 0; }}
</style>
</head>
<body>
<h1>HAI Distribution by Strain / Day</h1>
<p>Total plots: {total_plots} (groups with N &lt; {MIN_N} excluded)</p>
<div id="toc"><strong>Jump to:</strong>{toc_html}</div>
{''.join(body_sections)}
</body>
</html>
"""

filepath = os.path.join(outdir, "hai_report_by_strain_day.html")
with open(filepath, "w", encoding="utf-8") as f:
    f.write(html)

print(f"\nSaved combined report: {filepath}")


















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




OUTPUT_ROOT = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/final_hai_regression_results"
os.makedirs(OUTPUT_ROOT, exist_ok=True)

# ---------------------------------------------------------------
# 1. LOAD + CLEAN (same as original)
# ---------------------------------------------------------------
df = final_merge.copy()
df.columns = df.columns.str.strip()

COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_demo",
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










# =================================================================
# 9. BUILD SINGLE-PAGE HTML REPORT WITH STICKY TOC AND INDEX
#    Subtype-aware + includes diagnostics plots
# =================================================================

def img_to_base64(path):
    """Converts an image file to a base64 string for direct HTML embedding."""
    if not path or not os.path.exists(path):
        return ""
    with open(path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded_string}"


print("\nGenerating unified HTML report...")

# -----------------------------------------------------------------
# Build lookup tables for report content
# -----------------------------------------------------------------
# One descriptive plot per (Subtype, Virus, Day)
desc_lookup = {}
for key, path in desc_plot_paths.items():
    desc_lookup[key] = path

# One diagnostic plot per (Subtype, Virus, Day) and one chosen model row
diag_lookup = {}
perf_lookup = {}

for _, row in perf_df.iterrows():
    key = (row["Subtype"], row["Virus"], row["Day"])
    perf_lookup[key] = row

    # keep first valid diagnostic plot path for each group
    if key not in diag_lookup and pd.notna(row.get("Diagnostic_Plot", np.nan)):
        diag_lookup[key] = row["Diagnostic_Plot"]

# 1. Gather all global overview plot base64 strings
overview_b64 = []
for title, path in overview_plot_paths:
    if os.path.exists(path):
        overview_b64.append((title, img_to_base64(path)))

# 2. Build TOC structure grouped by Subtype -> Strain -> Day
tree = {}
for (subtype_name, virus_name, day_val), sub_df in clean.groupby(["subtype", "virus", "day"]):
    st = str(subtype_name)
    vi = str(virus_name)
    d = int(day_val)

    if st not in tree:
        tree[st] = {}
    if vi not in tree[st]:
        tree[st][vi] = {}

    key = (subtype_name, virus_name, day_val)
    perf_row = perf_lookup.get(key, None)

    tree[st][vi][d] = {
        "n_obs": len(sub_df),
        "vacc": sub_df["vaccinated"].iloc[0],
        "desc_img": img_to_base64(desc_lookup.get(key, "")),
        "diag_img": img_to_base64(diag_lookup.get(key, "")),
        "model_used": perf_row["Model"] if perf_row is not None else "N/A",
        "rmse": perf_row["RMSE"] if perf_row is not None else np.nan,
        "aic": perf_row["AIC"] if perf_row is not None else np.nan,
        "bic": perf_row["BIC"] if perf_row is not None else np.nan,
        "r2_marg": perf_row["Marginal_R2"] if perf_row is not None else np.nan,
        "r2_cond": perf_row["Conditional_R2"] if perf_row is not None else np.nan,
        "shapiro_p": perf_row["Shapiro_p"] if perf_row is not None else np.nan,
        "normal_ok": perf_row["Residuals_Normal_at_0.05"] if perf_row is not None else np.nan,
    }

# 3. Construct HTML Content
html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HAI Regression & Variance Analysis Report</title>
    <style>
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --bg-color: #f8f9fa;
            --sidebar-width: 320px;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            margin: 0;
            padding: 0;
            background-color: var(--bg-color);
            color: #333;
            display: flex;
        }}
        #sidebar {{
            width: var(--sidebar-width);
            height: 100vh;
            position: sticky;
            top: 0;
            background: #ffffff;
            border-right: 1px solid #e0e0e0;
            box-sizing: border-box;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            z-index: 100;
        }}
        #sidebar-header {{
            padding: 15px;
            background: var(--primary-color);
            color: white;
        }}
        #sidebar-header h2 {{
            margin: 0 0 10px 0;
            font-size: 1.1rem;
        }}
        #search-input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ccc;
            border-radius: 4px;
            box-sizing: border-box;
            font-size: 0.9rem;
        }}
        #toc-container {{
            flex: 1;
            overflow-y: auto;
            padding: 15px;
        }}
        .toc-subtype {{
            font-weight: bold;
            margin-top: 10px;
            font-size: 0.95rem;
            color: var(--primary-color);
        }}
        .toc-strain {{
            margin-left: 10px;
            font-size: 0.85rem;
            margin-top: 5px;
            color: #555;
        }}
        .toc-day-list {{
            margin-left: 15px;
            list-style-type: none;
            padding-left: 0;
            font-size: 0.8rem;
        }}
        .toc-day-list li {{
            margin: 3px 0;
        }}
        .toc-day-list a {{
            color: var(--secondary-color);
            text-decoration: none;
        }}
        .toc-day-list a:hover {{
            text-decoration: underline;
        }}
        #main-content {{
            flex: 1;
            padding: 30px;
            max-width: 1200px;
            box-sizing: border-box;
        }}
        h1 {{
            color: var(--primary-color);
            border-bottom: 2px solid var(--primary-color);
            padding-bottom: 10px;
        }}
        h2 {{
            color: var(--primary-color);
            margin-top: 40px;
            border-bottom: 1px solid #ccc;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #444;
            margin-top: 25px;
        }}
        .card {{
            background: white;
            border-radius: 6px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.08);
            padding: 20px;
            margin-bottom: 25px;
        }}
        .annotation-box {{
            background: #eef7fc;
            border-left: 4px solid var(--secondary-color);
            padding: 10px 15px;
            margin: 10px 0 20px 0;
            font-weight: 500;
            font-size: 0.9rem;
        }}
        .plot-img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #eee;
            border-radius: 4px;
            margin-top: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 0.85rem;
        }}
        th, td {{
            text-align: left;
            padding: 8px 12px;
            border: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .badge {{
            display: inline-block;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 0.75rem;
            font-weight: bold;
            color: white;
        }}
        .badge-lmm {{ background-color: #27ae60; }}
        .badge-ols {{ background-color: #e67e22; }}
        .badge-good {{ background-color: #2ecc71; }}
        .badge-bad {{ background-color: #e74c3c; }}
    </style>
</head>
<body>

    <div id="sidebar">
        <div id="sidebar-header">
            <h2>HAI Analysis Index</h2>
            <input type="text" id="search-input" onkeyup="filterTOC()" placeholder="Search Strain or Subtype...">
        </div>
        <div id="toc-container">
            <div class="toc-subtype"><a href="#overview" style="color: inherit; text-decoration: none;">📊 Global Summaries</a></div>
            <hr style="border: 0; border-top: 1px solid #eee; margin: 10px 0;">
"""

# Sidebar TOC
for subtype, strains in tree.items():
    sub_id = safe_name(subtype)
    html_content += f'<div class="toc-group"><div class="toc-subtype">Subtype: {subtype}</div>\n'
    for strain, days in strains.items():
        strain_id = safe_name(strain)
        html_content += f'<div class="toc-strain">🦠 {strain[:25]}...</div><ul class="toc-day-list">\n'
        for day in sorted(days.keys()):
            target_id = f"{sub_id}__{strain_id}__D{day}"
            html_content += f'<li><a href="#{target_id}">Day {day} (N={days[day]["n_obs"]})</a></li>\n'
        html_content += '</ul>\n'
    html_content += '</div>\n'

html_content += """
        </div>
    </div>

    <div id="main-content">
        <h1>HAI Regression & Variance Analysis</h1>

        <section id="overview" class="card">
            <h2>Global Subtype & Variance Summaries</h2>
            <p>Aggregated diagnostics and model fits across all evaluated viral subtypes and study cohorts.</p>
"""

for title, b64_str in overview_b64:
    html_content += f"""
            <h3>{title}</h3>
            <img class="plot-img" src="{b64_str}" alt="{title}">
"""

html_content += """
        </section>

        <h2>Subtype & Strain Diagnostics</h2>
"""

# Main content
for subtype, strains in tree.items():
    sub_id = safe_name(subtype)
    html_content += f'<div id="sec_{sub_id}">'

    html_content += f"""
        <div class="card">
            <h2>Subtype: {subtype}</h2>
        </div>
    """

    for strain, days in strains.items():
        strain_id = safe_name(strain)
        annot = strain_annotation_text(strain)

        html_content += f"""
        <div class="card">
            <h3>Strain: {strain} <small>({subtype})</small></h3>
            <div class="annotation-box">
                📌 <strong>Strain Annotation Summary:</strong> {annot}
            </div>
        """

        for day in sorted(days.keys()):
            target_id = f"{sub_id}__{strain_id}__D{day}"
            data = days[day]

            model_cls = "badge-lmm" if data["model_used"] == "LMM" else "badge-ols"
            normal_cls = "badge-good" if bool(data["normal_ok"]) else "badge-bad"

            html_content += f"""
            <div id="{target_id}" style="border-top: 1px dashed #ccc; padding-top: 15px; margin-top: 20px;">
                <h4>Day {day} Analysis</h4>
                <table>
                    <tr>
                        <th>Observation Count (N)</th>
                        <th>Vaccination Status</th>
                        <th>Model Selected</th>
                        <th>RMSE</th>
                        <th>AIC</th>
                        <th>BIC</th>
                        <th>Marginal R²</th>
                        <th>Conditional R²</th>
                        <th>Shapiro p</th>
                    </tr>
                    <tr>
                        <td>{data['n_obs']}</td>
                        <td>{data['vacc']}</td>
                        <td><span class="badge {model_cls}">{data['model_used']}</span></td>
                        <td>{f"{data['rmse']:.3f}" if pd.notna(data['rmse']) else "N/A"}</td>
                        <td>{f"{data['aic']:.2f}" if pd.notna(data['aic']) else "N/A"}</td>
                        <td>{f"{data['bic']:.2f}" if pd.notna(data['bic']) else "N/A"}</td>
                        <td>{f"{data['r2_marg']:.3f}" if pd.notna(data['r2_marg']) else "N/A"}</td>
                        <td>{f"{data['r2_cond']:.3f}" if pd.notna(data['r2_cond']) else "N/A"}</td>
                        <td><span class="badge {normal_cls}">{f"{data['shapiro_p']:.4f}" if pd.notna(data['shapiro_p']) else "N/A"}</span></td>
                    </tr>
                </table>
            """

            if data["desc_img"]:
                html_content += f'<h4>Descriptive Plots</h4><img class="plot-img" src="{data["desc_img"]}" alt="Descriptive Plots D{day}">'

            if data["diag_img"]:
                html_content += f'<h4>Model Diagnostics</h4><img class="plot-img" src="{data["diag_img"]}" alt="Diagnostics D{day}">'

            html_content += "</div>"

        html_content += "</div>"

    html_content += "</div>"

html_content += """
    </div>

    <script>
        function filterTOC() {
            var input = document.getElementById('search-input');
            var filter = input.value.toLowerCase();
            var groups = document.getElementsByClassName('toc-group');

            for (var i = 0; i < groups.length; i++) {
                var text = groups[i].innerText.toLowerCase();
                if (text.includes(filter)) {
                    groups[i].style.display = "";
                } else {
                    groups[i].style.display = "none";
                }
            }
        }
    </script>
</body>
</html>
"""

# 4. Save Final Report
html_file_path = os.path.join(OUTPUT_ROOT, "hai_variance_report.html")
with open(html_file_path, "w", encoding="utf-8") as f:
    f.write(html_content)

print(f"✅ HTML Report successfully built and saved to:\n   {html_file_path}")





####plots needed for lab notebook

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
AGE_COL = "Age Reported_demo"
GENDER_COL = "Gender_demo"
HAI_COL = "log2_HAI"
MIN_N = 10  # per strain/day slice; days below this are dropped before plotting


outdir = os.path.join(SCRATCH, "final_hai_strain_facet_png")
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
















##post confirmation regression analysis 

#1 filter final_merge for final strain list and only d0 and d28 for a single cohort 
strain_list = ['A/Solomon Islands/3/2006', 'A/California/7/2009', 'A/Perth/16/2009',
               'A/Victoria/361/2011', 'B/Malaysia/2506/2004', 'B/Brisbane/3/2007', 'B/Wisconsin/01/2010']

post_df = final_merge[(final_merge['Virus'].isin(strain_list)) & 
                      (final_merge['Study Time Collected'].isin([0, 28]))]


#h1n1
h1n1 = post_df[post_df['Virus'].isin(['A/California/7/2009', 'A/Solomon Islands/3/2006'])]
solomon = h1n1[h1n1['Cohort for regression'] == 'Older participants aged 60 to 89 years, vaccinated with Fluzone']
cali = h1n1[h1n1['Cohort for regression'] == '150 healthy adults, 50-74 yo']


#h3n2
h3n2 = post_df[(post_df['Virus'].isin(['A/Perth/16/2009', 'A/Victoria/361/2011']))]
perth = h3n2[(h3n2['Virus'] == 'A/Perth/16/2009') & (h3n2['Cohort for regression'] == '150 healthy adults, 50-74 yo')]
vic= h3n2[(h3n2['Virus'] == 'A/Victoria/361/2011') & (h3n2['Cohort for regression'] == 'Healthy Adults 2012 - 2013')]


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




OUTPUT_ROOT = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/replication_hai_regression_results"
os.makedirs(OUTPUT_ROOT, exist_ok=True)

# ---------------------------------------------------------------
# 1. LOAD + CLEAN (same as original)
# ---------------------------------------------------------------
df = rep_df.copy()
df.columns = df.columns.str.strip()

COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_demo",
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
strain_counts.to_csv(os.path.join(OUTPUT_ROOT, "replication_strain_level_counts.csv"), index=False)
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
summary_df.to_csv(os.path.join(OUTPUT_ROOT, "replication_hai_regression_summary.csv"), index=False)

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
choice_df.to_csv(os.path.join(OUTPUT_ROOT, "replication_hai_model_choice_log.csv"), index=False)

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

vdf.to_csv(os.path.join(OUTPUT_ROOT, "replication_variance_decomposition_table.csv"), index=False)

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
        OUTPUT_ROOT,
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
        OUTPUT_ROOT,
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
        OUTPUT_ROOT,
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
        OUTPUT_ROOT,
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


p5 = os.path.join(OUTPUT_ROOT, "05_replication_stddev_vs_unexplained_all.png")
plot_std_vs_unexplained(vdf, "Std Dev vs Unexplained Variance (color=strain, shape=day, size=N)", p5)
overview_plot_paths.append(("Std Dev vs Unexplained Variance — All Subtypes", p5))

for st in vdf["Subtype"].unique():
    sub = vdf[vdf["Subtype"] == st]
    p_st = os.path.join(OUTPUT_ROOT, f"05_replication_stddev_vs_unexplained_{safe_name(st)}.png")
    plot_std_vs_unexplained(sub, f"Std Dev vs Unexplained Variance — Subtype {st}", p_st)
    overview_plot_paths.append((f"Std Dev vs Unexplained Variance — Subtype {st}", p_st))




# =================================================================
# 8. DESCRIPTIVE PLOTS per (subtype, virus, day), saved into a
#    subtype/strain/day folder tree
# =================================================================
BASE_DESC_DIR = os.path.join(OUTPUT_ROOT, "replication_descriptive_by_subtype_strain_day")
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
        OUTPUT_ROOT,
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

BASE_DIAG_DIR = os.path.join(OUTPUT_ROOT, "replication_diagnostics_by_subtype_strain_day")
os.makedirs(BASE_DIAG_DIR, exist_ok=True)

diag_rows = []

def plot_ols_diagnostics(fit, sub, label, out_dir):
    resid = fit.resid
    fitted = fit.fittedvalues
    influence = OLSInfluence(fit)
    std_resid = influence.resid_studentized_internal
    leverage = influence.hat_matrix_diag
    cooks_d = influence.cooks_distance[0]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Residuals vs Fitted
    axes[0, 0].scatter(fitted, resid, alpha=0.5, color="#4472C4", edgecolor="black", linewidth=0.3)
    axes[0, 0].axhline(0, color="red", ls="--", lw=1)
    axes[0, 0].set_xlabel("Fitted values")
    axes[0, 0].set_ylabel("Residuals")
    axes[0, 0].set_title("Residuals vs Fitted")
    # lowess trend line if enough points
    if len(fitted) >= 10:
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            sm_fit = lowess(resid, fitted, frac=0.6)
            axes[0, 0].plot(sm_fit[:, 0], sm_fit[:, 1], color="orange", lw=1.5)
        except Exception:
            pass

    # 2. QQ plot of standardized residuals
    (osm, osr), (slope, intercept, r) = stats.probplot(std_resid, dist="norm")
    axes[0, 1].scatter(osm, osr, alpha=0.6, color="#ED7D31", edgecolor="black", linewidth=0.3)
    axes[0, 1].plot(osm, slope * osm + intercept, color="red", lw=1.5)
    axes[0, 1].set_xlabel("Theoretical quantiles")
    axes[0, 1].set_ylabel("Standardized residuals")
    axes[0, 1].set_title("Normal Q-Q")

    fig.suptitle(label, fontsize=11)
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
    fig_path = plot_ols_diagnostics(fit, sub, label, out_dir)

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
diag_df.to_csv(os.path.join(OUTPUT_ROOT, "replication_regression_diagnostics_summary.csv"), index=False)
print("\nDiagnostics summary:")
print(diag_df.drop(columns=["Diagnostic_Plot_Path"]).to_string(index=False))

# %%
# ---------------------------------------------------------------
# 9. BUILD SELF-CONTAINED HTML REPORT
#    Organized: Overview -> Subtype -> Strain -> Day
#    All images embedded as base64 so the report is a single
#    portable file with no dependency on relative image paths.
# ---------------------------------------------------------------
import base64
from datetime import datetime

def img_to_b64(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")

def img_tag(path, alt="", max_width="100%"):
    b64 = img_to_b64(path)
    if b64 is None:
        return f'<p class="missing">[missing: {alt}]</p>'
    return f'<img src="data:image/png;base64,{b64}" alt="{alt}" style="max-width:{max_width};height:auto;border:1px solid #ddd;border-radius:4px;">'

def df_to_html_table(d, float_cols=None):
    if d is None or len(d) == 0:
        return "<p class='missing'>No data.</p>"
    d = d.copy()
    return d.to_html(index=False, classes="stat-table", border=0, na_rep="—")

# --- QQ plot path lookup (per-group standalone, from section 7f-ii) ---
qq_plot_paths = {}
for (subtype_name, virus_name, day_val) in ols_results.keys():
    vacc_status = vacc_label[(subtype_name, virus_name, day_val)]
    p = os.path.join(
        BASE_QQ_DIR, safe_name(subtype_name), safe_name(virus_name),
        f"Day{int(day_val)}__{safe_name(vacc_status)}", "replication_qq_plot.png",
    )
    qq_plot_paths[(subtype_name, virus_name, day_val)] = p

# --- diagnostic plot path lookup ---
diag_plot_paths = dict(zip(
    zip(diag_df["Subtype"], diag_df["Virus"], diag_df["Day"]),
    diag_df["Diagnostic_Plot_Path"],
))
diag_lookup = diag_df.set_index(["Subtype", "Virus", "Day"]).to_dict("index")

CSS = """
<style>
  :root {
    --blue: #4472C4; --orange: #ED7D31; --green: #70AD47;
    --bg: #f7f8fa; --card: #ffffff; --border: #e2e5ea; --text: #1f2937;
  }
  * { box-sizing: border-box; }
  body {
    font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg); color: var(--text); margin: 0; padding: 0;
    line-height: 1.5;
  }
  header {
    background: linear-gradient(135deg, #2c3e6b, #4472C4);
    color: white; padding: 32px 40px;
  }
  header h1 { margin: 0 0 6px 0; font-size: 26px; }
  header p { margin: 0; opacity: 0.85; font-size: 14px; }
  nav.toc {
    background: var(--card); margin: 20px 40px; padding: 20px 28px;
    border: 1px solid var(--border); border-radius: 8px;
  }
  nav.toc h2 { margin-top: 0; font-size: 16px; text-transform: uppercase; letter-spacing: .05em; color: #555; }
  nav.toc ul { columns: 3; column-gap: 24px; padding-left: 18px; }
  nav.toc li { break-inside: avoid; margin-bottom: 4px; font-size: 13.5px; }
  nav.toc a { color: var(--blue); text-decoration: none; }
  nav.toc a:hover { text-decoration: underline; }
  section.block { margin: 0 40px 36px 40px; }
  h2.subtype-header {
    background: #2c3e6b; color: white; padding: 12px 20px; border-radius: 6px;
    font-size: 19px; margin-bottom: 18px; position: sticky; top: 0; z-index: 5;
  }
  h3.strain-header {
    background: #dde5f5; padding: 10px 18px; border-left: 5px solid var(--blue);
    border-radius: 4px; font-size: 16px; margin: 22px 0 12px 0;
  }
  .day-card {
    background: var(--card); border: 1px solid var(--border); border-radius: 8px;
    padding: 20px; margin-bottom: 20px;
  }
  .day-card h4 { margin: 0 0 12px 0; font-size: 15px; color: #2c3e6b; }
  .badge {
    display: inline-block; font-size: 11px; padding: 2px 8px; border-radius: 10px;
    margin-left: 8px; font-weight: 600;
  }
  .badge.vacc { background: #d9f2e3; color: #1e7a46; }
  .badge.notvacc { background: #f2e3d9; color: #a15c1e; }
  .badge.skipped { background: #f0f0f0; color: #888; }
  .badge.flag { background: #fde2e2; color: #b91c1c; }
  .badge.ok { background: #e2f5e9; color: #157347; }
  .img-grid {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
    gap: 14px; margin-top: 10px;
  }
  .img-grid figure { margin: 0; }
  .img-grid figcaption { font-size: 12px; color: #666; margin-top: 4px; text-align: center; }
  .stat-table { border-collapse: collapse; width: 100%; font-size: 12.5px; margin-top: 8px; }
  .stat-table th, .stat-table td { padding: 5px 10px; border-bottom: 1px solid #eee; text-align: right; }
  .stat-table th { background: #f2f4f8; text-align: right; font-weight: 600; }
  .stat-table td:first-child, .stat-table th:first-child { text-align: left; }
  .diag-summary { display: flex; flex-wrap: wrap; gap: 10px; margin: 10px 0; }
  .diag-pill { background: #f2f4f8; border: 1px solid var(--border); border-radius: 6px; padding: 6px 12px; font-size: 12.5px; }
  .missing { color: #aaa; font-style: italic; font-size: 13px; }
  .overview-grid {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
    gap: 18px;
  }
  .overview-grid figure { background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 12px; margin: 0; }
  .overview-grid figcaption { font-size: 13px; font-weight: 600; margin-bottom: 8px; color: #2c3e6b; }
  .back-to-top { font-size: 12px; }
  .back-to-top a { color: var(--blue); text-decoration: none; }
  footer { text-align: center; color: #999; font-size: 12px; padding: 30px; }
</style>
"""

html_parts = []
html_parts.append(f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>HAI Replication Regression Report</title>
{CSS}
</head><body>
<a id="top"></a>
<header>
  <h1>HAI Covariate Regression — Replication Analysis Report</h1>
  <p>Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp;|&nbsp; Model: log2_HAI ~ age + sex (OLS) &nbsp;|&nbsp; N_Obs threshold: 19</p>
</header>
""")

# ---------------- Table of Contents ----------------
subtypes_sorted = sorted(clean["subtype"].dropna().unique())
toc_items = ['<li><a href="#overview">Overview plots</a></li>']
for st in subtypes_sorted:
    toc_items.append(f'<li><a href="#subtype-{safe_name(st)}"><b>{st}</b></a></li>')
    viruses_this = sorted(clean[clean["subtype"] == st]["virus"].dropna().unique())
    for v in viruses_this:
        toc_items.append(f'<li style="margin-left:14px;">&mdash; <a href="#strain-{safe_name(st)}-{safe_name(v)}">{v}</a></li>')

html_parts.append(f"""
<nav class="toc">
  <h2>Contents</h2>
  <ul>{''.join(toc_items)}</ul>
</nav>
""")

# ---------------- Overview section ----------------
html_parts.append('<section class="block"><a id="overview"></a><h2 class="subtype-header">Overview Plots</h2>')
html_parts.append('<div class="overview-grid">')
for title, path in overview_plot_paths:
    html_parts.append(f'<figure><figcaption>{title}</figcaption>{img_tag(path, title)}</figure>')
html_parts.append('</div></section>')

# ---------------- Per Subtype -> Strain -> Day ----------------
for st in subtypes_sorted:
    html_parts.append(f'<section class="block"><a id="subtype-{safe_name(st)}"></a>')
    html_parts.append(f'<h2 class="subtype-header">Subtype: {st}</h2>')

    viruses_this = sorted(clean[clean["subtype"] == st]["virus"].dropna().unique())
    for virus_name in viruses_this:
        ann = strain_annotation_text(virus_name)
        html_parts.append(f'<a id="strain-{safe_name(st)}-{safe_name(virus_name)}"></a>')
        html_parts.append(f'<h3 class="strain-header">{virus_name} &nbsp; <span style="font-weight:400;font-size:13px;color:#555;">({ann})</span></h3>')

        days_this = sorted(clean[(clean["subtype"] == st) & (clean["virus"] == virus_name)]["day"].dropna().unique())
        for day_val in days_this:
            key = (st, virus_name, day_val)
            n_obs = len(clean[(clean["subtype"] == st) & (clean["virus"] == virus_name) & (clean["day"] == day_val)])
            vacc_status = vacc_label.get(key, "—")
            vacc_badge_cls = "vacc" if vacc_status == "Vaccinated" else "notvacc"
            has_model = key in ols_results
            model_badge = '<span class="badge ok">OLS fit</span>' if has_model else '<span class="badge skipped">Skipped (n_obs&lt;19)</span>'

            html_parts.append('<div class="day-card">')
            html_parts.append(
                f'<h4>Day {int(day_val)} '
                f'<span class="badge {vacc_badge_cls}">{vacc_status}</span> '
                f'{model_badge} '
                f'<span class="badge" style="background:#eee;color:#555;">N={n_obs}</span></h4>'
            )

            # --- images: descriptive / diagnostics / QQ ---
            html_parts.append('<div class="img-grid">')
            desc_path = desc_plot_paths.get(key)
            html_parts.append(f'<figure>{img_tag(desc_path, "Descriptive plots")}<figcaption>Descriptive (distribution / age vs HAI / sex)</figcaption></figure>')
            if has_model:
                diag_path = diag_plot_paths.get(key)
                qq_path = qq_plot_paths.get(key)
                html_parts.append(f'<figure>{img_tag(diag_path, "Diagnostic 4-panel")}<figcaption>Regression diagnostics (4-panel)</figcaption></figure>')
                html_parts.append(f'<figure>{img_tag(qq_path, "QQ plot")}<figcaption>Normal Q-Q (standalone)</figcaption></figure>')
            html_parts.append('</div>')

            # --- diagnostic stat pills ---
            if has_model and key in diag_lookup:
                d = diag_lookup[key]
                def flag(cond, text_ok, text_bad):
                    if pd.isna(cond):
                        return f'<span class="diag-pill">{text_ok.split(":")[0]}: n/a</span>'
                    cls = "ok" if cond else "flag"
                    txt = text_ok if cond else text_bad
                    return f'<span class="diag-pill badge {cls}">{txt}</span>'

                html_parts.append('<div class="diag-summary">')
                html_parts.append(f'<span class="diag-pill">Shapiro-Wilk p={d["Shapiro_p"]}</span>')
                html_parts.append(flag(d["Normal_Residuals"], "Residuals ~ Normal", "Residuals non-Normal"))
                html_parts.append(f'<span class="diag-pill">Breusch-Pagan p={d["BreuschPagan_p"]}</span>')
                html_parts.append(flag(d["Homoscedastic"], "Homoscedastic", "Heteroscedastic"))
                html_parts.append(f'<span class="diag-pill">Durbin-Watson={d["DurbinWatson"]}</span>')
                infl_cls = "flag" if d["N_Influential_CooksD"] > 0 else "ok"
                html_parts.append(f'<span class="diag-pill badge {infl_cls}">{d["N_Influential_CooksD"]} influential pts (Cook\'s D)</span>')
                html_parts.append('</div>')

            # --- regression coefficient table ---
            group_summary = summary_df[
                (summary_df["Subtype"] == st) & (summary_df["Virus"] == virus_name) & (summary_df["Day"] == day_val)
            ][["Term", "Estimate", "CI_low", "CI_high", "p_value"]]
            if len(group_summary) > 0:
                html_parts.append(df_to_html_table(group_summary))
            else:
                html_parts.append('<p class="missing">No model fit for this group (below N_Obs threshold).</p>')

            html_parts.append('<p class="back-to-top"><a href="#top">↑ back to top</a></p>')
            html_parts.append('</div>')  # day-card

    html_parts.append('</section>')

html_parts.append(f'<footer>HAI Replication Regression Report &middot; {len(ols_results)} models fit &middot; generated automatically</footer>')
html_parts.append('</body></html>')

report_path = os.path.join(OUTPUT_ROOT, "replication_hai_report.html")
with open(report_path, "w", encoding="utf-8") as f:
    f.write("\n".join(html_parts))

print(f"\nHTML report saved to: {report_path}")
print(f"Report size: {os.path.getsize(report_path) / 1e6:.1f} MB")

