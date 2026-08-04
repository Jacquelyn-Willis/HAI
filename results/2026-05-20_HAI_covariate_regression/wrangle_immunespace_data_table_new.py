
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
    "A/California/7/2009": "H1N1pdm09",
    "A/Perth/16/2009": "H3N2",
    "B/Brisbane/60/2008": "Victoria",
    "A/Indonesia/5/2005": "H5N1",
    "A/Brisbane/10/2007": "H3N2",
    "A/Victoria/361/2011": "H3N2",
    "B/Wisconsin/01/2010": "Yamagata",
    "B/Massachusetts/02/2012": "Yamagata",
    "A/Puerto Rico/8/1934": "H1N1",
    "A/Victoria/3/1975": "H3N2",
    "B/Lee/1940": "Pre-lineage",
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





#plot duistribution of values by strain and cohort of sex, age, and HAI titer values
import os
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---- config ----
STRAIN_COL = "Virus"
COHORT_COL = "Cohort"
DAY_COL = "Study Time Collected"
DAY_UNIT_COL = "Study Time Collected Unit"
SUBTYPE_COL = "subtype"
MIN_N = 10

outdir = os.path.join(SCRATCH, "hai_demo_plots_html")
os.makedirs(outdir, exist_ok=True)


def make_group_figure(df, title):
    """Build the 2x2 subplot figure for a single strain/day slice."""
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
                go.Box(y=gdf["log2_HAI"], name=str(gender), boxmean=True),
                row=2, col=1
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

toc_entries = []       # (anchor_id, display_label)
body_sections = []     # rendered HTML chunks
first_plot_written = False
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
            if SUBTYPE_COL in day_df.columns and day_df[SUBTYPE_COL].notna().any()
            else "NA"
        )
        day_unit = (
            day_df[DAY_UNIT_COL].dropna().iloc[0]
            if DAY_UNIT_COL in day_df.columns and day_df[DAY_UNIT_COL].notna().any()
            else ""
        )
        cohorts = ", ".join(sorted(day_df[COHORT_COL].dropna().unique().astype(str)))

        label = f"{strain} — Day {day} {day_unit} — Subtype: {subtype} — N={n}"
        print(f"\n=== {label} ===")
        print(day_df["log2_HAI"].describe())

        fig = make_group_figure(day_df, label)

        day_anchor = f"{strain_anchor}_day_{slugify(day)}"
        toc_entries.append((day_anchor, label))

        section_html = f"""
        <div id="{day_anchor}" class="plot-section">
            <h3>{label}</h3>
            <p class="meta">Cohort(s): {cohorts}</p>
            {fig.to_html(full_html=False, include_plotlyjs=(not first_plot_written))}
        </div>
        """
        first_plot_written = True
        strain_section_plots.append(section_html)
        total_plots += 1

    if strain_section_plots:
        body_sections.append(f'<h2 id="{strain_anchor}">{strain}</h2>')
        body_sections.extend(strain_section_plots)

# ---- assemble final HTML ----
toc_html = "<ul>" + "".join(
    f'<li><a href="#{anchor}">{label}</a></li>' for anchor, label in toc_entries
) + "</ul>"

html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>HAI Distribution Report by Strain and Day</title>
<style>
    body {{ font-family: Arial, sans-serif; margin: 40px; }}
    h1 {{ border-bottom: 2px solid #333; padding-bottom: 8px; }}
    h2 {{ margin-top: 50px; color: #1a5276; border-bottom: 1px solid #ccc; }}
    h3 {{ margin-top: 30px; color: #444; }}
    .meta {{ color: #666; font-size: 0.9em; }}
    .plot-section {{ margin-bottom: 40px; }}
    #toc {{ background: #f7f7f7; padding: 15px 25px; border-radius: 8px; margin-bottom: 30px; }}
    #toc a {{ text-decoration: none; color: #1a5276; }}
    #toc a:hover {{ text-decoration: underline; }}
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
with open(filepath, "w") as f:
    f.write(html)

print(f"\nSaved combined report: {filepath}")




#Regression models for HAI titer values with age and sex

  
    
"""

#########

"""


"""
HAI regression: auto LMM (random effect) or OLS, decided per STRAIN x DAY
==========================================================================
For EACH (virus strain, day) combination, the script automatically decides
which model to fit:

  - >=2 cohorts present for that strain-day group:
        LMM:  log2_HAI ~ age + sex + (1 | cohort)
        Cohort is a random intercept -- appropriate because there's more
        than one cohort to estimate variance across.

  - Only 1 cohort present for that strain-day group:
        OLS:  log2_HAI ~ age + sex
        Cohort can't be modeled (no variation across cohorts to estimate),
        so it's dropped entirely and only age + sex are covariates.

No model is fit comparing LMM vs OLS on the same group -- each strain-day
group gets exactly one model, chosen by its own cohort count.

VACCINATION STATUS
-------------------
Every row gets a "vaccinated" label derived from its "day" value:
    day  0  to < 14  ->  "Not Vaccinated"
    day 14  and beyond -> "Vaccinated"

Because a strain-day group is, by construction, entirely on one side of
that 14-day cutoff, "vaccinated" is CONSTANT within every group and is NOT
added as a term inside any model formula (a constant column has no
within-group variance to estimate a coefficient against, and would be
collinear with the intercept). Instead, "vaccinated" is carried through as
a LABEL on every output row/table/plot so it's easy to see, filter, and
compare Not-Vaccinated vs Vaccinated results.

HOW TO USE WITH YOUR REAL DATA
-------------------------------
This script runs on `final_merge`, which contains multiple visits per
participant/strain via the "Study Time Collected" column (values in Days,
per "Study Time Collected Unit"). Splitting by day only makes sense because
multiple days are present in this single merged frame.

If your real headers differ from the ones below (including the day column
name), run check_columns.py first and update the COLUMNS dict accordingly.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

OUTPUT_ROOT = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/hai_regression_results"
mkdirs = os.makedirs(OUTPUT_ROOT, exist_ok=True)   

# ---------------------------------------------------------------
# 1. LOAD DATA -- swap this line for your real merged (all-days) file
# ---------------------------------------------------------------
df = final_merge.copy()  # contains "Study Time Collected" (Days) across multiple visits

# ---------------------------------------------------------------
# 2. COLUMN MAP -- adjust right-hand values if your headers differ
#    NOTE: "day" must point at the column holding each row's study day
#    (the actual numeric day, e.g. 0, 3, 14, 28, 180 -- not a dataframe name).
# ---------------------------------------------------------------
df.columns = df.columns.str.strip()  # guard against stray whitespace in headers

COLUMNS = {
    "outcome": "log2_HAI",
    "age": "Age Reported_demo",
    "sex": "Gender_demo",
    "cohort": "Cohort",
    "virus": "Virus",
    "participant_id": "Participant ID",
    "day": "Study Time Collected",  # values are in Days (see "Study Time Collected Unit")
}

clean = df.rename(columns={
    COLUMNS["outcome"]: "log2_HAI",
    COLUMNS["age"]: "age",
    COLUMNS["sex"]: "sex",
    COLUMNS["cohort"]: "cohort",
    COLUMNS["virus"]: "virus",
    COLUMNS["participant_id"]: "participant_id",
    COLUMNS["day"]: "day",
})

missing = [c for c in ["log2_HAI", "age", "sex", "cohort", "virus", "day"] if c not in clean.columns]
if missing:
    raise KeyError(f"Rename didn't produce expected columns: {missing}. "
                    f"Available columns are: {df.columns.tolist()}")

# ---------------------------------------------------------------
# 3. BASIC CLEANING
# ---------------------------------------------------------------
needed = ["log2_HAI", "age", "sex", "cohort", "virus", "day"]
clean = clean.dropna(subset=needed).copy()

clean["sex"] = clean["sex"].astype("category")
clean["cohort"] = clean["cohort"].astype("category")
clean["age"] = pd.to_numeric(clean["age"], errors="coerce")
clean["day"] = pd.to_numeric(clean["day"], errors="coerce")
clean = clean.dropna(subset=["age", "day"])

# Safety check: "day" values are only comparable if they're all in the same
# unit. final_merge carries this in "Study Time Collected Unit" -- confirm
# it's "Days" everywhere before trusting the 14-day vaccination cutoff.
if "Study Time Collected Unit" in df.columns:
    bad_units = df.loc[clean.index, "Study Time Collected Unit"].dropna().unique()
    bad_units = [u for u in bad_units if str(u).strip().lower() != "days"]
    if bad_units:
        raise ValueError(
            f"'Study Time Collected' isn't uniformly in Days -- found units: {bad_units}. "
            f"Convert to days before running, or the day<14 cutoff will be wrong."
        )

# ---------------------------------------------------------------
# 3a. VACCINATION STATUS (derived from day, informational label only)
#     0  <= day < 14  -> "Not Vaccinated"
#     day >= 14        -> "Vaccinated"
# ---------------------------------------------------------------
clean["vaccinated"] = np.where(clean["day"] >= 14, "Vaccinated", "Not Vaccinated")
clean["vaccinated"] = clean["vaccinated"].astype("category")

# ---------------------------------------------------------------
# 3b. MARGINAL / CONDITIONAL R^2 FOR MIXED MODELS
#     (Nakagawa & Schielzeth 2013 method, theoretical/delta approach
#     for a linear mixed model with a single random intercept)
#
#     Marginal R^2   = Var(fixed)              / [Var(fixed) + Var(random) + Var(resid)]
#     Conditional R^2 = [Var(fixed) + Var(random)] / [Var(fixed) + Var(random) + Var(resid)]
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
# 4. PER-STRAIN-PER-DAY MODEL SELECTION + FIT
#    Decision rule (evaluated within each virus x day group):
#      n_cohorts >= 2  -> LMM with random intercept for cohort
#      n_cohorts == 1  -> OLS, age + sex only (no cohort term)
#      n_obs < 5        -> skipped
#
#    "vaccinated" is NOT in either formula (it's constant within a group,
#    since it's derived from day) -- it's recorded alongside each result
#    as a label instead.
# ---------------------------------------------------------------
lmm_results = {}       # (virus, day) -> fit
ols_results = {}       # (virus, day) -> fit
model_choice = {}      # (virus, day) -> "LMM" or "OLS"
r2_values = {}          # (virus, day) -> (marginal_r2, conditional_r2)
vacc_label = {}         # (virus, day) -> "Vaccinated" / "Not Vaccinated"

for (virus_name, day_val), sub in clean.groupby(["virus", "day"]):
    sub = sub.copy()
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)
    vacc_status = sub["vaccinated"].iloc[0]  # constant within group by construction
    key = (virus_name, day_val)
    vacc_label[key] = vacc_status

    print("=" * 70)
    print(f"Virus: {virus_name}   Day: {day_val}   [{vacc_status}]   (n_obs={n_obs}, n_cohorts={n_cohorts})")
    print("=" * 70)

    if n_obs < 20:
        print("  -> Skipped: too few observations to fit anything meaningful.\n")
        continue

    if n_cohorts >= 2:
        # ---- LMM: random intercept for cohort ----
        model_choice[key] = "LMM"
        print("  -> Decision: >=2 cohorts present -> fitting LMM (1 | cohort)")
        model = smf.mixedlm(
            "log2_HAI ~ age + sex",
            data=sub,
            groups=sub["cohort"],
            re_formula="1",
        )
        try:
            fit = model.fit(reml=True)
            lmm_results[key] = fit
            print("\n--- LMM: log2_HAI ~ age + sex + (1 | cohort) ---")
            print(fit.summary())

            marg_r2, cond_r2 = mixedlm_r2(fit)
            r2_values[key] = (marg_r2, cond_r2)
            print(f"\n  Marginal R^2 (fixed effects only):        {marg_r2:.4f}")
            print(f"  Conditional R^2 (fixed + random effects): {cond_r2:.4f}")
        except Exception as e:
            print(f"  -> LMM failed to converge: {e}")

    else:
        # ---- OLS: age + sex only, no cohort term (only 1 cohort present) ----
        model_choice[key] = "OLS"
        print("  -> Decision: only 1 cohort present -> fitting OLS (age + sex only)")
        try:
            ols_fit = smf.ols("log2_HAI ~ age + sex", data=sub).fit()
            ols_results[key] = ols_fit
            print("\n--- OLS: log2_HAI ~ age + sex ---")
            print(ols_fit.summary())

            # No random effect exists for OLS, so marginal R^2 == conditional R^2
            r2_values[key] = (ols_fit.rsquared, ols_fit.rsquared)
            print(f"\n  R^2 (= marginal = conditional, no random effect modeled): {ols_fit.rsquared:.4f}")
        except Exception as e:
            print(f"  -> OLS failed: {e}")

    print()

# ---------------------------------------------------------------
# 5. UNIFIED SUMMARY TABLE ACROSS ALL STRAIN-DAY GROUPS
#    (one row per term, regardless of which model type was used)
# ---------------------------------------------------------------
rows = []

for (virus_name, day_val), fit in lmm_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()
    for term in params.index:
        if term == "Group Var":
            continue
        rows.append({
            "Virus": virus_name, "Day": day_val, "Vaccinated": vacc_label[(virus_name, day_val)],
            "Model": "LMM", "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4),
            "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Virus": virus_name, "Day": day_val, "Vaccinated": vacc_label[(virus_name, day_val)],
        "Model": "LMM", "Term": "Cohort (random intercept) variance",
        "Estimate": round(fit.cov_re.iloc[0, 0], 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })
    rows.append({
        "Virus": virus_name, "Day": day_val, "Vaccinated": vacc_label[(virus_name, day_val)],
        "Model": "LMM", "Term": "Residual variance",
        "Estimate": round(fit.scale, 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })

for (virus_name, day_val), fit in ols_results.items():
    params, pvals, ci = fit.params, fit.pvalues, fit.conf_int()
    for term in params.index:
        rows.append({
            "Virus": virus_name, "Day": day_val, "Vaccinated": vacc_label[(virus_name, day_val)],
            "Model": "OLS", "Term": term,
            "Estimate": round(params[term], 4),
            "CI_low": round(ci.loc[term, 0], 4),
            "CI_high": round(ci.loc[term, 1], 4),
            "p_value": round(pvals[term], 4),
        })
    rows.append({
        "Virus": virus_name, "Day": day_val, "Vaccinated": vacc_label[(virus_name, day_val)],
        "Model": "OLS", "Term": "Residual variance",
        "Estimate": round(fit.mse_resid, 4),
        "CI_low": np.nan, "CI_high": np.nan, "p_value": np.nan,
    })
    # No "Cohort variance" row for OLS groups -- cohort wasn't modeled at all

summary_df = pd.DataFrame(rows)
print("=" * 70)
print("SUMMARY TABLE (all strain-day groups, model chosen automatically per group)")
print("=" * 70)
print(summary_df.to_string(index=False))


strains = final_merge[["Virus", "subtype"]]

subtype_map = strains.drop_duplicates("Virus").set_index("Virus")["subtype"]
summary_df["subtype"] = summary_df["Virus"].map(subtype_map)


summary_df.to_csv(os.path.join(OUTPUT_ROOT, "hai_regression_summary.csv"), index=False)
print(f"\nSaved: {os.path.join(OUTPUT_ROOT, 'hai_regression_summary.csv')}")

# ---------------------------------------------------------------
# 6. MODEL CHOICE LOG -- which strain-day groups got which model, and why
# ---------------------------------------------------------------
choice_rows = []
for (virus_name, day_val), sub in clean.groupby(["virus", "day"]):
    n_cohorts = sub["cohort"].nunique()
    n_obs = len(sub)
    key = (virus_name, day_val)
    marg_r2, cond_r2 = r2_values.get(key, (np.nan, np.nan))
    choice_rows.append({
        "Virus": virus_name,
        "Day": day_val,
        "Vaccinated": vacc_label.get(key, sub["vaccinated"].iloc[0]),
        "N_Obs": n_obs,
        "N_Cohorts": n_cohorts,
        "Model_Used": model_choice.get(key, "Skipped (n_obs<5)"),
        "Reason": ("Cohort modeled as random intercept (>=2 cohorts)" if n_cohorts >= 2
                   else "Cohort not modeled -- only 1 cohort present" if n_cohorts == 1
                   else "No cohort data"),
        "Marginal_R2": round(marg_r2, 4) if pd.notna(marg_r2) else np.nan,
        "Conditional_R2": round(cond_r2, 4) if pd.notna(cond_r2) else np.nan,
    })

choice_df = pd.DataFrame(choice_rows)
choice_df.to_csv(os.path.join(OUTPUT_ROOT, "hai_model_choice_log.csv"), index=False)
print("\n" + "=" * 70)
print("MODEL CHOICE LOG")
print("=" * 70)
print(choice_df.to_string(index=False))
print(f"\nSaved: {os.path.join(OUTPUT_ROOT, 'hai_model_choice_log.csv')}")


# =================================================================
# PLOTS: R^2 by strain-day group
# =================================================================
df_plot = choice_df.dropna(subset=["Model_Used"]).copy()
df_plot = df_plot[df_plot["Model_Used"] != "Skipped (n_obs<5)"]
df_plot = df_plot.sort_values("Conditional_R2", ascending=False).reset_index(drop=True)

def short_label(virus, day, maxlen=30):
    v = str(virus).replace("Influenza ", "")[:maxlen]
    return f"{v} (D{int(day)})"

df_plot["Group_short"] = df_plot.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)

# ---------------------------------------------------------------
# PLOT 1: Grouped bar chart -- Marginal vs Conditional R^2 per strain-day group
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(max(8, 1.8 * len(df_plot)), 7))

x = np.arange(len(df_plot))
width = 0.35

bars1 = ax.bar(x - width/2, df_plot["Marginal_R2"], width, label="Marginal R² (fixed effects only)",
               color="#4472C4", edgecolor="black", linewidth=0.6)
bars2 = ax.bar(x + width/2, df_plot["Conditional_R2"], width, label="Conditional R² (fixed + random effects)",
               color="#ED7D31", edgecolor="black", linewidth=0.6)

for b in bars1:
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.02, f"{b.get_height():.3f}",
            ha="center", va="bottom", fontsize=9)
for b in bars2:
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.02, f"{b.get_height():.3f}",
            ha="center", va="bottom", fontsize=9)

# annotate model type + vaccination status per group
for i, (model_used, vacc) in enumerate(zip(df_plot["Model_Used"], df_plot["Vaccinated"])):
    ax.text(i, -0.06, f"[{model_used} | {vacc}]", ha="center", va="top", fontsize=7.5,
            color="gray", style="italic")

ax.set_xlabel("Strain (Day)", fontsize=12, fontweight="bold")
ax.set_ylabel("R²", fontsize=12, fontweight="bold")
ax.set_title("Marginal vs Conditional R² by Strain-Day Group", fontsize=13, fontweight="bold", pad=15)
ax.set_xticks(x)
ax.set_xticklabels(df_plot["Group_short"], rotation=45, ha="right", fontsize=9)
ax.set_ylim([0, 1.0])
ax.legend(loc="upper right", fontsize=10)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_ROOT, "r2_marginal_vs_conditional.png"), dpi=300, bbox_inches="tight")
print("Saved: r2_marginal_vs_conditional.png")
plt.close(fig)

# ---------------------------------------------------------------
# PLOT 2: Gap plot -- how much cohort (random effect) adds, per strain-day group
# ---------------------------------------------------------------
df_plot["R2_gap"] = df_plot["Conditional_R2"] - df_plot["Marginal_R2"]
df_gap_sorted = df_plot.sort_values("R2_gap", ascending=True)

fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(df_plot))))

y_pos = np.arange(len(df_gap_sorted))
colors = ["#91bfdb" if g == 0 else "#4472C4" for g in df_gap_sorted["R2_gap"]]

bars = ax.barh(y_pos, df_gap_sorted["R2_gap"], color=colors, edgecolor="black", linewidth=0.7)

for i, (_, row) in enumerate(df_gap_sorted.iterrows()):
    label = f"{row['R2_gap']:.3f}" if row["R2_gap"] > 0 else "0 (OLS, no random effect)"
    ax.text(row["R2_gap"] + 0.01, i, label, va="center", fontsize=9)

ax.set_yticks(y_pos)
ax.set_yticklabels(df_gap_sorted["Group_short"], fontsize=9)
ax.set_xlabel("Conditional R² − Marginal R²\n(variance explained by cohort random effect)",
              fontsize=11, fontweight="bold")
ax.set_title("How Much Does Cohort (Random Effect) Add to R²?\nBy Strain-Day Group", fontsize=12, fontweight="bold", pad=15)
ax.set_xlim([0, max(0.1, df_plot["R2_gap"].max() * 1.3)])
ax.grid(axis="x", alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_ROOT, "r2_cohort_contribution_gap.png"), dpi=300, bbox_inches="tight")
print("Saved: r2_cohort_contribution_gap.png")
plt.close(fig)

print("\nDone. Groups with R2_gap = 0 used OLS (single cohort, no random effect to add).")


# =================================================================
# VARIANCE DECOMPOSITION PLOTS (per strain-day group)
# =================================================================
"""
NOTE ON LABELING
-----------------
The "residual" share is labeled "Unexplained / individual-level variance",
NOT "potential genetic signal". The model only adjusts for age, sex, and
cohort -- whatever variance is left over could reflect genetics, but could
equally reflect assay noise, batch effects, prior exposure history, timing
of sample collection, etc. The model has no way to distinguish these, so
the plot does not claim a genetic interpretation.
"""

var_rows = summary_df[summary_df["Term"].isin(
    ["Cohort (random intercept) variance", "Residual variance"]
)].copy()

wide_var = var_rows.pivot_table(
    index=["Virus", "Day", "Vaccinated"], columns="Term", values="Estimate"
).reset_index()
wide_var = wide_var.rename(columns={
    "Cohort (random intercept) variance": "Cohort_Variance",
    "Residual variance": "Residual_Variance",
})

# Groups fit with OLS have no Cohort_Variance column value -- treat as 0
if "Cohort_Variance" not in wide_var.columns:
    wide_var["Cohort_Variance"] = 0.0
wide_var["Cohort_Variance"] = wide_var["Cohort_Variance"].fillna(0.0)

wide_var["Total_Variance"] = wide_var["Cohort_Variance"] + wide_var["Residual_Variance"]
wide_var["Cohort_Variance_pct"] = 100 * wide_var["Cohort_Variance"] / wide_var["Total_Variance"]
wide_var["Residual_Variance_pct"] = 100 * wide_var["Residual_Variance"] / wide_var["Total_Variance"]

# ---------------------------------------------------------------
# Pull Std_HAI_log2 and N_Participants per strain-day group from `clean`
# ---------------------------------------------------------------
stats_by_group = (
    clean.groupby(["virus", "day"])["log2_HAI"]
    .agg(Std_HAI_log2="std", N_Participants="count")
    .reset_index()
    .rename(columns={"virus": "Virus", "day": "Day"})
)

vdf = wide_var.merge(stats_by_group, on=["Virus", "Day"], how="left")
vdf = vdf.sort_values("Residual_Variance_pct", ascending=False).reset_index(drop=True)
vdf["Group_short"] = vdf.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)

EXPLAINED_LABEL = "Explained by Age, Sex & Cohort"
RESIDUAL_LABEL = "Unexplained / individual-level variance"

# ---------------------------------------------------------------
# PLOT: Variance decomposition -- stacked bar, all strain-day groups
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(vdf)), 7))

x = np.arange(len(vdf))
explained_pct = vdf["Cohort_Variance_pct"]
residual_pct = vdf["Residual_Variance_pct"]

ax.bar(x, explained_pct, width=0.6, label="Cohort (random intercept) variance",
       color="#4472C4", alpha=0.9)
ax.bar(x, residual_pct, width=0.6, bottom=explained_pct, label=RESIDUAL_LABEL,
       color="#ED7D31", alpha=0.85)

for i, vacc in enumerate(vdf["Vaccinated"]):
    ax.text(i, -6, vacc, ha="center", va="top", fontsize=7, color="gray", style="italic", rotation=0)

ax.set_xlabel("Strain (Day)", fontsize=12, fontweight="bold")
ax.set_ylabel("Share of variance (%)", fontsize=12, fontweight="bold")
ax.set_title("HAI Titer Variance Decomposition by Strain-Day Group\n(Orange = variance not attributable to cohort, after age/sex are modeled)",
             fontsize=13, fontweight="bold", pad=20)
ax.set_xticks(x)
ax.set_xticklabels(vdf["Group_short"], rotation=45, ha="right", fontsize=9)
ax.set_ylim([0, 105])
ax.legend(loc="upper right", fontsize=10)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_ROOT, "01_variance_decomposition.png"), dpi=300, bbox_inches="tight")
print("Saved: 01_variance_decomposition.png")
plt.close(fig)

# ---------------------------------------------------------------
# PLOT: Ranked horizontal bar -- Residual Variance % by strain-day group
# ---------------------------------------------------------------
top_n = min(15, len(vdf))
top_groups = vdf.head(top_n).sort_values("Residual_Variance_pct")
y_pos = np.arange(len(top_groups))

colors = ["#d73027" if v > 60 else "#fee090" if v > 40 else "#91bfdb"
          for v in top_groups["Residual_Variance_pct"]]

fig, ax = plt.subplots(figsize=(12, max(5, 0.45 * len(top_groups))))
bars = ax.barh(y_pos, top_groups["Residual_Variance_pct"], color=colors,
                edgecolor="black", linewidth=0.8)

for i, (_, row) in enumerate(top_groups.iterrows()):
    ax.text(row["Residual_Variance_pct"] + 1, i, f"{row['Residual_Variance_pct']:.1f}%",
            va="center", fontsize=9, fontweight="bold")

ax.set_yticks(y_pos)
ax.set_yticklabels(top_groups["Group_short"], fontsize=9)
ax.set_xlabel("Unexplained / individual-level variance (%)", fontsize=11, fontweight="bold")
ax.set_title("Strain-Day Groups Ranked by Unexplained Variance\nRed >60%, Yellow 40-60%, Blue <40%",
             fontsize=12, fontweight="bold", pad=15)
ax.set_xlim([0, 100])
ax.axvline(x=60, color="red", linestyle="--", linewidth=1.5, alpha=0.5)
ax.axvline(x=40, color="orange", linestyle="--", linewidth=1.5, alpha=0.5)
ax.grid(axis="x", alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_ROOT, "02_residual_variance_ranking.png"), dpi=300, bbox_inches="tight")
print("Saved: 02_residual_variance_ranking.png")
plt.close(fig)

# ---------------------------------------------------------------
# PLOT: Scatter -- Residual Variance % vs raw titer variability
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

scatter = ax.scatter(vdf["Std_HAI_log2"], vdf["Residual_Variance_pct"],
                      s=vdf["N_Participants"] * 4, alpha=0.65,
                      c=vdf["Residual_Variance_pct"], cmap="YlOrRd",
                      edgecolors="black", linewidth=0.5)

ax.set_xlabel("HAI Titer Variability (Std of log2 titers)", fontsize=11, fontweight="bold")
ax.set_ylabel("Unexplained / individual-level variance (%)", fontsize=11, fontweight="bold")
ax.set_title("Unexplained Variance vs. Raw HAI Titer Variability\n(Bubble size = N participants, per strain-day group)",
             fontsize=12, fontweight="bold")
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Unexplained variance (%)", fontsize=10)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_ROOT, "03_residual_vs_variability.png"), dpi=300, bbox_inches="tight")
print("Saved: 03_residual_vs_variability.png")
plt.close(fig)

print("\n" + "=" * 70)
print("DONE -- 3 variance-decomposition plots generated")
print("=" * 70)


# =================================================================
# POST-RESULTS ANALYSIS: descriptive plots per strain-day group
# =================================================================
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "hai_analysis_by_group")

def safe_name(s):
    """Make a string filesystem-safe for use as a folder/file name."""
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s

def plot_group(sub_df, label, out_dir, folder_name=None):
    print(f"\n=== {label} (n={len(sub_df)}) ===")
    print(sub_df['log2_HAI'].describe())

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # log2 HAI distribution
    axes[0, 0].hist(sub_df['log2_HAI'], bins=10, alpha=0.7)
    axes[0, 0].set_xlabel("Value (log2)")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("HAI distribution")

    # age vs HI
    axes[0, 1].scatter(sub_df["age"], sub_df["log2_HAI"], alpha=0.3)
    axes[0, 1].set_xlabel("Age")
    axes[0, 1].set_ylabel("log2 HAI")
    axes[0, 1].set_title("Age vs HAI")

    # sex vs HI
    if sub_df["sex"].nunique() > 1:
        sub_df.boxplot(column="log2_HAI", by="sex", ax=axes[1, 0])
        axes[1, 0].set_title("HAI by sex")
        axes[1, 0].set_xlabel("sex")
        axes[1, 0].set_ylabel("log2 HAI")
    else:
        axes[1, 0].axis("off")

    axes[1, 1].axis("off")

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
    sub_df['log2_HAI'].describe().to_csv(stats_path)

    print(f"Saved -> {fig_path}")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- per strain-day group ---
for (strain, day_val), sub_df in clean.groupby(["virus", "day"]):
    vacc_status = sub_df["vaccinated"].iloc[0]
    label = f"Strain: {strain} | Day: {day_val} | {vacc_status}"
    folder_name = f"{safe_name(strain)}__Day{int(day_val)}__{safe_name(vacc_status)}"
    plot_group(sub_df, label, OUTPUT_DIR, folder_name=folder_name)
    
    
    
    
 """
per-group `plot_group(...)` loop finishes). It reuses `summary_df`,
`choice_df`, `OUTPUT_ROOT`, `OUTPUT_DIR`, and `clean` that the script
already built. Produces a single self-contained HTML report --
hai_regression_report.html -- with:
  - the regression term summary table (sortable, filterable by model type)
  - the model choice log table
  - R^2 plots (marginal vs conditional, cohort contribution gap)
  - variance decomposition plots (stacked bar, ranking, scatter)
  - every per-group descriptive panel (summary_plots.png), embedded inline
    as base64, with that group's log2_HAI describe() stats alongside it
"""

import base64

def img_to_base64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")

def fmt(v, nd=4):
    if pd.isna(v):
        return "—"
    return f"{v:.{nd}f}" if isinstance(v, (float, np.floating)) else str(v)

def img_tag(path, alt):
    if path and os.path.exists(path):
        return f'<img src="data:image/png;base64,{img_to_base64(path)}" alt="{alt}">'
    return f'<p class="missing">Plot not found: {os.path.basename(path) if path else alt}</p>'

# ---------------------------------------------------------------
# Regression term summary table rows
# ---------------------------------------------------------------
term_rows_html = []
for _, r in summary_df.iterrows():
    sig = pd.notna(r["p_value"]) and r["p_value"] < 0.05
    row_class = "sig" if sig else ""
    term_rows_html.append(f"""
    <tr class="{row_class}" data-model="{r['Model']}">
        <td>{r['Virus']}</td>
        <td>{r['Day']}</td>
        <td>{r['Vaccinated']}</td>
        <td>{r['Model']}</td>
        <td>{r['Term']}</td>
        <td>{fmt(r['Estimate'])}</td>
        <td>{fmt(r['CI_low'])}</td>
        <td>{fmt(r['CI_high'])}</td>
        <td>{fmt(r['p_value'])}</td>
    </tr>""")

# ---------------------------------------------------------------
# Model choice log table rows
# ---------------------------------------------------------------
choice_rows_html = []
for _, r in choice_df.iterrows():
    model_class = {"LMM": "lmm", "OLS": "ols"}.get(r["Model_Used"], "skipped")
    choice_rows_html.append(f"""
    <tr class="{model_class}">
        <td>{r['Virus']}</td>
        <td>{r['Day']}</td>
        <td>{r['Vaccinated']}</td>
        <td>{r['N_Obs']}</td>
        <td>{r['N_Cohorts']}</td>
        <td>{r['Model_Used']}</td>
        <td>{r['Reason']}</td>
        <td>{fmt(r['Marginal_R2'])}</td>
        <td>{fmt(r['Conditional_R2'])}</td>
    </tr>""")

# ---------------------------------------------------------------
# Per-group descriptive panels (summary_plots.png + summary_stats.csv)
# ---------------------------------------------------------------
def safe_name(s):
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s

group_panels_html = []
for (strain, day_val), sub_df in clean.groupby(["virus", "day"]):
    vacc_status = sub_df["vaccinated"].iloc[0]
    folder_name = f"{safe_name(strain)}__Day{int(day_val)}__{safe_name(vacc_status)}"
    group_folder = os.path.join(OUTPUT_DIR, folder_name)
    fig_path = os.path.join(group_folder, "summary_plots.png")

    desc = sub_df["log2_HAI"].describe()
    stats_line = " | ".join(f"{k}={v:.3f}" for k, v in desc.items())

    group_panels_html.append(f"""
    <div class="panel-card">
        <h3>{strain} &middot; Day {int(day_val)} &middot; {vacc_status} <span class="tag">n={len(sub_df)}</span></h3>
        <p class="stats-line">{stats_line}</p>
        {img_tag(fig_path, f"{strain} day {day_val}")}
    </div>""")

# ---------------------------------------------------------------
# Top-level plot paths
# ---------------------------------------------------------------
r2_bar_path = os.path.join(OUTPUT_ROOT, "r2_marginal_vs_conditional.png")
r2_gap_path = os.path.join(OUTPUT_ROOT, "r2_cohort_contribution_gap.png")
var_decomp_path = os.path.join(OUTPUT_ROOT, "01_variance_decomposition.png")
var_rank_path = os.path.join(OUTPUT_ROOT, "02_residual_variance_ranking.png")
var_scatter_path = os.path.join(OUTPUT_ROOT, "03_residual_vs_variability.png")

# ---------------------------------------------------------------
# Assemble final HTML
# ---------------------------------------------------------------
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HAI Regression Report</title>
<style>
    :root {{
        --sig: #fee090; --lmm: #4472C4; --ols: #91bfdb; --skipped: #d9d9d9;
        --bg: #f7f8fa; --card-bg: #ffffff; --border: #dfe3e8; --text: #1f2937; --muted: #6b7280;
    }}
    * {{ box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        background: var(--bg); color: var(--text); margin: 0; padding: 0 0 60px 0;
    }}
    header {{ background: #1f2937; color: white; padding: 28px 40px; }}
    header h1 {{ margin: 0 0 4px 0; font-size: 22px; }}
    header p {{ margin: 0; color: #b7c0cc; font-size: 14px; }}
    nav {{
        position: sticky; top: 0; z-index: 10; background: white; border-bottom: 1px solid var(--border);
        padding: 10px 40px; display: flex; gap: 20px; font-size: 14px; flex-wrap: wrap;
    }}
    nav a {{ color: #2563eb; text-decoration: none; font-weight: 500; }}
    nav a:hover {{ text-decoration: underline; }}
    section {{ padding: 30px 40px; }}
    h2 {{ font-size: 18px; border-bottom: 2px solid var(--border); padding-bottom: 8px; }}
    .legend {{ display: flex; gap: 18px; font-size: 13px; margin-bottom: 14px; color: var(--muted); flex-wrap: wrap; }}
    .swatch {{ display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 6px; vertical-align: middle; }}
    table {{
        border-collapse: collapse; width: 100%; background: var(--card-bg);
        font-size: 13px; box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    th, td {{ padding: 8px 10px; border-bottom: 1px solid var(--border); text-align: right; white-space: nowrap; }}
    th {{ text-align: right; background: #f0f2f5; cursor: pointer; user-select: none; position: sticky; top: 45px; }}
    td:first-child, th:first-child, td:nth-child(3), th:nth-child(3),
    td:nth-child(5), th:nth-child(5), td:nth-child(7), th:nth-child(7) {{ text-align: left; }}
    tr.sig {{ background: rgba(254,224,144,0.35); }}
    tr.lmm {{ background: rgba(68,114,196,0.10); }}
    tr.ols {{ background: rgba(145,191,219,0.15); }}
    tr.skipped {{ background: rgba(217,217,217,0.3); color: var(--muted); }}
    th:after {{ content: " \\2195"; color: #9ca3af; font-size: 11px; }}
    .filter-bar {{ margin-bottom: 10px; font-size: 13px; }}
    .filter-bar button {{
        padding: 5px 12px; margin-right: 6px; border: 1px solid var(--border); border-radius: 6px;
        background: white; cursor: pointer; font-size: 12px;
    }}
    .filter-bar button.active {{ background: #2563eb; color: white; border-color: #2563eb; }}
    .panel-card {{
        background: var(--card-bg); border: 1px solid var(--border); border-radius: 8px;
        padding: 16px 20px; margin-bottom: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    .panel-card h3 {{ margin: 0 0 4px 0; font-size: 15px; }}
    .panel-card .tag {{
        font-size: 11px; font-weight: 600; background: #eef2ff; color: #4338ca;
        padding: 2px 8px; border-radius: 10px; margin-left: 6px;
    }}
    .stats-line {{ font-size: 12px; color: var(--muted); margin: 0 0 10px 0; word-break: break-word; }}
    .panel-card img, .plot-block img {{ width: 100%; max-width: 900px; border: 1px solid var(--border); border-radius: 4px; }}
    .plot-block {{ margin-bottom: 30px; }}
    .plot-block h3 {{ font-size: 14px; margin-bottom: 8px; }}
    .missing {{ color: #d73027; font-size: 12px; font-style: italic; }}
</style>
</head>
<body>

<header>
    <h1>HAI Regression Report</h1>
    <p>{clean['virus'].nunique()} strains &middot; {summary_df[['Virus','Day']].drop_duplicates().shape[0]} strain-day groups modeled &middot;
       auto LMM (&ge;2 cohorts) / OLS (1 cohort)</p>
</header>

<nav>
    <a href="#terms">Regression Terms</a>
    <a href="#choice">Model Choice Log</a>
    <a href="#r2">R&sup2; Plots</a>
    <a href="#variance">Variance Decomposition</a>
    <a href="#groups">Per-Group Descriptives</a>
</nav>

<section id="terms">
    <h2>Regression Term Summary</h2>
    <div class="legend">
        <span><span class="swatch" style="background:var(--sig)"></span>p &lt; 0.05</span>
    </div>
    <div class="filter-bar">
        Filter: <button class="active" onclick="filterModel('all', this)">All</button>
        <button onclick="filterModel('LMM', this)">LMM only</button>
        <button onclick="filterModel('OLS', this)">OLS only</button>
    </div>
    <table id="term-table">
        <thead>
            <tr>
                <th>Virus</th><th>Day</th><th>Vaccinated</th><th>Model</th><th>Term</th>
                <th>Estimate</th><th>CI low</th><th>CI high</th><th>p-value</th>
            </tr>
        </thead>
        <tbody>
            {''.join(term_rows_html)}
        </tbody>
    </table>
</section>

<section id="choice">
    <h2>Model Choice Log</h2>
    <div class="legend">
        <span><span class="swatch" style="background:var(--lmm)"></span>LMM (&ge;2 cohorts)</span>
        <span><span class="swatch" style="background:var(--ols)"></span>OLS (1 cohort)</span>
        <span><span class="swatch" style="background:var(--skipped)"></span>Skipped</span>
    </div>
    <table id="choice-table">
        <thead>
            <tr>
                <th>Virus</th><th>Day</th><th>Vaccinated</th><th>N Obs</th><th>N Cohorts</th>
                <th>Model Used</th><th>Reason</th><th>Marg. R&sup2;</th><th>Cond. R&sup2;</th>
            </tr>
        </thead>
        <tbody>
            {''.join(choice_rows_html)}
        </tbody>
    </table>
</section>

<section id="r2">
    <h2>R&sup2; Plots</h2>
    <div class="plot-block">
        <h3>Marginal vs Conditional R&sup2; by Strain-Day Group</h3>
        {img_tag(r2_bar_path, "R2 marginal vs conditional")}
    </div>
    <div class="plot-block">
        <h3>Cohort (Random Effect) Contribution Gap</h3>
        {img_tag(r2_gap_path, "R2 cohort contribution gap")}
    </div>
</section>

<section id="variance">
    <h2>Variance Decomposition</h2>
    <div class="plot-block">
        <h3>Variance Decomposition by Strain-Day Group</h3>
        {img_tag(var_decomp_path, "variance decomposition")}
    </div>
    <div class="plot-block">
        <h3>Ranked by Unexplained Variance</h3>
        {img_tag(var_rank_path, "residual variance ranking")}
    </div>
    <div class="plot-block">
        <h3>Unexplained Variance vs Raw Titer Variability</h3>
        {img_tag(var_scatter_path, "residual vs variability scatter")}
    </div>
</section>

<section id="groups">
    <h2>Per-Group Descriptive Panels</h2>
    {''.join(group_panels_html)}
</section>

<script>
document.querySelectorAll("table").forEach(table => {{
    table.querySelectorAll("th").forEach((th, idx) => {{
        th.addEventListener("click", () => {{
            const tbody = table.querySelector("tbody");
            const rows = Array.from(tbody.querySelectorAll("tr"));
            const asc = th.dataset.asc !== "true";
            table.querySelectorAll("th").forEach(h => h.dataset.asc = "");
            th.dataset.asc = asc;
            rows.sort((a, b) => {{
                let x = a.children[idx].innerText.trim();
                let y = b.children[idx].innerText.trim();
                const nx = parseFloat(x), ny = parseFloat(y);
                if (!isNaN(nx) && !isNaN(ny)) {{ x = nx; y = ny; }}
                if (x < y) return asc ? -1 : 1;
                if (x > y) return asc ? 1 : -1;
                return 0;
            }});
            rows.forEach(r => tbody.appendChild(r));
        }});
    }});
}});

function filterModel(model, btn) {{
    document.querySelectorAll(".filter-bar button").forEach(b => b.classList.remove("active"));
    btn.classList.add("active");
    document.querySelectorAll("#term-table tbody tr").forEach(row => {{
        row.style.display = (model === "all" || row.dataset.model === model) ? "" : "none";
    }});
}}
</script>

</body>
</html>
"""

report_path = os.path.join(OUTPUT_ROOT, "hai_regression_report.html")
with open(report_path, "w", encoding="utf-8") as f:
    f.write(html)

print(f"Saved: {report_path}")
print("Open this single file in a browser -- all images are embedded, no other files needed.")   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

"""
HAI model diagnostics: QQ plots, residual checks, and performance metrics
============================================================================
Run this AFTER your per-strain-per-day regression script in the same
session -- it reuses `clean`, `lmm_results`, `ols_results`, `model_choice`,
`r2_values`, and `vacc_label` that script already built. Nothing here
refits models; it only diagnoses the ones already fit.

Models there are keyed by (virus, day) tuples, e.g.:
    lmm_results[("A/South Dakota/06/2007", 0.0)]

For EACH (strain, day) group that got a model (LMM or OLS), this produces
a 2x2 diagnostic panel:
    1. Residuals vs Fitted   -- look for curvature (non-linearity) or a
                                 funnel shape (heteroscedasticity)
    2. Normal Q-Q            -- points should hug the diagonal; systematic
                                 curvature = non-normal residuals
    3. Scale-Location         -- sqrt(|standardized resid|) vs fitted;
                                 should be roughly flat (constant variance)
    4. Histogram of residuals -- with a fitted normal curve overlay

It also builds ONE summary table (model_performance_summary.csv) across
all strain-day groups with AIC, BIC, RMSE, MAE, R^2, and a Shapiro-Wilk
normality test (W statistic + p-value) on each model's residuals.

NOTE ON "standardized residuals" for LMM
------------------------------------------
Proper studentized/leverage-adjusted residuals require a hat matrix, which
isn't well-defined the same way for mixed models. Here, residuals are
standardized simply as (resid - mean) / std(resid). This is fine for
visually checking variance patterns and normality, but isn't a leverage
diagnostic -- there's no per-group Cook's distance / leverage panel here.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

DIAG_DIR = os.path.join(OUTPUT_ROOT, "model_diagnostics")
os.makedirs(DIAG_DIR, exist_ok=True)

def safe_name(s):
    s = str(s).strip()
    s = re.sub(r'[\\/*?:"<>|]', "_", s)
    s = re.sub(r'\s+', "_", s)
    return s

def get_resid_fitted(fit):
    """Works for both OLS and MixedLM results."""
    resid = np.asarray(fit.resid)
    fitted = np.asarray(fit.fittedvalues)
    return resid, fitted

def diagnostic_plots(virus_name, day_val, vacc_status, fit, model_type, out_dir):
    resid, fitted = get_resid_fitted(fit)
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

    # 4. Histogram of residuals + normal curve
    ax = axes[1, 1]
    ax.hist(resid, bins=15, density=True, alpha=0.7, color="#4472C4", edgecolor="black")
    xs = np.linspace(resid.min(), resid.max(), 200)
    ax.plot(xs, stats.norm.pdf(xs, resid.mean(), resid.std(ddof=1)), color="red", linewidth=2)
    ax.set_xlabel("Residuals")
    ax.set_ylabel("Density")
    ax.set_title("Residual Distribution")
    ax.grid(alpha=0.3)

    plt.suptitle(f"{virus_name}  |  Day {day_val}  |  {vacc_status}  [{model_type}]",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()

    folder_name = f"{safe_name(virus_name)}__Day{day_val}__{safe_name(vacc_status)}"
    folder = os.path.join(out_dir, folder_name)
    os.makedirs(folder, exist_ok=True)
    fig_path = os.path.join(folder, "diagnostics.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return fig_path

# ---------------------------------------------------------------
# Run diagnostics for every fitted model, LMM and OLS alike
# keys are (virus_name, day_val) tuples
# ---------------------------------------------------------------
perf_rows = []

all_fits = [(key, f, "LMM") for key, f in lmm_results.items()] + \
           [(key, f, "OLS") for key, f in ols_results.items()]

for (virus_name, day_val), fit, model_type in all_fits:
    vacc_status = vacc_label.get((virus_name, day_val), "Unknown")

    resid, fitted = get_resid_fitted(fit)
    n_obs = int(fit.nobs) if hasattr(fit, "nobs") else len(resid)
    n_params = len(fit.params)

    rmse = float(np.sqrt(np.mean(resid ** 2)))
    mae = float(np.mean(np.abs(resid)))
    aic = getattr(fit, "aic", np.nan)
    bic = getattr(fit, "bic", np.nan)

    # Shapiro-Wilk normality test (unreliable above ~5000 obs, unusable below 3)
    if 3 <= n_obs <= 5000:
        try:
            sw_stat, sw_p = stats.shapiro(resid)
        except Exception:
            sw_stat, sw_p = np.nan, np.nan
    else:
        sw_stat, sw_p = np.nan, np.nan

    marg_r2, cond_r2 = r2_values.get((virus_name, day_val), (np.nan, np.nan))

    fig_path = diagnostic_plots(virus_name, day_val, vacc_status, fit, model_type, DIAG_DIR)

    perf_rows.append({
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
    print(f"Diagnostics done: {virus_name} | Day {day_val} | {vacc_status} [{model_type}]  "
          f"RMSE={rmse:.3f}  AIC={aic if pd.notna(aic) else 'NA'}  "
          f"Shapiro_p={sw_p if pd.notna(sw_p) else 'NA'}")

perf_df = pd.DataFrame(perf_rows).sort_values(["Model", "AIC"]).reset_index(drop=True)
perf_df.to_csv(os.path.join(OUTPUT_ROOT, "model_performance_summary.csv"), index=False)

print("\n" + "=" * 70)
print("MODEL PERFORMANCE SUMMARY")
print("=" * 70)
print(perf_df.to_string(index=False))
print(f"\nSaved: {os.path.join(OUTPUT_ROOT, 'model_performance_summary.csv')}")


# ---------------------------------------------------------------
# Overview plot: RMSE per strain-day group, colored by whether residuals
# passed the Shapiro-Wilk normality check
# ---------------------------------------------------------------
def short_label(virus, day, maxlen=30):
    v = str(virus).replace("Influenza ", "")[:maxlen]
    return f"{v} (D{day})"

perf_df["Group_short"] = perf_df.apply(lambda r: short_label(r["Virus"], r["Day"]), axis=1)
perf_sorted = perf_df.sort_values("RMSE")



colors = ["#91bfdb" if v else "#d73027" for v in perf_sorted["Residuals_Normal_at_0.05"].fillna(False)]

fig, ax = plt.subplots(figsize=(10, max(4, 0.4 * len(perf_sorted))))
y_pos = np.arange(len(perf_sorted))
ax.barh(y_pos, perf_sorted["RMSE"], color=colors, edgecolor="black", linewidth=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(perf_sorted["Group_short"], fontsize=9)
ax.set_xlabel("RMSE (log2 HAI units)", fontsize=11, fontweight="bold")
ax.set_title("Model RMSE by Strain-Day Group\nBlue = residuals pass Shapiro-Wilk normality (p>0.05), Red = fail",
             fontsize=12, fontweight="bold", pad=15)
ax.grid(axis="x", alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_ROOT, "05_rmse_by_group.png"), dpi=300, bbox_inches="tight")
print("Saved: 05_rmse_by_group.png")
plt.close(fig)

print("\nDone. Per-group 2x2 diagnostic panels saved under:")
print(f"  {DIAG_DIR}/<strain>__Day<day>__<vaccinated>/diagnostics.png")






"""
RMSE overview plot are built). It produces a single self-contained HTML
report -- model_diagnostics_report.html -- with:
  - the summary table (sortable via click, color-coded by Shapiro-Wilk pass/fail)
  - the RMSE overview bar chart
  - every 2x2 diagnostic panel, embedded inline as base64 (no external image
    files needed to view it -- you can email/share the single .html file)
"""

import base64

def img_to_base64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")

def fmt(v, nd=4):
    if pd.isna(v):
        return "—"
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)

# ---------------------------------------------------------------
# Build summary table rows
# ---------------------------------------------------------------
table_rows_html = []
for _, r in perf_df.iterrows():
    passed = r["Residuals_Normal_at_0.05"]
    if pd.isna(passed):
        row_class = "unknown"
    elif passed:
        row_class = "pass"
    else:
        row_class = "fail"

    table_rows_html.append(f"""
    <tr class="{row_class}">
        <td>{r['Virus']}</td>
        <td>{r['Day']}</td>
        <td>{r['Vaccinated']}</td>
        <td>{r['Model']}</td>
        <td>{r['N_Obs']}</td>
        <td>{fmt(r['AIC'], 2)}</td>
        <td>{fmt(r['BIC'], 2)}</td>
        <td>{fmt(r['RMSE'])}</td>
        <td>{fmt(r['MAE'])}</td>
        <td>{fmt(r['Marginal_R2'])}</td>
        <td>{fmt(r['Conditional_R2'])}</td>
        <td>{fmt(r['Shapiro_W'])}</td>
        <td>{fmt(r['Shapiro_p'])}</td>
        <td>{'Yes' if passed is True else ('No' if passed is False else '—')}</td>
    </tr>""")

# ---------------------------------------------------------------
# Build per-group diagnostic panel sections, embedded as base64
# ---------------------------------------------------------------
panel_sections_html = []
for _, r in perf_df.sort_values(["Virus", "Day"]).iterrows():
    b64 = img_to_base64(r["Diagnostic_Plot"])
    panel_sections_html.append(f"""
    <div class="panel-card">
        <h3>{r['Virus']} &middot; Day {r['Day']} &middot; {r['Vaccinated']} <span class="tag">{r['Model']}</span></h3>
        <p class="stats-line">
            RMSE={fmt(r['RMSE'])} | MAE={fmt(r['MAE'])} | AIC={fmt(r['AIC'],2)} | BIC={fmt(r['BIC'],2)} |
            Shapiro p={fmt(r['Shapiro_p'])} ({'normal' if r['Residuals_Normal_at_0.05'] is True else 'non-normal' if r['Residuals_Normal_at_0.05'] is False else 'n/a'})
        </p>
        <img src="data:image/png;base64,{b64}" alt="diagnostics for {r['Virus']} day {r['Day']}">
    </div>""")

rmse_overview_path = os.path.join(OUTPUT_ROOT, "05_rmse_by_group.png")
rmse_b64 = img_to_base64(rmse_overview_path) if os.path.exists(rmse_overview_path) else None

# ---------------------------------------------------------------
# Assemble final HTML
# ---------------------------------------------------------------
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HAI Model Diagnostics Report</title>
<style>
    :root {{
        --pass: #91bfdb;
        --fail: #d73027;
        --bg: #f7f8fa;
        --card-bg: #ffffff;
        --border: #dfe3e8;
        --text: #1f2937;
        --muted: #6b7280;
    }}
    * {{ box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        background: var(--bg);
        color: var(--text);
        margin: 0;
        padding: 0 0 60px 0;
    }}
    header {{
        background: #1f2937;
        color: white;
        padding: 28px 40px;
    }}
    header h1 {{ margin: 0 0 4px 0; font-size: 22px; }}
    header p {{ margin: 0; color: #b7c0cc; font-size: 14px; }}
    nav {{
        position: sticky; top: 0; z-index: 10;
        background: white; border-bottom: 1px solid var(--border);
        padding: 10px 40px; display: flex; gap: 20px; font-size: 14px;
    }}
    nav a {{ color: #2563eb; text-decoration: none; font-weight: 500; }}
    nav a:hover {{ text-decoration: underline; }}
    section {{ padding: 30px 40px; }}
    h2 {{ font-size: 18px; border-bottom: 2px solid var(--border); padding-bottom: 8px; }}
    .legend {{ display: flex; gap: 18px; font-size: 13px; margin-bottom: 14px; color: var(--muted); }}
    .swatch {{ display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 6px; vertical-align: middle; }}
    table {{
        border-collapse: collapse; width: 100%; background: var(--card-bg);
        font-size: 13px; box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    th, td {{ padding: 8px 10px; border-bottom: 1px solid var(--border); text-align: right; white-space: nowrap; }}
    th {{ text-align: right; background: #f0f2f5; cursor: pointer; user-select: none; position: sticky; top: 45px; }}
    td:first-child, th:first-child {{ text-align: left; }}
    td:nth-child(3), th:nth-child(3) {{ text-align: left; }}
    td:nth-child(4), th:nth-child(4) {{ text-align: left; }}
    tr.pass {{ background: rgba(145,191,219,0.15); }}
    tr.fail {{ background: rgba(215,48,39,0.10); }}
    th:after {{ content: " \\2195"; color: #9ca3af; font-size: 11px; }}
    .panel-card {{
        background: var(--card-bg); border: 1px solid var(--border); border-radius: 8px;
        padding: 16px 20px; margin-bottom: 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    .panel-card h3 {{ margin: 0 0 4px 0; font-size: 15px; }}
    .panel-card .tag {{
        font-size: 11px; font-weight: 600; background: #eef2ff; color: #4338ca;
        padding: 2px 8px; border-radius: 10px; margin-left: 6px;
    }}
    .stats-line {{ font-size: 12px; color: var(--muted); margin: 0 0 10px 0; }}
    .panel-card img {{ width: 100%; max-width: 900px; border: 1px solid var(--border); border-radius: 4px; }}
    #overview img {{ max-width: 900px; width: 100%; border: 1px solid var(--border); border-radius: 4px; }}
</style>
</head>
<body>

<header>
    <h1>HAI Model Diagnostics Report</h1>
    <p>{len(perf_df)} strain-day models &middot; generated from per-strain-per-day LMM/OLS fits</p>
</header>

<nav>
    <a href="#summary">Summary Table</a>
    <a href="#overview">RMSE Overview</a>
    <a href="#panels">Diagnostic Panels</a>
</nav>

<section id="summary">
    <h2>Model Performance Summary</h2>
    <div class="legend">
        <span><span class="swatch" style="background:var(--pass)"></span>Residuals normal (Shapiro p &gt; 0.05)</span>
        <span><span class="swatch" style="background:var(--fail)"></span>Residuals non-normal (p &le; 0.05)</span>
    </div>
    <table id="perf-table">
        <thead>
            <tr>
                <th>Virus</th><th>Day</th><th>Vaccinated</th><th>Model</th><th>N</th>
                <th>AIC</th><th>BIC</th><th>RMSE</th><th>MAE</th>
                <th>Marg. R²</th><th>Cond. R²</th><th>Shapiro W</th><th>Shapiro p</th><th>Normal?</th>
            </tr>
        </thead>
        <tbody>
            {''.join(table_rows_html)}
        </tbody>
    </table>
</section>

<section id="overview">
    <h2>RMSE by Strain-Day Group</h2>
    {f'<img src="data:image/png;base64,{rmse_b64}" alt="RMSE overview">' if rmse_b64 else '<p>Overview plot not found.</p>'}
</section>

<section id="panels">
    <h2>Per-Group Diagnostic Panels</h2>
    {''.join(panel_sections_html)}
</section>

<script>
// Simple click-to-sort on the summary table
document.querySelectorAll("#perf-table th").forEach((th, idx) => {{
    th.addEventListener("click", () => {{
        const table = th.closest("table");
        const tbody = table.querySelector("tbody");
        const rows = Array.from(tbody.querySelectorAll("tr"));
        const asc = th.dataset.asc !== "true";
        table.querySelectorAll("th").forEach(h => h.dataset.asc = "");
        th.dataset.asc = asc;
        rows.sort((a, b) => {{
            let x = a.children[idx].innerText.trim();
            let y = b.children[idx].innerText.trim();
            const nx = parseFloat(x), ny = parseFloat(y);
            if (!isNaN(nx) && !isNaN(ny)) {{ x = nx; y = ny; }}
            if (x < y) return asc ? -1 : 1;
            if (x > y) return asc ? 1 : -1;
            return 0;
        }});
        rows.forEach(r => tbody.appendChild(r));
    }});
}});
</script>

</body>
</html>
"""

report_path = os.path.join(OUTPUT_ROOT, "model_diagnostics_report.html")
with open(report_path, "w", encoding="utf-8") as f:
    f.write(html)

print(f"Saved: {report_path}")
print("Open this single file in a browser -- all images are embedded, no other files needed.")




#####choosing a strain 

#summary_df 

strains = final_merge[["Virus", "subtype"]]

subtype_map = strains.drop_duplicates("Virus").set_index("Virus")["subtype"]
summary_df["subtype"] = summary_df["Virus"].map(subtype_map)

