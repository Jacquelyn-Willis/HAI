import pandas as pd 
import os
import numpy as np

#directories 

DATA = "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/"

data = "/sc/arion/work/willij115/projects/HAI/data/2026-05-20_HAI_covariate_regression"
scratch = "/sc/arion/scratch/willij115/projects/HAI/2026-05-20_HAI_covariate_regression"
results = "/sc/arion/work/willij115/projects/HAI/results/2026-05-20_HAI_covariate_regression"

#upload immunespace data tables
 
studies = pd.read_csv(os.path.join(DATA,'immunespaceHAI_studies_tables.csv'), sep = ',', header = 0)
arms = pd.read_csv(os.path.join(DATA,'immunespaceHAI_arms_tables.csv'), sep = ',', header = 0)
participants = pd.read_csv(os.path.join(DATA,'immunespaceHAI_participants_tables.csv'), sep = ',', header = 0)
events = pd.read_csv(os.path.join(DATA,'immunespaceHAI_events_tables.csv'), sep = ',', header = 0)
assays = pd.read_csv(os.path.join(DATA,'immunespaceHAI_assays_tables.csv'), sep = ',', header = 0)