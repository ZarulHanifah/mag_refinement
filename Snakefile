
import sys , os , glob , re
import pandas as pd
import numpy as np
from pathlib import Path

from snakemake.exceptions import WorkflowError

from magrefine.mags import Mag
from magrefine.sessionmanager import SessionManager

if not os.path.exists("logs/cluster"):
    os.makedirs("logs/cluster")

SUMMARY_PATH = Path(config["summary_path"])
MAG_DIR = Path(config["mag_dir"])
ABUND_DIR  = Path(config["abund_dir"])

sesh = SessionManager(
    summary_path=SUMMARY_PATH,
    mag_dir=MAG_DIR,
    abund_dir=ABUND_DIR
)

# configfile: "config_snakefile.yaml"

results_path    = config["results_path"]
temp_path       = config["temp_path"]
samples_dorado  = config["dorado7"]
samples_list    = sorted(list(samples_dorado.keys()))

uniref_db       = config["uniref_db"]
checkm1_db      = config["checkm1_db"]
checkm2_db      = config["checkm2_db"]

medaka_model    = "r1041_e82_400bps_sup_v5.0.0"

# --- MODIFIED SECTION: MAG SELECTION LOGIC ---
if "mags_to_process" in config:
    # If 'mags_to_process' is provided via the config (e.g., from execute_workflow.py)
    mags = config["mags_to_process"]
else:
    raise WorkflowError("The 'mags_to_process' list must be provided via Snakemake config (from execute_workflow.py)")
# --- END MODIFIED SECTION ---

include: "rules/herro.smk"
include: "rules/assembly.smk"
include: "rules/coverm.smk"
include: "rules/ntlink.smk"
include: "rules/checkm2.smk"

wildcard_constraints:
    mag="[a-zA-Z0-9_\.]+"

rule all:
    input:
        expand(rules.flye_fq.output.assem, mag=mags),
        expand(rules.hifiasm_fq.output.assem, mag=mags),
        expand(rules.myloasm_fq.output.assem, mag=mags),
        expand(rules.longstitch.output, mag=mags),
        expand(rules.checkm2_original.output.report, mag=mags),
        expand(rules.checkm2_flye_fq.output.report, mag=mags),
        expand(rules.checkm2_hifiasm_fq.output.report, mag=mags),
        expand(rules.checkm2_myloasm_fq.output.report, mag=mags),
        expand(rules.checkm2_longstitch.output.report, mag=mags),
        expand(rules.checkm1_original.output.report, mag=mags),
        expand(rules.magpurify2_coverage.output, mag=mags)

