# Summary: A Manifest-Driven Workflow for MAG Refinement

This document outlines a proposed software design to evolve the current MAG processing pipeline into a more flexible, scalable, and robust system for iterative MAG refinement.

## 1. The Core Problem

The current system has a **tight coupling** between a MAG's identity and its file path. The MAG name (e.g., `C1A3_A_metabat.872`) is used directly to find its corresponding FASTA file (`C1A3_A_metabat.872.fasta`) and to derive the path to its abundance data.

This works for a single round of analysis but becomes brittle when we consider re-assembly and refinement:
- **Naming is Restrictive**: We cannot give a memorable, scientific name (e.g., `Asg1`) to a MAG without breaking the file-finding logic.
- **Tracking is Difficult**: If we re-assemble `Asg1` with two different tools (e.g., Flye and Hifiasm), creating `Asg1_flye` and `Asg1_hifiasm`, the relationship between the parent and its children is not formally tracked.
- **File Paths are Inflexible**: All MAGs must reside in a single directory.

## 2. The Solution: A Manifest-Driven Workflow

The proposed solution is to **decouple** the MAG's logical identity from its physical storage by introducing a central **manifest file**: `mags_manifest.tsv`.

This file will act as the **single source of truth** for the entire workflow. Instead of scanning directories and relying on naming conventions, our scripts will now read the manifest to discover MAGs, their properties, and their file locations.

## 3. The Manifest File Specification

The `mags_manifest.tsv` will be a tab-separated file with the following columns:

| Column Name | Purpose | Example |
| :--- | :--- | :--- |
| `mag_id` | **Primary Key**. The unique, stable ID for a MAG. Can be user-defined (e.g., `Asg1`) or default to the original filename. **MUST NOT be empty.** | `Asg1` |
| `parent_id` | The `mag_id` of the parent MAG. This is the key to tracking lineage. Empty for initial MAGs. | `Asg1` |
| `fasta_path` | The relative or absolute path to the MAG's FASTA file. | `data/original_mags/C1A3_A.fasta` |
| `abundance_path`| The path to the abundance/depth file this MAG should be compared against. | `data/abundance/myloasm__C001_A3.tsv` |
| `summary_name` | The ID used to find the original stats (completeness, etc.) in the initial `summary.tsv` file. | `C1A3_A_metabat.872` |

### Addressing Key Concerns:
1.  **"What if I'm not ready to assign a custom ID?"**
    *   The `mag_id` will simply default to the original filename stem (e.g., `C1A3_A_metabat.872`). This allows the system to work immediately without requiring any manual naming. You can update the `mag_id` to a custom name later.
2.  **"I can't change a mag_id halfway through."**
    *   This is correct. Once a `mag_id` has been used as a `parent_id` for child MAGs, it should be considered **immutable**. It is the primary key that links the family tree together.

## 4. The New Workflow Lifecycle

The workflow becomes a cyclical process of refinement, where the manifest is both an input and an output.

### Step 1: Bootstrap
A helper script, `create_manifest_from_directory.py`, is run once to generate the initial `mags_manifest.tsv` from the existing file structure. This makes adoption painless.

### Step 2: Selection & Execution
`execute_workflow.py` is modified to read from the manifest. The user selects MAGs to process based on their `mag_id`. The script then uses the `fasta_path` and `abundance_path` columns to run the Snakemake workflow.

### Step 3: Refinement & Analysis
The Snakemake workflow performs re-assembly, generating new FASTA files (e.g., `results/reassembly/Asg1/flye.fasta`). It then runs quality control (CheckM2, etc.) on these new assemblies.

### Step 4: Augmentation
A final step in the workflow appends the results of the refinement back to the manifest file as new rows. For example:

| mag_id | parent_id | fasta_path | abundance_path | summary_name |
| :--- | :--- | :--- | :--- | :--- |
| ... | ... | ... | ... | ... |
| `Asg1_flye`| `Asg1` | `results/Asg1/flye.fasta` | `data/abund/myloasm_C001_A3.tsv`| `Asg1_flye` |

The `abundance_path` is inherited from the parent, as the underlying read data remains the same. The `summary_name` now points to the new MAG's own ID, as its quality stats will be generated and stored under this new name.

## 5. Visualizing Progress

This design allows for the intuitive, hierarchical progress tracking you envisioned. A script could parse the manifest and generate a report like this, clearly showing the results of each refinement step:

```
                       Completeness    Contamination    Genome size
Asg1                       92.5            1.5              3.2M
├─ Asg1_flye               95.8            1.2              3.3M
│  └─ Asg1_flye_pilon      96.1            1.1              3.3M
└─ Asg1_hifiasm            94.2            2.1              3.1M
```

## 6. Advantages
- **Flexibility**: MAG names and file locations are fully independent.
- **Traceability**: The `parent_id` column provides a complete audit trail of how MAGs were improved.
- **Scalability**: The system is no longer limited by naming conventions and can manage thousands of MAGs across multiple refinement iterations.
- **Explicitness**: The manifest is a clear, human-readable, and machine-parsable declaration of all data dependencies.
