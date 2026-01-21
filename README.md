#  `execute_workflow.py` - Flexible Snakemake Workflow Orchestrator

This Python script provides a powerful and flexible way to orchestrate your Snakemake workflow, allowing you to dynamically select Metagenome-Assembled Genomes (MAGs) based on custom Python logic, load them in parallel, and execute the Snakemake pipeline in manageable chunks.

---

## 🚀 Features

*   **Python-driven MAG Selection:** Define complex filtering criteria using native Python code on `Mag` object attributes.
*   **Two-Step Filtering:** Utilizes an efficient pre-filter (pandas query) followed by a flexible post-filter (Python logic on `Mag` objects).
*   **Parallel MAG Loading:** Instantiates `Mag` objects in parallel using `ProcessPoolExecutor` for faster processing.
*   **Chunked Workflow Execution:** Run the Snakemake workflow on batches of MAGs, reducing resource overhead and improving fault tolerance.
<!-- *   **Snakemake API Integration:** Seamlessly calls Snakemake as a Python library. -->

---

## 📁 Project Structure

For `execute_workflow.py` to work correctly, your project should follow this structure:

```
the_workflow/
├── config/
│   ├── config.yaml
│   └── ...
├── rules/
│   ├── checkm2.smk
│   └── ...
├── src/
│   ├── models.py         # Contains Mag, SessionManager, etc.
│   └── usecase.py        # (Optional, if still in use)
├── execute_workflow.py   # This script
└── Snakefile             # Your main Snakemake workflow
```

---

## 🛠️ Setup

1.  **Dependencies:** Ensure you have the required Python packages installed:
    ```bash
    pip install snakemake pandas tqdm
    ```

2.  **`config.yaml` Configuration**: Before running the script, ensure that the paths in `config/config.yaml` are set correctly. The script relies on this file to find the necessary input data. Key paths to verify include:
    *   `summary_path`: Path to the MAGs summary TSV file.
    *   `mag_dir`: Directory containing the dereplicated MAG FASTA files.
    *   `abund_dir`: Directory containing the abundance/depth TSV files.

3.  **`Snakefile` Configuration:** Make sure your `Snakefile` is configured to accept the `mags_to_process` variable from the configuration. It should look something like this:

    ```python
    # Inside your Snakefile
    # ...
    if "mags_to_process" in config:
        mags = config["mags_to_process"]
    else:
        # Fallback to discovering MAGs from the filesystem if not provided
        # Or you can raise an error to enforce usage of execute_workflow.py
        print("INFO: No 'mags_to_process' config found. Discovering MAGs from filesystem.")
        mags = [re.sub(".fasta", "", i) for i in os.listdir("input_folder/dereplicated_genomes") if i.endswith(".fasta")]
    # ... rest of your Snakefile
    ```

4.  **`models.py`:** Ensure your `models.py` in `src/` contains the `SessionManager` with `get_mags_by_query` and the `Mag` class, as developed previously.

---

## 🚀 Usage

Navigate to your `the_workflow/` directory in the terminal and run the script. The script now supports various command-line arguments to control filtering and execution.

```bash
python execute_workflow.py --help
```

**Examples:**

*   **Run the workflow with default settings:**
    ```bash
    python execute_workflow.py
    ```

*   **List selected MAGs without running the workflow:**
    ```bash
    python execute_workflow.py --list-mags
    ```

*   **Filter with custom completeness and contamination:**
    ```bash
    python execute_workflow.py --completeness 80 --contamination 15
    ```

*   **Set the number of Snakemake jobs and chunk size:**
    ```bash
    python execute_workflow.py --jobs 50 --chunk-size 5
    ```

*   **Perform a Snakemake dry run:**
    ```bash
    python execute_workflow.py --dry-run
    ```

*   **Pass additional Snakemake arguments directly (e.g., to specify cores for local execution, or specific rules):**
    ```bash
    python execute_workflow.py --jobs 8 --cores 4  # --cores is a snakemake arg, not our script's
    python execute_workflow.py --dry-run -p # equivalent to --dry-run and printshellcmds
    python execute_workflow.py all_MAGs_output_file # run a specific target
    ```
    Note: Arguments after `--` are passed directly to Snakemake.

*   **Use a different configuration file:**
    ```bash
    python execute_workflow.py --configfile my_other_config.yaml
    ```

*   **Override data paths from the command line:**
    ```bash
    python execute_workflow.py --summary-path /data/new/summary.tsv --mag-dir /data/new/mags
    ```

---

## ⚙️ Adjusting Filtering Logic

The core flexibility of `execute_workflow.py` lies in the `apply_custom_filters` function. This is where you define **which MAGs you want to process**. While the initial `pre_filter_query` is excellent for broad-stroke filtering on summary statistics, the real power comes from performing nuanced filtering on the `Mag` objects themselves. This allows you to inspect properties, like the underlying contigs, that are not available in the summary file.

Open `execute_workflow.py` and locate the `apply_custom_filters` function. Here are some examples, from basic to advanced.

```python
# Inside execute_workflow.py

def apply_custom_filters(mag: Mag) -> bool:
    """
    Apply custom filtering logic to a single Mag object.
    <<< THIS IS THE SECTION YOU WILL MODIFY OFTEN FOR YOUR FILTERING CRITERIA. >>>
    Returns True if the MAG passes the filters, False otherwise.
    """
    # === Basic Filtering (based on summary data also available in the pre-filter) ===

    # --- Example 1: High-quality MAGs ---
    # Uses the built-in convenience method for Completeness >= 90 and Contamination <= 5.
    if mag.is_high_quality():
        return True

    # --- Example 2: Filter by assembler and a minimum GC content ---
    if mag.assembler == "myloasm" and mag.gc_content > 65:
        return True

    # === Advanced Filtering (inspecting contigs) ===
    # The following examples are only possible here, not in the pre-filter.

    # --- Example 3: Select complete, single-contig MAGs (e.g., plasmids, phages) ---
    # These are often high-value targets.
    if mag.total_contigs == 1 and mag.contigids[0].is_circular():
        return True

    # --- Example 4: Select HQ MAGs that have a large circular contig (>1 Mbp) ---
    # This could indicate a complete chromosome within a fragmented assembly.
    if mag.is_high_quality() and any(c.is_circular() and c.length > 1_000_000
                                     for c in mag.contigids):
        return True

    # --- Example 5: Find MAGs that are not HQ but contain a circular contig ---
    # Useful for finding interesting elements (e.g., viruses, plasmids) in lower-quality bins.
    if not mag.is_high_quality() and mag.contamination < 20 and any(c.is_circular()
                                                                    for c in mag.contigids):
        return True


    # Add as many conditions as you need. Use 'return True' for MAGs that pass.

    # If the MAG doesn't meet any of the 'return True' conditions above, it fails the filters.
    return False

# You can also adjust the pre-filter query in the select_mags function:
def select_mags(max_workers=4) -> list[str]:
    # ...
    pre_filter_query = "Completeness > 70 and Contamination < 20" # Adjust this for initial broad filtering
    # ...
```

**To adjust your filters:**
1. **FIRST PASS, for speed**: Modify the `pre_filter_query` string in the `select_mags` function for initial broad filtering using `pandas.query` syntax.
1.  **FINER FILTERING, including contig level**: Modify the conditions within `apply_custom_filters` using any attributes or methods available on the `Mag` object.

---

## 🔧 Configuring Snakemake Execution

The `execute_workflow.py` script now utilizes command-line arguments (parsed by `argparse`) to control Snakemake execution. This provides a more flexible way to configure your runs without modifying the script code directly.

**Key Command-Line Arguments:**

*   **`--configfile`**: Specifies the path to the base YAML configuration file. Defaults to `config/config.yaml`.
*   **`--summary-path`**: Overrides the `summary_path` from the config file.
*   **`--mag-dir`**: Overrides the `mag_dir` from the config file.
*   **`--abund-dir`**: Overrides the `abund_dir` from the config file.
*   **`--chunk-size`**: Controls how many MAGs are processed in each Snakemake call. Defaults to `10`.
    *   Example: `--chunk-size 50`
*   **`-j` or `--jobs`**: Sets the number of parallel jobs for Snakemake. This overrides the `njobs` setting from your `config.yaml` if specified.
    *   Example: `--jobs 100`
*   **`-n` or `--dry-run`**: Performs a Snakemake dry run, showing what would be executed without actually running anything. The script will still perform the unlock step before the dry run.
    *   Example: `--dry-run`
*   **Pass-through Snakemake Arguments**: Any arguments not recognized by `execute_workflow.py` will be passed directly to the underlying `snakemake` command. This allows you to use standard Snakemake flags (e.g., `--cores`, `--printshellcmds`, `--unlock` if you only want to unlock and not run the workflow, target rules, etc.).
    *   Example: `python execute_workflow.py --cores 8 --printshellcmds`

**Important Notes:**
*   The `--unlock` step is always performed before any Snakemake run to prevent issues with locked working directories.
*   The `njobs` value from `config.yaml` will be used as a fallback if `--jobs` is not specified on the command line.
*   The `chunk_size` is now exclusively controlled via the command-line argument `--chunk-size` (default `10`).

For a full list of available arguments, run `python execute_workflow.py --help`.

This setup provides you with a robust and highly customizable framework for managing your Snakemake-based metagenomics workflows.
