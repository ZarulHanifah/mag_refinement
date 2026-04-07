
import mmap
from dataclasses import field
from pathlib import Path


class AbundanceDB:
    """
    A memory-efficient, database-like interface for a single abundance TSV file,
    using memory mapping (mmap) to avoid loading the entire file into memory.
    This class is a context manager to ensure proper resource cleanup.
    """
    def __init__(self, depth_file_path: Path):
        self.depth_file_path = depth_file_path
        self._file = None
        self._mm: mmap.mmap = field(init=False)

        self._file = open(self.depth_file_path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

        # Read and process header once.
        self._mm.seek(0)
        header_line = self._mm.readline().decode().strip()
        self.header_fields = header_line.split('\t')

        original_var_cols = {h for h in self.header_fields if h.endswith("-var")}
        cleaned_var_names = {h.removesuffix("-var") for h in original_var_cols}

        self.original_keys_to_extract = []
        self.target_data_cols = []

        for hf in self.header_fields:
            if hf not in original_var_cols and hf in cleaned_var_names:
                self.original_keys_to_extract.append(hf)
                self.target_data_cols.append(hf.removesuffix("_merged.bam"))

        if not self.target_data_cols:
            raise ValueError(f"No valid sample columns found in header of '{self.depth_file_path}'.")

    def get_abund_for_contigs(self, contig_names_to_find: set[str]) -> dict[str, dict[str, float]]:
        """
        Finds and parses abundance data for a specific set of contigs in a single
        pass over the mmap'd file, using a substring search for each contig name.
        """
        results = {}
        contigs_left = contig_names_to_find.copy()

        self._mm.seek(0)
        _ = self._mm.readline()  # Skip header

        while contigs_left and (line_bytes := self._mm.readline()):
            found_name = None
            # Check if any of the required contig names exist as a substring in the current line.
            for name in contigs_left:
                if name.encode() in line_bytes:
                    found_name = name
                    break # Found a match, stop checking names for this line

            if found_name:
                try:
                    line_str = line_bytes.decode('utf-8').strip()
                    line_values = line_str.split('\t')
                    full_row_dict = dict(zip(self.header_fields, line_values))

                    result_dict = {
                        target_key: float(full_row_dict[original_key])
                        for target_key, original_key in zip(self.target_data_cols, self.original_keys_to_extract)
                        if original_key in full_row_dict
                    }

                    results[found_name] = result_dict
                    contigs_left.remove(found_name)
                except (UnicodeDecodeError, IndexError, ValueError):
                    # If line is malformed, skip it and continue.
                    continue

        return results

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self._mm and not self._mm.closed:
            self._mm.close()
        if self._file and not self._file.closed:
            self._file.close()

