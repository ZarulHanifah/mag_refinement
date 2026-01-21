
from pathlib import Path

import pandas as pd
import pytest

from models import (AbundanceDB, ContigID, FilesystemLocator, Mag,
                    SessionManager, SummaryRepository)


@pytest.fixture(scope="session")
def temp_project(tmp_path_factory):
    """
    Create a fake project directory ONCE per session with:
    - MAG FASTA file
    - summary.tsv
    - depth/abundance file
    """
    tmp_path = tmp_path_factory.mktemp("data")
    sample_mag_name = "C1A3_A_metabat.872"
    sample_mag_names = ["C1A3_A_metabat.872", "C1A3_M_semibin.175"]
    mag_dir = tmp_path / "dereplicated_genomes"
    abund_dir = tmp_path / "depths"
    mag_dir.mkdir()
    abund_dir.mkdir()

    # --- fake summary table ---
    real_summary_path = Path("./tests/input_folder/summary.tsv")
    df = pd.read_csv(real_summary_path, sep="\t", index_col=0)
    df = df.loc[sample_mag_names, :]

    summary_path = tmp_path / "summary.tsv"
    df.to_csv(summary_path, sep="\t")

    # --- fake MAG FASTA ---
    for n in sample_mag_names:
        real_fasta_path = Path(f"./tests/input_folder/dereplicated_genomes/{n}.fasta")
        fasta_path = mag_dir / f"{n}.fasta"
        with open(real_fasta_path) as source:
            with open(fasta_path, "w") as target:
                target.write(source.read())

    # --- fake abundance/depth files ---
    real_depth_file_path = Path("./tests/input_folder/mtbt_gen_depth/myloasm__C001_A3.tsv")
    with open(real_depth_file_path) as source:
        depth_content = source.read()

    # Create both required depth files for the tests, using the same content.
    myloasm_depth_path = abund_dir / "myloasm__C001_A3.tsv"
    medaka_depth_path = abund_dir / "medaka__C001_A3.tsv"
    myloasm_depth_path.write_text(depth_content)
    medaka_depth_path.write_text(depth_content)

    summary_repo = SummaryRepository(_summary_path=summary_path)
    locator = FilesystemLocator(mag_dir=mag_dir, abund_dir=abund_dir)

    return {
        "tmpdir": tmp_path,
        "summary_repo": summary_repo,
        "locator": locator,
        "summary_path": summary_path,
        "depth_file_path": myloasm_depth_path, # Keep returning one for other fixtures
        "mag_dir": mag_dir,
        "abund_dir": abund_dir,
    }


@pytest.fixture(scope="session")
def abundance_db(temp_project):
    """Provides an AbundanceDB instance based on the temp project's depth file."""
    return AbundanceDB(temp_project["depth_file_path"])


@pytest.fixture(scope="session")
def mock_session_manager(temp_project):
    return SessionManager(
        summary_path=temp_project["summary_path"],
        mag_dir=temp_project["mag_dir"],
        abund_dir=temp_project["abund_dir"],
    )


@pytest.fixture(scope="session")
def mock_mag(mock_session_manager: SessionManager):
    return mock_session_manager.get_mag("C1A3_A_metabat.872")


@pytest.fixture(scope="session")
def mock_contig(mock_mag: Mag):
    # This fixture provides a real ContigID from the mock Mag object.
    # To create a ContigID manually for testing, you would need a valid
    # header string and an abundance_db fixture instance, like so:
    #
    # def test_my_func(abundance_db):
    #   header = ">mycontig_len-100_circular-no_depth-1.0_duplicated-no"
    #   my_contig = ContigID(header, abundance_db)
    #   ...
    return mock_mag.contigids[0]

