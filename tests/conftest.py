
from pathlib import Path

import pandas as pd
import pytest

from magrefine.abundancedb import AbundanceDB
from magrefine.mags import Mag
from magrefine.sessionmanager import (FilesystemLocator, SessionManager,
                                      SummaryRepository)

# from models import (AbundanceDB, ContigID, FilesystemLocator, Mag,
                    # SessionManager, SummaryRepository)


@pytest.fixture(scope="session")
def temp_project(tmp_path_factory):
    """
    Create a fake project directory ONCE per session with:
    - MAG FASTA file
    - summary.tsv
    - depth/abundance file
    """
    tmp_path = tmp_path_factory.mktemp("data")
    # sample_mag_name = "C1A3_A_metabat.872"
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

@pytest.fixture
def mock_checkm1_data():
    return "C1E5_M_metabat.1297	{'marker lineage': 'k__Archaea', '# genomes': 207, '# markers': 149, '# marker sets': 107, '0': 28, '1': 108, '2': 12, '3': 1, '4': 0, '5+': 0, 'Completeness': 78.43278217109992, 'Contamination': 11.214953271028037, 'GC': 0.4702310878198455, 'GC std': 0.014870769377286222, 'Genome size': 2347852, '# ambiguous bases': 0, '# scaffolds': 145, '# contigs': 145, 'Longest scaffold': 154524, 'Longest contig': 154524, 'N50 (scaffolds)': 31122, 'N50 (contigs)': 31122, 'Mean scaffold length': 16192.082758620689, 'Mean contig length': 16192.082758620689, 'Coding density': 0.8149998381499345, 'Translation table': 11, '# predicted genes': 2931, 'GCN0': ['PF01090', 'TIGR03683', 'PF01984', 'PF00867', 'PF01269', 'PF00736', 'TIGR00419', 'PF00749', 'TIGR00389', 'PF03950', 'PF01798', 'PF01192', 'PF01981', 'PF00164', 'TIGR00057', 'TIGR00522', 'PF01198', 'TIGR00344', 'PF01912', 'PF04560', 'PF04127', 'TIGR00422', 'PF04010', 'PF01982', 'TIGR00432', 'TIGR02076', 'PF00752', 'PF00832'], 'GCN1': ['TIGR02153', 'PF00181', 'TIGR00064', 'PF04563', 'TIGR00442', 'PF00237', 'TIGR00329', 'PF09249', 'PF00189', 'PF03947', 'PF00687', 'PF01201', 'PF08068', 'PF09377', 'PF00900', 'PF01864', 'PF00380', 'PF00562', 'PF04983', 'PF02005', 'TIGR03677', 'TIGR00468', 'PF00411', 'PF00935', 'TIGR03685', 'TIGR00289', 'PF00276', 'PF00177', 'PF05000', 'PF03874', 'PF03439', 'PF04919', 'PF03484', 'PF04561', 'PF13656', 'PF00623', 'TIGR00408', 'PF00861', 'TIGR01213', 'PF01280', 'PF01194', 'TIGR00134', 'PF00416', 'PF04567', 'TIGR01080', 'PF06418', 'PF00410', 'PF08071', 'PF03719', 'PF00466', 'PF00281', 'PF01282', 'PF01866', 'PF01287', 'PF05670', 'PF01868', 'PF01351', 'TIGR00336', 'PF04565', 'PF02978', 'PF00298', 'PF03876', 'PF00958', 'PF01193', 'PF00827', 'PF04566', 'PF01200', 'PF00573', 'PF01092', 'PF00750', 'PF00203', 'PF00831', 'PF01849', 'PF00347', 'TIGR03679', 'PF00572', 'TIGR01046', 'PF01780', 'PF11987', 'PF09173', 'PF04997', 'TIGR01018', 'PF01157', 'PF00673', 'PF01922', 'PF01246', 'PF01655', 'PF00252', 'PF00297', 'PF00398', 'PF00318', 'TIGR00425', 'PF00238', 'PF01191', 'PF01667', 'PF00366', 'PF04019', 'TIGR00392', 'PF00327', 'TIGR02338', 'PF01172', 'PF07541', 'TIGR00549', 'PF01000', 'PF00333', 'TIGR02389', 'PF03764', 'TIGR00270'], 'GCN2': ['PF00833', 'TIGR03665', 'PF01725', 'PF06026', 'PF00312', 'PF08069', 'PF13685', 'PF05221', 'PF00679', 'PF01015', 'TIGR03724', 'PF05746'], 'GCN3': ['TIGR00670'], 'GCN4': [], 'GCN5+': []}"

@pytest.fixture
def mock_checkm1_file(tmp_path, mock_checkm1_data):
    datafile = tmp_path / "mock_checkm1.tsv"
    datafile.write_text(mock_checkm1_data)
    return Path(datafile)
