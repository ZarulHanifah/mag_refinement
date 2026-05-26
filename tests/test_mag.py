from pathlib import Path

import pytest

from magrefine.contigids import GenericContigID
from magrefine.mags import RefinedMag
from magrefine.sessionmanager import SessionManager


def test_mag_initialization_builds_contigs(mock_session_manager: SessionManager):
    mag_name = "C1A3_A_metabat.872"
    mag = mock_session_manager.get_mag(mag_name)

    assert len(mag.contigids) == 3
    assert mag.sample == "C1A3"
    assert mag.assembler == "myloasm"
    assert mag.total_contigs == 3


def test_mag_average_coverage_total(mock_session_manager: SessionManager):
    mag_name = "C1A3_A_metabat.872"
    mag = mock_session_manager.get_mag(mag_name)

    cov = mag.average_coverage_total

    assert cov == pytest.approx(12.2464, rel=1e-2)


def test_mag_is_high_quality_true(mock_session_manager: SessionManager):
    mq_mag_name = "C1A3_A_metabat.872"
    hq_mag_name = "C1A3_M_semibin.175"

    mq_mag = mock_session_manager.get_mag(mq_mag_name)
    hq_mag = mock_session_manager.get_mag(hq_mag_name)

    assert mq_mag.is_high_quality() is False
    assert hq_mag.is_high_quality() is True


def test_refinedmag_initialization(mock_session_manager: SessionManager):
    """Test that a RefinedMag can be initialized correctly without external _data dictionaries."""
    mag_name = "C1A3_A_metabat.872"
    # Using existing MAG as a template
    original_mag = mock_session_manager.get_mag(mag_name)

    # Simulate a scaffolded new version
    refined_mag = RefinedMag(
        name="C1A3_A_metabat.872_scaffolded",
        _fp=Path("/fake/path/scaffolded.fasta"),
        contigids=original_mag.contigids[
            :2
        ],  # simulate merging contigs, so total_contigs is less
        completeness=96.0,
        contamination=1.0,
    )

    assert refined_mag.completeness == 96.0
    assert refined_mag.contamination == 1.0
    assert refined_mag.total_contigs == 2
    # assert refined_mag.parent is original_mag


def test_refinedmag_init_from_fasta_and_checkm2qualdir():
    mag_name = "C1E5_M_metabat.1297"

    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqualdir = Path("./tests/input_folder/checkm2/C1E5_M_metabat.1297")

    # mag = RefinedMag.from_checkm2qual(mag_name, _fp, checkmqualdir, parent=None)
    mag = RefinedMag.from_checkm2qual(mag_name, _fp, checkmqualdir)

    assert mag.completeness == 77.4
    assert mag.contamination == 7.54
    assert mag.total_contigs == 145
    # assert mag.parent is None


def test_refinedmag_init_from_fasta_and_checkm2qual():
    mag_name = "C1E5_M_metabat.1297"

    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqual = Path(
        "./tests/input_folder/checkm2/C1E5_M_metabat.1297/quality_report.tsv"
    )

    mag = RefinedMag.from_checkm2qual(mag_name, _fp, checkmqual)

    assert mag.completeness == 77.4
    assert mag.contamination == 7.54
    assert mag.total_contigs == 145
    # assert mag.parent is None


def test_refinedmag_init_checkm2_missing_file():
    mag_name = "C1E5_M_metabat.1297"
    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqual = Path("./non_existent_file.tsv")

    with pytest.raises(FileNotFoundError, match="CheckM2 quality path does not exist"):
        RefinedMag.from_checkm2qual(mag_name, _fp, checkmqual)


def test_refinedmag_init_checkm2_missing_in_dir(tmp_path):
    mag_name = "C1E5_M_metabat.1297"
    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    empty_dir = tmp_path / "empty_checkm2"
    empty_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="Could not find quality_report.tsv"):
        RefinedMag.from_checkm2qual(mag_name, _fp, empty_dir)


def test_refinedmag_init_from_fasta_and_checkm1qual(mock_checkm1_file):
    mag_name = "C1E5_M_metabat.1297"

    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")

    mag = RefinedMag.from_checkm1qual(mag_name, _fp, mock_checkm1_file)

    assert mag.completeness == pytest.approx(78.43, abs=0.1)
    assert mag.contamination == pytest.approx(11.21, abs=0.1)
    assert mag.total_contigs == pytest.approx(145)
    # assert mag.parent is None


def test_manual_mag_instantiation():
    """Demonstrates how to explicitly construct an original CheckM-style Mag without using fixtures."""
    from magrefine.mags import Mag

    mag_name = "MockSample_A_metabat.123"
    fake_fasta_path = Path("/path/to/mock.fasta")

    # These match the columns normally pulled from summary.tsv by SessionManager
    checkm_data = {
        "Completeness": 98.5,
        "Contamination": 1.2,
        "Total_Contigs": 4,
        "Contig_N50": 50000,
        "GC_Content": 45.5,
        "Max_Contig_Length": 100000,
        "genome_size": 250000,
        "average_coverage": 40.0,
        "tRNA counts": 42,
        "classification": "d__Bacteria;c__MockClass",
        "red_value": 0.99,
        "16S_rRNA": True,
        "23S_rRNA": False,
        "5S_rRNA": True,
    }

    contigs = [
        GenericContigID(">ctg1", {"sample1": 10.0}, _provided_length=100),
        GenericContigID(">ctg2", {"sample1": 20.0}, _provided_length=200),
    ]

    # Instantiation
    mag = Mag.from_summary_data(mag_name, checkm_data, fake_fasta_path, contigs)

    # Validating derived stats work when manually instantiated
    assert mag.sample == "MockSample"
    assert mag.assembler == "myloasm"
    assert mag.completeness == 98.5
    assert mag.total_contigs == 4
    assert mag.is_high_quality() is True
