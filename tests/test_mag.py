
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
        contigids=original_mag.contigids[:2],  # simulate merging contigs, so total_contigs is less
        completeness=96.0,
        contamination=1.0,
        parent=original_mag
    )

    assert refined_mag.completeness == 96.0
    assert refined_mag.contamination == 1.0
    assert refined_mag.total_contigs == 2
    assert refined_mag.parent is original_mag

def test_refinedmag_init_from_fasta_and_checkm2qualdir():
    mag_name = "C1E5_M_metabat.1297"

    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqualdir = Path("./tests/input_folder/checkm2/C1E5_M_metabat.1297")

    mag = RefinedMag.from_checkm2qual(mag_name, _fp, checkmqualdir, parent=None)

    assert mag.completeness == 77.4
    assert mag.contamination == 7.54
    assert mag.total_contigs == 145
    assert mag.parent is None 

def test_refinedmag_init_from_fasta_and_checkm2qual():
    mag_name = "C1E5_M_metabat.1297"

    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqual = Path("./tests/input_folder/checkm2/C1E5_M_metabat.1297/quality_report.tsv")

    mag = RefinedMag.from_checkm2qual(mag_name, _fp, checkmqual, parent=None)

    assert mag.completeness == 77.4
    assert mag.contamination == 7.54
    assert mag.total_contigs == 145
    assert mag.parent is None 

def test_refinedmag_init_checkm2_missing_file():
    mag_name = "C1E5_M_metabat.1297"
    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqual = Path("./non_existent_file.tsv")

    with pytest.raises(FileNotFoundError, match="CheckM2 quality path does not exist"):
        RefinedMag.from_checkm2qual(mag_name, _fp, checkmqual, parent=None)

def test_refinedmag_init_checkm2_missing_in_dir(tmp_path):
    mag_name = "C1E5_M_metabat.1297"
    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    empty_dir = tmp_path / "empty_checkm2"
    empty_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="Could not find quality_report.tsv"):
        RefinedMag.from_checkm2qual(mag_name, _fp, empty_dir, parent=None)

def test_refinedmag_init_from_fasta_and_checkm1qual(mock_checkm1_file):
    mag_name = "C1E5_M_metabat.1297"

    _fp = Path("./tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")

    mag = RefinedMag.from_checkm1qual(mag_name, _fp, mock_checkm1_file, parent=None)

    assert mag.completeness == pytest.approx(78.43, abs=0.1)
    assert mag.contamination == pytest.approx(11.21, abs=0.1)
    assert mag.total_contigs == pytest.approx(145)
    assert mag.parent is None 


def test_refinedmag_improvement_report():
    """Test the comparison logic and metric diffing between a RefinedMag and its original parent."""
    # Create a synthetic parent MAG
    parent_mag = RefinedMag(
        name="Parent",
        _fp=Path("parent.fa"),
        contigids=[ GenericContigID(str(i)) for i in range(4) ], # Dummy values just for calculating len(contigids) = 4
        completeness=90.0,
        contamination=5.0
    )

    # Create a child MAG simulating an improvement
    child_mag = RefinedMag(
        name="Child",
        _fp=Path("child.fa"),
        contigids=[ GenericContigID(str(i)) for i in range(2) ], # Decreased count from 4 to 2 (better contiguity)
        completeness=95.0, # Improved by +5.0%
        contamination=1.0, # Improved by -4.0%
        parent=parent_mag
    )

    report = child_mag.improvement_report()

    assert report is not None
    assert report["completeness_diff"] == 5.0
    assert report["contamination_diff"] == -4.0
    assert report["total_contigs_diff"] == -2


def test_refinedmag_lineage_report():
    parent_mag = RefinedMag(
        name="Root_MAG",
        _fp=Path("parent.fa"),
        contigids=[ GenericContigID(str(i)) for i in range(4) ], 
        completeness=90.0,
        contamination=5.0
    )

    child_mag = RefinedMag(
        name="Scaffolded_MAG",
        _fp=Path("child.fa"),
        contigids=[ GenericContigID(str(i)) for i in range(2) ], 
        completeness=95.0, 
        contamination=1.0, 
        parent=parent_mag
    )

    report_string = child_mag.lineage_report()

    assert "└── Root_MAG (Comp: 90.00% | Contam: 5.00% | Contigs: 4)" in report_string
    assert "    ↓ (comp: [green]+5.00%[/green], contam: [green]-4.00%[/green], contigs: [green]-2[/green])" in report_string
    assert "    └── Scaffolded_MAG (Comp: 95.00% | Contam: 1.00% | Contigs: 2)" in report_string

def test_refinedmag_tree_report():
    root_mag = RefinedMag(
        name="Root",
        _fp=Path(""),
        contigids=[GenericContigID('1')]*40,
        completeness=90.0,
        contamination=5.0
    )

    branch_a = RefinedMag(
        name="Branch_A",
        _fp=Path(""),
        contigids=[GenericContigID('1')]*30,
        completeness=94.5,
        contamination=4.0,
        parent=root_mag,
    )

    branch_b = RefinedMag(
        name="Branch_B",
        _fp=Path(""),
        contigids=[GenericContigID('1')]*32,
        completeness=95.0,
        contamination=4.5,
        parent=root_mag,
    )

    final_a_node = RefinedMag(
        name="Final_A",
        _fp=Path(""),
        contigids=[GenericContigID('1')]*28,
        completeness=95.5,
        contamination=3.5,
        parent=branch_a,
    )

    report_string = root_mag.tree_report()

    assert "└── Root (Comp: 90.00% | Contam: 5.00% | Contigs: 40)" in report_string
    assert "    ├── Branch_A (Comp: 94.50% | Contam: 4.00% | Contigs: 30)" in report_string
    assert "    │   └── Final_A (Comp: 95.50% | Contam: 3.50% | Contigs: 28)" in report_string
    assert "    └── Branch_B (Comp: 95.00% | Contam: 4.50% | Contigs: 32)" in report_string

def test_manual_mag_instantiation():
    """Demonstrates how to explicitly construct an original CheckM-style Mag without using fixtures."""
    from magrefine.mags import Mag

    mag_name = "MockSample_A_metabat.123"
    fake_fasta_path = Path("/path/to/mock.fasta")

    # These match the columns normally pulled from summary.tsv by SessionManager
    checkm_data = {
        'Completeness': 98.5,
        'Contamination': 1.2,
        'Total_Contigs': 4,
        'Contig_N50': 50000,
        'GC_Content': 45.5,
        'Max_Contig_Length': 100000,
        'genome_size': 250000,
        'average_coverage': 40.0,
        'tRNA counts': 42,
        'classification': "d__Bacteria;c__MockClass",
        'red_value': 0.99,
        '16S_rRNA': True,
        '23S_rRNA': False,
        '5S_rRNA': True
    }

    contigs = [
        GenericContigID(">ctg1", {"sample1": 10.0}, _provided_length=100),
        GenericContigID(">ctg2", {"sample1": 20.0}, _provided_length=200)
    ]

    # Instantiation
    mag = Mag.from_summary_data(mag_name, checkm_data, fake_fasta_path, contigs)

    # Validating derived stats work when manually instantiated
    assert mag.sample == "MockSample"
    assert mag.assembler == "myloasm"
    assert mag.completeness == 98.5
    assert mag.total_contigs == 4
    assert mag.is_high_quality() is True
