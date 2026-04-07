import pytest

from magrefine.contigids import MyloContigID


def test_contigid_parses_header_fields(mock_contig):
    """Tests that properties are correctly parsed from the header."""
    contig = mock_contig
    assert contig.name == "u3558093ctg"
    assert contig.length == 256827
    assert contig.mylo_depth == "5-5-3"
    assert contig.is_circular is False
    assert contig.is_duplicated is False


def test_contigid_depth_sum(mock_contig):
    """
    Tests that the sum of depths is calculated correctly from abundance info.
    This also implicitly tests that abundance info is loaded correctly for the mock contig.
    """
    contig = mock_contig
    assert contig.depth_from_all_samples == pytest.approx(8.6927, rel=1e-2)


def test_contigid_with_no_abundance_data():
    """
    Tests the behavior of a ContigID created without abundance data,
    simulating a contig not found in the abundance file.
    """
    header = ">notpresent_len-10000_circular-yes_depth-10_duplicated-no"
    # A contig not found in the abundance file gets an empty dict.
    contig = MyloContigID.from_header(header, {}) 

    assert contig.get_abund_info() == {}
    assert contig.depth_from_all_samples == 0.0

from magrefine.contigids import GenericContigID


def test_manual_mylo_contig_instantiation():
    """Demonstrates how to manually instantiate a MyloContigID and what it extracts."""
    header = ">my_sample_contig1_len-2500_circular-yes_depth-15.5_duplicated-no"
    abundance_data = {"sample1": 10.0, "sample2": 5.5}

    # Initialize directly
    mylo_contig = MyloContigID.from_header(header, abundance_data)

    assert mylo_contig.name == "my_sample_contig1"
    assert mylo_contig.length == 2500
    assert mylo_contig.is_circular is True
    assert mylo_contig.mylo_depth == "15.5"
    assert mylo_contig.depth_from_all_samples == 15.5


def test_manual_generic_contig_instantiation():
    """Demonstrates how to manually instantiate a GenericContigID with explicit attribute assignment."""
    header = ">NODE_1_length_5000_cov_30.5" 
    abundance_data = {"sample1": 30.0}

    generic_contig = GenericContigID.from_header(header, abundance_data, provided_length=5000)

    assert generic_contig.name == "NODE_1_length_5000_cov_30.5"
    assert generic_contig.length == 5000
    assert generic_contig.depth_from_all_samples == 30.0
