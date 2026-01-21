import pytest

from models import ContigID


def test_contigid_parses_header_fields(mock_contig):
    """Tests that properties are correctly parsed from the header."""
    contig = mock_contig
    assert contig.name == "u3558093ctg"
    assert contig.length == 256827
    assert contig.mylo_depth == "5-5-3"
    assert contig.is_circular() is False
    assert contig.is_duplicated() is False


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
    contig = ContigID(header, {}) 

    assert contig.get_abund_info() == {}
    assert contig.depth_from_all_samples == 0.0


