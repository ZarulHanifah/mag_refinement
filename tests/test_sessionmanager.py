import pytest

from magrefine.mags import Mag
from magrefine.sessionmanager import SessionManager


def test_session_manager_initialization(temp_project):
    """
    Tests that the SessionManager can be initialized correctly and
    that it holds the necessary file paths.
    """
    session = SessionManager(
        summary_path=temp_project["summary_path"],    # /tmp/summary.tsv
        mag_dir=temp_project["mag_dir"],              # /tmp/dereplicated_genomes/
        abund_dir=temp_project["abund_dir"],          # /tmp/depths/
    )

    # Check that the internal repositories/locators are created (optional but good practice)
    assert session.summary_repo is not None
    assert session.locator is not None

def test_session_manager_get_mag(mock_session_manager):
    """
    Tests the core functionality of the SessionManager: creating a valid Mag object.
    """

    mag_name = "C1A3_A_metabat.872"
    mag = mock_session_manager.get_mag(mag_name)

    # Verify the created object is a Mag instance
    assert isinstance(mag, Mag)

    # Verify the Mag object is correctly populated, same as in test_mag.py
    assert mag.name == mag_name
    assert len(mag.contigids) == 3
    assert mag.sample == "C1A3"
    assert mag.assembler == "myloasm"
    assert mag.total_contigs == 3
    assert mag.completeness == 63.39 # Check a value from the summary table

def test_session_manager_get_mag_raises_for_nonexistent_mag(mock_session_manager):
    """
    Tests that the session manager correctly raises an error when a MAG
    that doesn't exist in the summary is requested.
    """

    # This should raise a KeyError or a custom exception because the MAG
    # is not in the summary.tsv created by the fixture.
    with pytest.raises(KeyError):
        mock_session_manager.get_mag("nonexistent_mag_name")
