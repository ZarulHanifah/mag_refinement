import pytest

from models import Mag, SessionManager


def test_session_manager_initialization(temp_project):
    """
    Tests that the SessionManager can be initialized correctly and
    that it holds the necessary file paths.
    """
    session = SessionManager(
        summary_path=temp_project["summary_path"],
        mag_dir=temp_project["mag_dir"],
        abund_dir=temp_project["abund_dir"],
    )

    # Check that the internal repositories/locators are created (optional but good practice)
    assert session.summary_repo is not None
    assert session.locator is not None

def test_session_manager_get_mag(temp_project):
    """
    Tests the core functionality of the SessionManager: creating a valid Mag object.
    """
    session = SessionManager(
        summary_path=temp_project["summary_path"],
        mag_dir=temp_project["mag_dir"],
        abund_dir=temp_project["abund_dir"],
    )

    mag_name = "C1A3_A_metabat.872"
    mag = session.get_mag(mag_name)

    # Verify the created object is a Mag instance
    assert isinstance(mag, Mag)

    # Verify the Mag object is correctly populated, same as in test_mag.py
    assert mag.name == mag_name
    assert len(mag.contigids) == 3
    assert mag.sample == "C1A3"
    assert mag.assembler == "myloasm"
    assert mag.total_contigs == 3
    assert mag.completeness == 63.39 # Check a value from the summary table

def test_session_manager_get_mag_raises_for_nonexistent_mag(temp_project):
    """
    Tests that the session manager correctly raises an error when a MAG
    that doesn't exist in the summary is requested.
    """
    session = SessionManager(
        summary_path=temp_project["summary_path"],
        mag_dir=temp_project["mag_dir"],
        abund_dir=temp_project["abund_dir"],
    )

    # This should raise a KeyError or a custom exception because the MAG
    # is not in the summary.tsv created by the fixture.
    with pytest.raises(KeyError):
        session.get_mag("nonexistent_mag_name")
