
import pytest

from models import SessionManager


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
