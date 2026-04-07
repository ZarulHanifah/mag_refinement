# from models import FilesystemLocator, SummaryRepository


def test_summary_repository(temp_project):
    summary_repo = temp_project["summary_repo"]
    assert summary_repo.summary is not None
    assert "C1A3_A_metabat.872" in summary_repo.summary.index
    mag_data = summary_repo.get_mag_data("C1A3_A_metabat.872")
    assert mag_data["Completeness"] == 63.39

def test_filesystem_locator(temp_project):
    locator = temp_project["locator"]
    mag_name = "C1A3_A_metabat.872"
    mag_file = locator.find_mag_file(mag_name)
    assert mag_file.exists()
    assert mag_file.name == f"{mag_name}.fasta"
