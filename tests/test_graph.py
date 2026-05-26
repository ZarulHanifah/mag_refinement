from pathlib import Path

from magrefine.contigids import GenericContigID
from magrefine.graph import RefinementGraph
from magrefine.mags import RefinedMag


def test_refinedmag_improvement_report():
    """Test the comparison logic and metric diffing between a RefinedMag and its original parent."""
    # Create graph

    graph = RefinementGraph()

    # Create a synthetic parent MAG
    parent_mag = RefinedMag(
        name="Parent",
        _fp=Path("parent.fa"),
        contigids=[
            GenericContigID(str(i)) for i in range(4)
        ],  # Dummy values just for calculating len(contigids) = 4
        completeness=90.0,
        contamination=5.0,
    )
    graph.add_mag(parent_mag)

    # Create a child MAG simulating an improvement
    child_mag = RefinedMag(
        name="Child",
        _fp=Path("child.fa"),
        contigids=[
            GenericContigID(str(i)) for i in range(2)
        ],  # Decreased count from 4 to 2 (better contiguity)
        completeness=95.0,  # Improved by +5.0%
        contamination=1.0,  # Improved by -4.0%
    )
    graph.add_derivation(parent=parent_mag, child=child_mag)
    report = graph.improvement_report(child=child_mag)

    # report = child_mag.compare_to(parent_mag)

    # report = graph.improvement_report()

    assert report is not None
    assert report["completeness_diff"] == 5.0
    assert report["contamination_diff"] == -4.0
    assert report["total_contigs_diff"] == -2


def test_refinedmag_lineage_report():
    graph = RefinementGraph()
    parent_mag = RefinedMag(
        name="Root_MAG",
        _fp=Path("parent.fa"),
        contigids=[GenericContigID(str(i)) for i in range(4)],
        completeness=90.0,
        contamination=5.0,
    )
    graph.add_mag(parent_mag)

    child_mag = RefinedMag(
        name="Scaffolded_MAG",
        _fp=Path("child.fa"),
        contigids=[GenericContigID(str(i)) for i in range(2)],
        completeness=95.0,
        contamination=1.0,
        # parent=parent_mag,
    )
    graph.add_derivation(parent=parent_mag, child=child_mag)
    report_string = graph.lineage_report(child_mag)

    # report_string = child_mag.lineage_report()

    assert "└── Root_MAG (Comp: 90.00% | Contam: 5.00% | Contigs: 4)" in report_string
    assert (
        "    ↓ (comp: [green]+5.00%[/green], contam: [green]-4.00%[/green], contigs: [green]-2[/green])"
        in report_string
    )
    assert (
        "    └── Scaffolded_MAG (Comp: 95.00% | Contam: 1.00% | Contigs: 2)"
        in report_string
    )


def test_refinedmag_tree_report():
    graph = RefinementGraph()
    root_mag = RefinedMag(
        name="Root",
        _fp=Path(""),
        contigids=[GenericContigID("1")] * 40,
        completeness=90.0,
        contamination=5.0,
    )
    graph.add_mag(root_mag)

    branch_a = RefinedMag(
        name="Branch_A",
        _fp=Path(""),
        contigids=[GenericContigID("1")] * 30,
        completeness=94.5,
        contamination=4.0,
    )
    graph.add_derivation(root_mag, branch_a)

    branch_b = RefinedMag(
        name="Branch_B",
        _fp=Path(""),
        contigids=[GenericContigID("1")] * 32,
        completeness=95.0,
        contamination=4.5,
    )
    graph.add_derivation(root_mag, branch_b)

    final_a_node = RefinedMag(
        name="Final_A",
        _fp=Path(""),
        contigids=[GenericContigID("1")] * 28,
        completeness=95.5,
        contamination=3.5,
    )
    graph.add_derivation(branch_a, final_a_node)

    report_string = graph.tree_report(root_mag)

    assert "└── Root (Comp: 90.00% | Contam: 5.00% | Contigs: 40)" in report_string
    assert (
        "    ├── Branch_A (Comp: 94.50% | Contam: 4.00% | Contigs: 30)" in report_string
    )
    assert (
        "    │   └── Final_A (Comp: 95.50% | Contam: 3.50% | Contigs: 28)"
        in report_string
    )
    assert (
        "    └── Branch_B (Comp: 95.00% | Contam: 4.50% | Contigs: 32)" in report_string
    )
