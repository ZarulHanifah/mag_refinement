from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional

from magrefine.mags import BaseMag


@dataclass(frozen=True)
class Metric:
    attribute: str  # attr name of BaseMag instance
    label: str  # text displayed in reports
    diff_key: str  # key inside chil.compare_to(parent)
    fmt: str = ".2f%"
    higher_is_better: bool = True

    def format_value(self, mag: BaseMag) -> str:
        """Extracts and formats the metric value from a MAG instance"""
        val = getattr(mag, self.attribute, None)
        if val is None:
            return f"{self.label}: N/A"
        suffix = "%" if self.fmt.endswith("%") else ""
        clean_fmt = self.fmt.rstrip("%")

        raw_str = f"{val:{clean_fmt}}{suffix}"
        return f"{self.label}: {raw_str}"

    def format_delta(self, diffs: Dict[str, float]) -> Optional[str]:
        """Calculates color and formatting for the change between parent and child"""
        val = diffs.get(self.diff_key)
        if val is None:
            return None

        if val == 0:
            color = "white"
        else:
            is_improvement = (val > 0) == self.higher_is_better
            color = "green" if is_improvement else "red"

        suffix = "%" if self.fmt.endswith("%") else ""
        clean_fmt = self.fmt.rstrip("%")

        # formatted_num = f"{val:+{clean_fmt}}{suffix}"
        formatted_val = f"{val:+{clean_fmt}}{suffix}"
        colored_val = f"[{color}]{formatted_val}[/{color}]"

        return f"{self.label.lower()}: {colored_val}"


DEFAULT_METRICS = [
    Metric(
        attribute="completeness",
        label="Comp",
        diff_key="completeness_diff",
        fmt=".2f%",
        higher_is_better=True,
    ),
    Metric(
        attribute="contamination",
        label="Contam",
        diff_key="contamination_diff",
        fmt=".2f%",
        higher_is_better=False,
    ),
    Metric(
        attribute="total_contigs",
        label="Contigs",
        diff_key="total_contigs_diff",
        fmt="d",
        higher_is_better=False,
    ),
]


class RefinementGraph:
    """
    Manages the lineage, derivation paths, and visualization reports
    for MAGs, fully decoupling relationships from the physical MAG properties.
    """

    def __init__(self, metrics: Optional[List[Metric]] = None):
        self.nodes: Dict[str, BaseMag] = {}
        # Adjacency list: parent_name -> list of child_names
        self.adjacency: Dict[str, List[str]] = defaultdict(list)
        # Parent lookup: child_name -> parent_name
        self.parents: Dict[str, str] = {}

        self.metrics = metrics if metrics is not None else DEFAULT_METRICS

    def add_mag(self, mag: BaseMag):
        """Registers a MAG node in the graph."""
        self.nodes[mag.name] = mag

    def add_derivation(self, parent: BaseMag, child: BaseMag):
        """Registers a parent -> child refinement step."""
        self.add_mag(parent)
        self.add_mag(child)

        if child.name not in self.adjacency[parent.name]:
            self.adjacency[parent.name].append(child.name)
        self.parents[child.name] = parent.name

    def get_parent(self, mag_name: str) -> Optional[BaseMag]:
        """Looks up the parent of a MAG name if it exists in the graph."""
        parent_name = self.parents.get(mag_name)
        return self.nodes.get(parent_name) if parent_name else None

    def get_children(self, mag_name: str) -> List[BaseMag]:
        """Gets all direct child nodes derived from a parent MAG name."""
        return [self.nodes[name] for name in self.adjacency[mag_name]]

    def improvement_report(self, child: BaseMag) -> Optional[Dict[str, float]]:
        """
        Returns quality metric differences between a child and its parent.
        """
        parent_name = self.parents.get(child.name)
        if not parent_name:
            return None
        parent = self.nodes[parent_name]
        return child.compare_to(parent)

    # --- SIMPLIFIED HELPERS ---

    def _get_metric_summary(self, mag: BaseMag) -> str:
        """Delegates formatting rules straight to the metric objects."""
        return " | ".join(metric.format_value(mag) for metric in self.metrics)

    def _get_delta_summary(self, child: BaseMag, parent: BaseMag) -> str:
        """Delegates delta calculations cleanly to the metric objects"""
        diffs = child.compare_to(parent) or {}
        delta_strings = []

        for metric in self.metrics:
            delta_str = metric.format_delta(diffs)
            if delta_str:
                delta_strings.append(delta_str)

        return f"↓ ({', '.join(delta_strings)})"

    # --- REPORTING METHODS ---

    def tree_report(self, root: BaseMag) -> str:
        """Generates a top-down tree showing this MAG and all descendant derivations."""
        lines = []
        straight_down = "│"

        def recurse(mag: BaseMag, prefix: str, is_last: bool):
            summary = self._get_metric_summary(mag)

            if not prefix:
                lines.append(f"└── {mag.name} ({summary})")
                child_prefix = "    "
            else:
                branch = "└──" if is_last else "├──"
                lines.append(f"{prefix}{branch} {mag.name} ({summary})")
                child_prefix = prefix + ("    " if is_last else "│   ")

            children_names = self.adjacency.get(mag.name, [])
            for idx, child_name in enumerate(children_names):
                child = self.nodes[child_name]
                child_is_last = idx == len(children_names) - 1

                # diffs = child.compare_to(mag)
                # comp_diff = diffs["completeness_diff"]
                # contam_diff = diffs["contamination_diff"]
                # contigs_diff = diffs["total_contigs_diff"]
                #
                # comp_color = (
                #     "green" if comp_diff > 0 else ("red" if comp_diff < 0 else "white")
                # )
                # contam_color = (
                #     "red"
                #     if contam_diff > 0
                #     else ("green" if contam_diff < 0 else "white")
                # )
                # contigs_color = (
                #     "red"
                #     if contigs_diff > 0
                #     else ("green" if contigs_diff < 0 else "white")
                # )
                #
                # comp_str = f"[{comp_color}]{comp_diff:+.2f}%[/{comp_color}]"
                # contam_str = f"[{contam_color}]{contam_diff:+.2f}%[/{contam_color}]"
                # contigs_str = f"[{contigs_color}]{contigs_diff:+d}[/{contigs_color}]"

                lines.append(f"{child_prefix}{straight_down}")
                lines.append(f"{child_prefix}{self._get_delta_summary(child, mag)}")
                # lines.append(
                #     f"{child_prefix}↓ (comp: {comp_str}, contam: {contam_str}, contigs: {contigs_str})"
                # )
                recurse(child, child_prefix, child_is_last)

        recurse(root, "", True)
        return "\n".join(lines)

    def lineage_report(self, leaf: BaseMag) -> str:
        """
        Generates a tree-like string showing the history of improvements
        from the original root MAG to this leaf MAG.
        """
        chain = []
        current = leaf
        while current is not None:
            chain.append(current)
            parent_name = self.parents.get(current.name)
            current = self.nodes.get(parent_name) if parent_name else None

        chain.reverse()

        lines = []
        for i, mag in enumerate(chain):
            indent = "    " * i

            # comp = mag.completeness
            # contam = mag.contamination
            # contigs = mag.total_contigs
            # mag_str = f"{indent}└── {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})"

            summary = self._get_metric_summary(mag)

            if i > 0:
                # prev_mag = chain[i - 1]
                # diffs = mag.compare_to(prev_mag)
                # comp_diff = diffs["completeness_diff"]
                # contam_diff = diffs["contamination_diff"]
                # contigs_diff = diffs["total_contigs_diff"]
                #
                # comp_color = (
                #     "green" if comp_diff > 0 else ("red" if comp_diff < 0 else "white")
                # )
                # contam_color = (
                #     "red"
                #     if contam_diff > 0
                #     else ("green" if contam_diff < 0 else "white")
                # )
                # contigs_color = (
                #     "red"
                #     if contigs_diff > 0
                #     else ("green" if contigs_diff < 0 else "white")
                # )
                #
                # comp_str = f"[{comp_color}]{comp_diff:+.2f}%[/{comp_color}]"
                # contam_str = f"[{contam_color}]{contam_diff:+.2f}%[/{contam_color}]"
                # contigs_str = f"[{contigs_color}]{contigs_diff:+d}[/{contigs_color}]"
                #
                # arrow_indent = "    " * i
                # delta_str = f"{arrow_indent}↓ (comp: {comp_str}, contam: {contam_str}, contigs: {contigs_str})"
                # lines.append(delta_str)
                lines.append(f"{indent}{self._get_delta_summary(mag, chain[i - 1])}")

            lines.append(f"{indent}└── {mag.name} ({summary})")

            # lines.append(mag_str)

        return "\n".join(lines)
