from collections import defaultdict
from typing import Dict, List, Optional
from magrefine.mags import BaseMag

class RefinementGraph:
    """
    Manages the lineage, derivation paths, and visualization reports 
    for MAGs, fully decoupling relationships from the physical MAG properties.
    """
    def __init__(self):
        self.nodes: Dict[str, BaseMag] = {}
        # Adjacency list: parent_name -> list of child_names
        self.adjacency: Dict[str, List[str]] = defaultdict(list)
        # Parent lookup: child_name -> parent_name
        self.parents: Dict[str, str] = {}

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

    def tree_report(self, root: BaseMag) -> str:
        """Generates a top-down tree showing this MAG and all descendant derivations."""
        lines = []
        straight_down = "│"

        def recurse(mag: BaseMag, prefix: str, is_last: bool):
            comp, contam, contigs = mag.completeness, mag.contamination, mag.total_contigs

            if not prefix:
                lines.append(f"└── {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})")
                child_prefix = "    "
            else:
                branch = "└──" if is_last else "├──"
                lines.append(f"{prefix}{branch} {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})")
                child_prefix = prefix + ("    " if is_last else "│   ")

            children_names = self.adjacency.get(mag.name, [])
            for idx, child_name in enumerate(children_names):
                child = self.nodes[child_name]
                child_is_last = (idx == len(children_names) - 1)

                diffs = child.compare_to(mag)
                comp_diff = diffs["completeness_diff"]
                contam_diff = diffs["contamination_diff"]
                contigs_diff = diffs["total_contigs_diff"]

                comp_color = "green" if comp_diff > 0 else ("red" if comp_diff < 0 else "white")
                contam_color = "red" if contam_diff > 0 else ("green" if contam_diff < 0 else "white")
                contigs_color = "red" if contigs_diff > 0 else ("green" if contigs_diff < 0 else "white")

                comp_str = f"[{comp_color}]{comp_diff:+.2f}%[/{comp_color}]"
                contam_str = f"[{contam_color}]{contam_diff:+.2f}%[/{contam_color}]"
                contigs_str = f"[{contigs_color}]{contigs_diff:+d}[/{contigs_color}]"

                lines.append(f"{child_prefix}{straight_down}")
                lines.append(f"{child_prefix}↓ (comp: {comp_str}, contam: {contam_str}, contigs: {contigs_str})")
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

            comp = mag.completeness
            contam = mag.contamination
            contigs = mag.total_contigs
            mag_str = f"{indent}└── {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})"

            if i > 0:
                prev_mag = chain[i - 1]
                diffs = mag.compare_to(prev_mag)
                comp_diff = diffs["completeness_diff"]
                contam_diff = diffs["contamination_diff"]
                contigs_diff = diffs["total_contigs_diff"]

                comp_color = "green" if comp_diff > 0 else ("red" if comp_diff < 0 else "white")
                contam_color = "red" if contam_diff > 0 else ("green" if contam_diff < 0 else "white")
                contigs_color = "red" if contigs_diff > 0 else ("green" if contigs_diff < 0 else "white")

                comp_str = f"[{comp_color}]{comp_diff:+.2f}%[/{comp_color}]"
                contam_str = f"[{contam_color}]{contam_diff:+.2f}%[/{contam_color}]"
                contigs_str = f"[{contigs_color}]{contigs_diff:+d}[/{contigs_color}]"

                arrow_indent = "    " * i
                delta_str = f"{arrow_indent}↓ (comp: {comp_str}, contam: {contam_str}, contigs: {contigs_str})"
                lines.append(delta_str)

            lines.append(mag_str)

        return "\n".join(lines)
