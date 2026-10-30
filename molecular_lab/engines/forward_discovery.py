"""
Integrated forward molecular discovery workflow.

Version: 0.1
Rev.: 00
Status: First end-to-end composition-to-graph pipeline.

This workflow combines bounded composition enumeration, advisory diagnostics,
local reference lookup, restricted graph proposals, and optional structure
identification. It does not perform geometry optimization or stability tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

from molecular_lab.chemistry.plausibility import (
    PlausibilityAssessment,
    assess_composition,
)
from molecular_lab.core.models import (
    Composition,
    CompositionKind,
    MolecularGraph,
    SearchConstraints,
    SearchDomain,
)
from molecular_lab.database.reference_molecules import (
    MatchLevel,
    MoleculeLookupResult,
    identify_candidate,
)
from molecular_lab.engines.composition_enumerator import iter_compositions
from molecular_lab.engines.molecular_graphs import iter_saturated_acyclic_graphs
from molecular_lab.engines.structure_identification import (
    GraphIdentificationResult,
    identify_graph_against_catalogue,
)


@dataclass(frozen=True, slots=True)
class ForwardDiscoveryResult:
    """All currently available evidence for one enumerated composition."""

    composition: Composition
    diagnostics: PlausibilityAssessment
    reference_lookup: MoleculeLookupResult
    graphs: tuple[MolecularGraph, ...]
    graph_identifications: tuple[GraphIdentificationResult, ...]
    graph_generation_status: str
    graph_generation_note: str


def iter_forward_discovery(
    constraints: SearchConstraints,
    *,
    required_elements: tuple[str, ...] | set[str] | frozenset[str] = (),
    max_compositions: int = 1000,
    max_graphs_per_composition: int = 100,
    max_heavy_atoms: int = 8,
) -> Iterator[ForwardDiscoveryResult]:
    """
    Run the first integrated forward-discovery workflow.

    max_compositions bounds the number of compositions processed, not the
    number of raw combinations silently skipped. Every yielded result carries
    diagnostics and catalogue status, even if no graph can be proposed.

    The current graph prototype is molecular-only, neutral, saturated,
    acyclic, and limited to H/C/N/O/F/Cl. Out-of-scope compositions remain
    visible in the output with an explanatory graph-generation status.
    """
    for name, value in (
        ("max_compositions", max_compositions),
        ("max_graphs_per_composition", max_graphs_per_composition),
        ("max_heavy_atoms", max_heavy_atoms),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")

    if constraints.domain == SearchDomain.PERIODIC:
        raise ValueError(
            "Periodic-only discovery is not yet supported by the molecular forward pipeline"
        )

    processed = 0
    for composition in iter_compositions(
        constraints,
        required_elements=required_elements,
    ):
        if processed >= max_compositions:
            break
        processed += 1

        diagnostics = assess_composition(composition)
        if composition.kind == CompositionKind.MOLECULAR:
            reference_lookup = identify_candidate(composition)
        else:
            reference_lookup = MoleculeLookupResult(
                match_level=MatchLevel.NOT_APPLICABLE,
                matches=(),
                explanation=(
                    "The local reference catalogue contains molecular records only; "
                    "formula-unit compositions require a crystal/material reference database."
                ),
            )

        if composition.kind != CompositionKind.MOLECULAR:
            graphs: tuple[MolecularGraph, ...] = ()
            identifications: tuple[GraphIdentificationResult, ...] = ()
            status = "not_attempted_non_molecular"
            note = (
                "This composition is represented as a formula unit, not a "
                "molecule. Crystal/periodic structure generation is not implemented."
            )
        elif composition.net_charge != 0:
            graphs = ()
            identifications = ()
            status = "not_attempted_charged"
            note = (
                "The current graph prototype handles neutral species only. "
                "This charged composition has not been assessed structurally."
            )
        else:
            graphs = tuple(
                iter_saturated_acyclic_graphs(
                    composition,
                    max_candidates=max_graphs_per_composition,
                    max_heavy_atoms=max_heavy_atoms,
                )
            )
            identifications = tuple(
                identify_graph_against_catalogue(graph) for graph in graphs
            )
            if graphs:
                status = "graphs_proposed"
                note = (
                    "One or more connectivity hypotheses were generated within "
                    "the restricted saturated acyclic graph model."
                )
            else:
                status = "no_graph_within_prototype_scope"
                note = (
                    "No graph was generated by the current restricted model. "
                    "This is not evidence that the composition is impossible."
                )

        yield ForwardDiscoveryResult(
            composition=composition,
            diagnostics=diagnostics,
            reference_lookup=reference_lookup,
            graphs=graphs,
            graph_identifications=identifications,
            graph_generation_status=status,
            graph_generation_note=note,
        )


def run_forward_discovery(
    constraints: SearchConstraints,
    **kwargs: object,
) -> tuple[ForwardDiscoveryResult, ...]:
    """Convenience wrapper returning a bounded result tuple."""
    return tuple(iter_forward_discovery(constraints, **kwargs))
