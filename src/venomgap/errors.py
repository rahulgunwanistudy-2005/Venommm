"""Typed errors. Every failure path in VenomGap raises one of these, never a bare Exception."""

from __future__ import annotations


class VenomGapError(Exception):
    """Base class for every VenomGap failure."""


class DataValidationError(VenomGapError):
    """A curated data file violated its contract (bad schema, composition not summing to 1, ...)."""


class ProvenanceError(VenomGapError):
    """A curated row is missing its DOI, table reference or access date."""


class HoldoutLeakError(VenomGapError):
    """A fitting routine was handed a row labelled `holdout=True`.

    This is the guard that makes the blinded retrodiction study meaningful. It is raised, never
    warned, and never caught internally.
    """


class MissingCacheError(VenomGapError):
    """Offline mode was requested but a required network resource is not cached."""


class SequenceUnavailableError(VenomGapError):
    """No toxin sequences are available for a (species, family) pair, even after fallback."""


class CalibrationError(VenomGapError):
    """The parameter fit failed to converge or hit a degenerate solution."""


class OptimisationError(VenomGapError):
    """The siting optimiser was given an infeasible or ill-posed instance."""


class SpatialModelError(VenomGapError):
    """Spatial interpolation was asked for something it cannot answer."""
