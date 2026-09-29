"""Every constant in VenomGap lives here.

Nothing in `model/`, `optimize/` or `validate/` may hard-code a numeric parameter. Fitted
parameters are written to and read from `data/derived/fitted_parameters.json`; everything that is
fixed a priori is in this module, with its justification beside it.

The severity weights `SIGMA` are frozen before any fitting. See `experiments/preregistration.md`.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Final

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------

PACKAGE_ROOT: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = PACKAGE_ROOT.parent.parent

DATA_DIR: Final[Path] = PROJECT_ROOT / "data"
RAW_DIR: Final[Path] = DATA_DIR / "raw"
CURATED_DIR: Final[Path] = DATA_DIR / "curated"
DERIVED_DIR: Final[Path] = DATA_DIR / "derived"
RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
FIGURES_DIR: Final[Path] = PROJECT_ROOT / "figures"

COMPOSITIONS_CSV: Final[Path] = CURATED_DIR / "compositions.csv"
SEVERITY_CSV: Final[Path] = CURATED_DIR / "severity_weights.csv"
CALIBRATION_CSV: Final[Path] = CURATED_DIR / "antivenomics_calibration.csv"
HOLDOUT_CSV: Final[Path] = CURATED_DIR / "holdout_targets.csv"
BURDEN_CSV: Final[Path] = CURATED_DIR / "state_snakebite_mortality.csv"
DISTRICTS_CSV: Final[Path] = CURATED_DIR / "districts.csv"
BITE_SHARE_CSV: Final[Path] = CURATED_DIR / "bite_attribution.csv"
RANGE_SPECIES_CSV_FALLBACK: Final[Path] = CURATED_DIR / "species_range.csv"

FITTED_PARAMS_JSON: Final[Path] = DERIVED_DIR / "fitted_parameters.json"
SEQUENCE_CACHE: Final[Path] = RAW_DIR / "uniprot"
OCCURRENCE_CACHE: Final[Path] = RAW_DIR / "gbif"

# Set VENOMGAP_OFFLINE=1 to forbid all network access. The demo and CI run with this set;
# every fetcher must then be satisfiable from cache or raise.
OFFLINE: Final[bool] = os.environ.get("VENOMGAP_OFFLINE", "0") == "1"

RANDOM_SEED: Final[int] = 20260929

# The four studies underlying the pre-registered holdout findings R1-R4. Their *compositions* are
# legitimate model inputs -- the model cannot predict anything about a population without knowing
# what its venom contains -- but none of their antivenom, immunorecognition or neutralisation
# results may enter a parameter fit. `model/calibrate.py` raises on any calibration row citing one.
HOLDOUT_STUDY_DOIS: Final[frozenset[str]] = frozenset(
    {
        "10.3390/toxins18010054",          # R1: E. c. sochureki, north-west India
        "10.3389/fphar.2021.768210",       # R2: Naja sagittifera, Andaman and Nicobar
        "10.1371/journal.pntd.0007899",    # R3: beyond the big four, North Indian B. caeruleus
        "10.1371/journal.pntd.0009659",    # R4: pathology-specific experimental antivenoms
    }
)

# --------------------------------------------------------------------------------------
# Toxin family vocabulary
# --------------------------------------------------------------------------------------

FAMILIES: Final[tuple[str, ...]] = (
    "3FTx",
    "PLA2",
    "SVMP",
    "SVSP",
    "LAAO",
    "CRISP",
    "CTL",
    "KSPI",
    "DIS",
    "NP",
    "VEGF",
    "NUC",
    "PDE",
    "AChE",
    "CVF",
    "other",
)

FAMILY_INDEX: Final[dict[str, int]] = {f: i for i, f in enumerate(FAMILIES)}

FAMILY_LONG_NAME: Final[dict[str, str]] = {
    "3FTx": "Three-finger toxin",
    "PLA2": "Phospholipase A2",
    "SVMP": "Snake venom metalloproteinase",
    "SVSP": "Snake venom serine proteinase",
    "LAAO": "L-amino acid oxidase",
    "CRISP": "Cysteine-rich secretory protein",
    "CTL": "C-type lectin / snaclec",
    "KSPI": "Kunitz-type serine proteinase inhibitor",
    "DIS": "Disintegrin",
    "NP": "Natriuretic peptide",
    "VEGF": "Snake venom VEGF",
    "NUC": "5'-nucleotidase / nucleotidase",
    "PDE": "Phosphodiesterase",
    "AChE": "Acetylcholinesterase",
    "CVF": "Cobra venom factor (complement C3 homologue)",
    "other": "Other / unassigned",
}

# --------------------------------------------------------------------------------------
# sigma_f — clinical severity weights. FIXED A PRIORI. NEVER FITTED. NEVER TUNED.
# --------------------------------------------------------------------------------------
# Each weight is the relative clinical consequence of leaving unit mass of that family
# unneutralised, on a 0-1 scale anchored at 1.0 = post-synaptic respiratory paralysis.
# The rationale column is reproduced verbatim into data/curated/severity_weights.csv and into
# the paper's parameter table. These values were written before any composition data was
# ingested and before any holdout target was inspected.
SIGMA: Final[dict[str, float]] = {
    "3FTx": 1.00,
    "PLA2": 0.95,
    "SVMP": 0.90,
    "SVSP": 0.75,
    "CTL": 0.55,
    "NP": 0.35,
    "DIS": 0.30,
    "LAAO": 0.25,
    "CRISP": 0.20,
    "VEGF": 0.20,
    "KSPI": 0.15,
    "CVF": 0.15,
    "NUC": 0.10,
    "PDE": 0.10,
    "AChE": 0.10,
    "other": 0.10,
}

SIGMA_RATIONALE: Final[dict[str, str]] = {
    "3FTx": (
        "Post-synaptic alpha-neurotoxins produce flaccid paralysis and respiratory arrest, the "
        "proximate cause of death in elapid envenoming. Anchor of the scale at 1.00."
    ),
    "PLA2": (
        "Pre-synaptic beta-neurotoxicity (krait), systemic myotoxicity and anticoagulant action. "
        "Irreversible pre-synaptic damage is not reversed by later antivenom, so near-anchor."
    ),
    "SVMP": (
        "Principal driver of viper pathology: systemic haemorrhage, venom-induced consumption "
        "coagulopathy, local necrosis and tissue loss. Leading cause of amputation and of "
        "haemorrhagic death."
    ),
    "SVSP": (
        "Thrombin-like and factor-activating enzymes causing defibrinogenation and incoagulable "
        "blood. Serious but more often correctable than SVMP-driven vessel-wall damage."
    ),
    "CTL": (
        "Snaclecs cause thrombocytopenia and platelet dysfunction, aggravating bleeding without "
        "being its primary cause."
    ),
    "NP": "Natriuretic and bradykinin-potentiating peptides contributing to early hypotension.",
    "DIS": "Disintegrins block platelet aggregation via integrins, contributing to bleeding.",
    "LAAO": "Local oedema, apoptosis and oxidative tissue injury; limited systemic consequence.",
    "CRISP": "Smooth-muscle and ion-channel effects; minor recorded clinical contribution.",
    "VEGF": "Increases vascular permeability, contributing to oedema and hypotension.",
    "KSPI": (
        "Mostly protease inhibitors; potentiating roles reported but little independent pathology "
        "in the Indian species of interest."
    ),
    "CVF": "Complement depletion; immunologically striking but of limited acute clinical weight.",
    "NUC": "Nucleotidases act as adjuvants to other toxins rather than causing primary pathology.",
    "PDE": "Purine-liberating enzyme; adjuvant role only.",
    "AChE": "Abundant in some elapids with no established independent clinical syndrome.",
    "other": "Unassigned or minor components; conservative low weight.",
}

# --------------------------------------------------------------------------------------
# kappa_f prior — antibody mass required per unit venom mass, before fitting.
# --------------------------------------------------------------------------------------
# Mechanistic prior from molar equivalence. Neutralisation needs antibody molecules, not antibody
# mass: to bind one mole of toxin you need of order one mole of F(ab')2. So the antibody *mass*
# needed per unit toxin *mass* scales as M_ab / M_toxin. A 7 kDa three-finger toxin therefore
# costs roughly seven times more antibody mass per milligram than a 50 kDa PIII metalloproteinase.
# This is why a cobra bite is antibody-expensive despite a smaller venom yield.
F_AB2_MASS_KDA: Final[float] = 100.0

FAMILY_MASS_KDA: Final[dict[str, float]] = {
    "3FTx": 7.0,
    "PLA2": 14.0,
    "SVMP": 47.0,
    "SVSP": 28.0,
    "LAAO": 55.0,
    "CRISP": 25.0,
    "CTL": 28.0,
    "KSPI": 7.0,
    "DIS": 7.5,
    "NP": 3.0,
    "VEGF": 27.0,
    "NUC": 60.0,
    "PDE": 100.0,
    "AChE": 60.0,
    "CVF": 150.0,
    "other": 30.0,
}

# Molar binding ratio: how many F(ab')2 are needed per toxin molecule for functional
# neutralisation. Held at 1.0 for all families a priori; the fitted per-family multiplier in
# calibration is what absorbs any real departure from 1:1.
KAPPA_PRIOR: Final[dict[str, float]] = {
    f: F_AB2_MASS_KDA / FAMILY_MASS_KDA[f] for f in FAMILIES
}

# Families with enough antivenomics support in the calibration split to justify a free
# per-family multiplier. Everything else keeps the prior exactly.
KAPPA_FITTED_FAMILIES: Final[tuple[str, ...]] = (
    "3FTx",
    "PLA2",
    "SVMP",
    "SVSP",
    "CTL",
    "LAAO",
    "CRISP",
    "KSPI",
)

# L2 penalties pulling the fitted values back toward their a priori settings.
#
# kappa_scale is a nuisance scale: it absorbs the venom-specific antibody fraction of Indian
# polyvalent antivenom, which is genuinely unknown per batch and commonly reported at only 10-20%
# of total protein. Pinning it near 1 would be pinning it to a number nobody has measured, so its
# prior is deliberately weak.
KAPPA_L2_PENALTY: Final[float] = 0.05

# theta_f is a physical recognition threshold with a real prior (conformational epitopes lose
# recognition below roughly 60-70% identity), so its prior is kept firm.
THETA_L2_PENALTY: Final[float] = 1.0

# --------------------------------------------------------------------------------------
# theta_f — sequence-identity recognition floor. Initialised a priori, fitted on calibration.
# --------------------------------------------------------------------------------------
# Conformational antibody epitopes span roughly 15-20 residues, so recognition degrades sharply
# below about 60-70% identity across the epitope. 0.65 is the prior for every family; the bounds
# keep the fit inside a physically meaningful range.
THETA_PRIOR: Final[float] = 0.65
THETA_BOUNDS: Final[tuple[float, float]] = (0.40, 0.85)
KAPPA_SCALE_BOUNDS: Final[tuple[float, float]] = (0.05, 20.0)
KAPPA_LOG_MULT_BOUNDS: Final[tuple[float, float]] = (-1.5, 1.5)

# A fitted theta_f whose perturbation by +/- this much changes the calibration objective by less
# than THETA_IDENTIFIABILITY_TOL is reported as *not identified* by the available data, and keeps
# its a priori value. Saying so is the point: an unidentified parameter that silently sits at its
# prior looks exactly like a tuned one unless the report distinguishes them.
THETA_PROBE_DELTA: Final[float] = 0.05
THETA_IDENTIFIABILITY_TOL: Final[float] = 1e-6

# --------------------------------------------------------------------------------------
# B and D — antibody budget and delivered venom dose. FIXED, swept in sensitivity.
# --------------------------------------------------------------------------------------
# Indian polyvalent antivenom is a lyophilised equine F(ab')2 preparation. A Bharat Serums vial
# contains 390 mg of lyophilised powder (95% CI 155-510 mg) of which 25.2% is protein, i.e. 98 mg
# per vial (95% CI 39-125 mg) -- doi:10.3390/toxins18010054, Discussion, citing its reference 29.
# The standard Indian national protocol initial dose is 10 vials. Only part of that protein is
# venom-specific; the specific fraction is absorbed into the fitted kappa_scale rather than
# guessed here, which is why B is a scale and not a claim about absolute potency.
VIAL_PROTEIN_MG: Final[float] = 98.0
DEFAULT_VIALS: Final[int] = 10
DEFAULT_B_MG: Final[float] = VIAL_PROTEIN_MG * DEFAULT_VIALS

# Delivered (injected) venom mass in a systemically envenomed bite, distinct from the much larger
# milked yield. Taken at 40 mg as a cross-species central estimate for a significant envenoming.
DEFAULT_D_MG: Final[float] = 40.0

# Sensitivity sweep grids.
B_SWEEP_MG: Final[tuple[float, ...]] = (200.0, 500.0, 1000.0, 2000.0, 4000.0)
D_SWEEP_MG: Final[tuple[float, ...]] = (10.0, 20.0, 40.0, 60.0, 100.0)

# --------------------------------------------------------------------------------------
# Spatial model
# --------------------------------------------------------------------------------------
# Candidate kernel length scales in km for the leave-one-population-out CV grid. The grid runs well
# past India's ~3000 km north-south extent deliberately: if the CV optimum sat at the top of a
# short grid it would be an artefact of the grid, and the only way to tell a real plateau from a
# truncation is to look beyond it.
ELL_GRID_KM: Final[tuple[float, ...]] = (
    50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 650.0, 800.0, 1000.0,
    1200.0, 1600.0, 2000.0, 3000.0, 5000.0, 10000.0,
)

# A length scale whose CV error is within this relative tolerance of the minimum is inside the
# plateau. Reporting the plateau alongside the minimiser is the honest way to present a
# flat-bottomed CV curve, and it is what stops a shallow minimum being quoted as a precise number.
ELL_PLATEAU_TOLERANCE: Final[float] = 0.01

# A district whose nearest same-species sampled population is further than this is reported as
# `unknown`, not as covered. Expressed as a multiple of the fitted length scale so that it
# follows the data rather than being a second free knob.
UNKNOWN_CUTOFF_ELL_MULTIPLE: Final[float] = 2.5

# Absolute backstop: no extrapolation beyond this distance regardless of fitted ell.
UNKNOWN_CUTOFF_MAX_KM: Final[float] = 1200.0

EARTH_RADIUS_KM: Final[float] = 6371.0088

# --------------------------------------------------------------------------------------
# Immunogen sets
# --------------------------------------------------------------------------------------
# The Big Four as actually used by Indian manufacturers: venom pooled from Naja naja,
# Bungarus caeruleus, Daboia russelii and Echis carinatus, sourced predominantly from the
# Irula Snake Catchers' Industrial Cooperative Society, Chengalpattu district, Tamil Nadu.
IRULA_LAT: Final[float] = 12.6819
IRULA_LON: Final[float] = 80.0000

# No complete family-level proteome is published for any of these species AT the locality the
# immunogen venom is collected from (preregistration amendment A1.3). The immunogen composition is
# therefore interpolated at the Irula coordinates by the same spatial model used for districts, so
# it carries the same fitted length scale, the same uncertainty and the same imputation flags.
# These are the four species whose venom is pooled, not four curated population rows.
BIG_FOUR_SPECIES: Final[tuple[str, str, str, str]] = (
    "Naja naja",
    "Bungarus caeruleus",
    "Daboia russelii",
    "Echis carinatus",
)

# Synthetic population ids for the interpolated immunogen, one per Big Four species.
BIG_FOUR_POP_IDS: Final[tuple[str, ...]] = tuple(
    f"{name.replace(' ', '_')}__IrulaTamilNadu" for name in BIG_FOUR_SPECIES
)

# Equal mixture weights, matching the standard description of the Indian polyvalent immunisation
# protocol. Non-equal weights are explored by the optimiser, not assumed here.
BIG_FOUR_WEIGHTS: Final[tuple[float, float, float, float]] = (0.25, 0.25, 0.25, 0.25)

# --------------------------------------------------------------------------------------
# Optimisation
# --------------------------------------------------------------------------------------
MAX_K: Final[int] = 6
EXACT_INSTANCE_POPULATIONS: Final[int] = 20
EXACT_INSTANCE_MAX_K: Final[int] = 4
LOCAL_SEARCH_MAX_ITERS: Final[int] = 200

# --------------------------------------------------------------------------------------
# Numerical tolerances
# --------------------------------------------------------------------------------------
SIMPLEX_TOL: Final[float] = 1e-6
EPS: Final[float] = 1e-12

# Families below this medically weighted share are ignored when naming the *worst* limiting
# family, so that a trace family cannot win the argument on rounding noise.
LIMITING_FAMILY_MIN_WEIGHTED_SHARE: Final[float] = 0.02

# --------------------------------------------------------------------------------------
# Disclaimer carried by every figure, API response and web screen
# --------------------------------------------------------------------------------------
DISCLAIMER: Final[str] = (
    "Research model. Not clinical guidance. Predictions are computational and require "
    "experimental validation."
)
