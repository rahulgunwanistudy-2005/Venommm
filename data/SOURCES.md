# Data sources

Every source used, with its licence, the date it was accessed, and the exact query, table or DOI the
numbers came from. Raw downloads are gitignored and reproducible by the fetch scripts; the derived
tables in `data/curated/` are committed.

Accessed: **2026-09-29** unless stated otherwise.

---

## 1. Venom proteomes — `data/curated/compositions.csv`

29 population-level proteomes, 9 species, 14 states or regions, 6 studies. 22 are Indian.
Every row carries a DOI or PMC id, a table reference, an access date and a proteomic method.

| Study | DOI | Licence | What was taken | Rows |
|---|---|---|---|---|
| Tasoulis & Isbister 2017, *Toxins* 9(9):290 | `10.3390/toxins9090290` | CC BY 4.0 | Tables 1 and 2, family-level relative abundances for elapids and Viperinae | 9 |
| Senji Laxme et al. 2021, *PLoS NTD*, pan-India *Naja naja* | `10.1371/journal.pntd.0009150` | CC BY 4.0 | S2A–S2C Tables, full proteomes for Punjab, West Bengal and Rajasthan | 3 |
| Senji Laxme et al. 2021, *PLoS NTD*, pan-India *Daboia russelii* | `10.1371/journal.pntd.0009247` | CC BY 4.0 | S3 and S4 Tables (Western Ghats, Deccan plateau) plus the four abundances stated numerically in the Results for the Gangetic plain | 3 |
| Senji Laxme et al. 2019, *PLoS NTD*, beyond the big four | `10.1371/journal.pntd.0007899` | CC BY 4.0 | S2A–S2H Tables, eight full proteomes | 8 |
| Bhatia & Sunagar 2021, *Front. Pharmacol.*, Andaman cobra | `10.3389/fphar.2021.768210` | CC BY 4.0 | Figure 3B and the Results text, 3FTx and PLA2 for *N. sagittifera* and mainland *N. naja* | 2 |
| Kumar et al. 2026, *Toxins* 18(1):54, *E. c. sochureki* north-west India | `10.3390/toxins18010054` | CC BY 4.0 | Results §2.3 and Figure 3, three site-level proteomes (Sam, Barmer, Pokhran) | 3 |

### How the numbers were handled

- **Family mapping.** Source labels are mapped onto the 16-family vocabulary in
  `ingest/compositions.py`. Neurotoxic and cytotoxic 3FTx are summed. Snaclec and C-type lectin are
  summed into `CTL`. SVMP subclasses PI, PII and PIII are summed into `SVMP`. Vespryn, cystatin,
  NGF, hyaluronidase, phospholipase B, serpin, calreticulin, aminopeptidase and SVMP-inhibitor have
  no vocabulary slot and fold into `other`.
- **Unreported residual.** Where a study characterised less than all of whole venom, the shortfall
  is assigned to `other`. It is **not** redistributed across the measured families, because that
  would invent abundance for families the study did not see. Four rows carry a material residual
  and are automatically flagged `partial_table_renormalised`:
  `Naja_naja__Karnataka` (33%), `Naja_sagittifera__Andaman` (28%), `Naja_naja__EasternIndia` (18%),
  `Daboia_russelii__WestBengal` (13%).
- **Rejection rule.** A row characterising under 55% of whole venom is rejected outright.
- **Coordinates.** Where a study reports a state but not a point, the state or named-town centroid
  is used and the row's `note` says so. `Echis_carinatus_sochureki__{Sam,Barmer,Pokhran}` use the
  named town or district centroid; the article reports no lat/lon.

### Exclusions, stated rather than silent

- **Tamil Nadu *Echis carinatus*.** The only published composition for the immunogen-source
  locality gives two families (SVMP 23.06%, snaclec 32.98%; `10.3390/toxins18010054`, Discussion,
  citing its reference 2), covering 56% of whole venom. This is outside the 45% residual limit, so
  it is excluded rather than used. It is the reason pre-registration amendment A1.3 interpolates
  the immunogen composition rather than curating a proxy row.
- **Tasoulis *E. c. sochureki* row.** Excluded because the compilation does not state a sampling
  locality, and inventing coordinates for it would corrupt the spatial model.
- **Gangetic plain *D. russelii* S2 Table.** Not released with the article. The four abundances
  stated numerically in the Results text are used and the row is flagged.

### The known disagreement

Two studies report Rajasthan *E. c. sochureki* and disagree sharply: SVMP at 55–67% in the 2026
proteome against 3% in the 2019 one, with PLA2 correspondingly reversed. The 2026 paper flags the
discrepancy against the 2019 study itself. **Both are curated**, a test asserts the conflict still
exists, and pre-registration amendment A1.1 resolves how R1 is evaluated across them.

---

## 2. Toxin sequences — UniProt / ToxProt

REST API, no key required. Cached to `data/raw/uniprot/`, summarised in
`data/derived/sequence_cache_summary.json`.

```
https://rest.uniprot.org/uniprotkb/search
  ?query=organism_name:"<species>" AND keyword:KW-0800
  &format=tsv
  &fields=accession,id,protein_name,protein_families,organism_name,organism_id,length,sequence,reviewed,keywordid
  &size=500          # cursor-paginated via the Link header
```

`KW-0800` is UniProt's **Toxin** keyword, the annotation the ToxProt programme applies. Without it
an organism query returns the whole proteome. Species with fewer than 8 keyword-tagged
classifiable entries fall back to a whole-organism query; the fallback is recorded as
`broadened_query` in the cache file so the change of method stays visible. It applied to
*Naja sagittifera* and *Bungarus caeruleus*.

Licence: UniProt is CC BY 4.0.

| Species | Sequences kept | Families |
|---|---|---|
| *Naja naja* | 144 | 9 |
| *Daboia russelii* | 66 | 10 |
| *Echis carinatus* | 47 | 5 |
| *Naja kaouthia* | 30 | 6 |
| *Bungarus fasciatus* | 24 | 6 |
| *Echis carinatus sochureki* | 15 | 5 |
| *Bungarus caeruleus* | 7 | 2 |
| *Naja sagittifera* | 7 | 2 |
| *Bungarus sindanus* | 0 | — falls back to congeners |

Cross-reactivity is computed as mean pairwise identity from global Needleman–Wunsch alignment,
BLOSUM62, gap open −11, gap extend −1, identity over aligned columns. 82 between-species values, of
which 41 rest on two sequences or fewer and are flagged thin. Thin estimates are kept rather than
discarded: *N. sagittifera*'s 3FTx is a single entry and accounts for 70% of its venom.

---

## 3. Species occurrence — GBIF

Free API, no key. Cached to `data/raw/gbif/`.

```
https://api.gbif.org/v1/occurrence/search
  ?scientificName=<species>&country=IN&hasCoordinate=true&hasGeospatialIssue=false&limit=300
```

**Occurrence counts are never used as abundance.** They are presence-only and effort-biased. They
are used in one direction only — to decide whether a species is plausibly present in a district —
and the count of nearby records is discarded so it cannot leak into a weighting. A test asserts
this.

| Species | Georeferenced Indian records |
|---|---|
| *Naja naja* | 1749 |
| *Daboia russelii* | 1020 |
| *Bungarus caeruleus* | 657 |
| *Echis carinatus* | 651 |
| *Bungarus fasciatus* | 116 |
| *Naja kaouthia* | 109 |
| *Echis carinatus sochureki* | 43 |
| *Bungarus sindanus* | 17 |
| *Naja sagittifera* | 5 |

Species below 25 records fall back to a published range description in
`data/curated/species_range.csv`, which is why the Andaman endemic *N. sagittifera* is not shrunk to
five points. GBIF data are CC BY / CC0 depending on publisher; the aggregated occurrence search is
used under GBIF's terms of use.

---

## 4. Snakebite burden — `data/curated/state_snakebite_mortality.csv`

Suraweera W, Warrell D, Whitaker R, et al. *Trends in snakebite deaths in India from 2000 to 2019 in
a nationally representative mortality study.* eLife 2020;9:e54076. `10.7554/eLife.54076`.
Open access, CC BY 4.0. **Table 3**, state-level age-standardised annual death rates per 100,000
(2010–2014 column) and estimated deaths 2001–14 in thousands.

Two rows in the source are aggregates — "Northeastern states" and "All other states" — and districts
in those states inherit the aggregate rate, flagged `burden_is_aggregate`.

---

## 5. District boundaries and the burden allocator — `data/curated/districts.csv`

594 district polygons from the `geohacker/india` open GeoJSON dataset
(`https://raw.githubusercontent.com/geohacker/india/master/district/india_district.geojson`),
derived from GADM. Centroids and areas are computed from the polygon rings by
`data/curated` build; the raw GeoJSON is not committed.

**The within-state allocator is a stated assumption, not a measurement.** A state's estimated deaths
are split among its districts in proportion to **district area**, used as a proxy for rural
population, because snakebite mortality in India is overwhelmingly rural and 2011 census
district-level rural populations were not retrievable here without an API key. The proxy is declared
on every row (`rural_population_basis`) and the sensitivity analysis re-runs the optimiser under a
uniform-within-state split to measure how much the choice matters. Nothing in this project claims
district-level snakebite mortality is known.

---

## 6. Antivenom calibration data — `data/curated/antivenomics_calibration.csv`

Only non-holdout studies. `model/calibrate.py` raises `HoldoutLeakError` on any row labelled
`holdout=True` **or** citing one of the four holdout DOIs.

| Study | DOI | Table | What was taken |
|---|---|---|---|
| Senji Laxme 2021 *Naja* | `10.1371/journal.pntd.0009150` | S3B Table | Murine ED50 neutralising potencies (mg/mL) for Punjab, West Bengal and the complete Rajasthan failure, recorded as a censored bound |
| Senji Laxme 2021 *Daboia* | `10.1371/journal.pntd.0009247` | S6 Table | Murine ED50 potencies for West Bengal, Maharashtra and Madhya Pradesh |
| Senji Laxme 2021 *Daboia* | `10.1371/journal.pntd.0009247` | S1 Table | Premium Serums antivenom protein content, 26.2 mg/mL, batch ASVS(I)-Lyo013 |
| Kalita, Patra & Mukherjee 2017, *Sci. Rep.* | `10.1038/s41598-017-17227-y` | Table 3 | Which toxin families are poorly recognised by commercial Indian polyvalent antivenoms, used as an ordinal constraint and never as a numeric recognition fraction |

---

## 7. Fixed model constants

| Constant | Value | Source |
|---|---|---|
| Vial protein content | 98 mg (95% CI 39–125) | `10.3390/toxins18010054`, Discussion, citing its ref 29: a Bharat Serums vial holds 390 mg lyophilisate at 25.2% protein |
| Standard initial dose | 10 vials | Indian national snakebite protocol |
| Delivered venom dose `D` | 40 mg, swept 10–100 mg | Central estimate for a significant envenoming, distinct from the much larger milked yield |
| Immunogen source | 12.6819 N, 80.0000 E | Irula Snake Catchers' Industrial Cooperative Society, Chengalpattu district, Tamil Nadu |
| `sigma_f` | `data/curated/severity_weights.csv` | Fixed a priori from the toxin-pathology literature with a per-family written rationale. Never fitted. |
| `kappa_f` shape | `M_ab / M_family` | Molar equivalence. Family masses in `config.FAMILY_MASS_KDA`. Only a single global scale is fitted. |

---

## 8. Holdout sources — used for compositions, never for fitting

These four studies supply venom compositions, which are legitimate model inputs, but none of their
antivenom results may enter a parameter fit. The guard is an assertion in code with tests.

| Target | DOI |
|---|---|
| R1 | `10.3390/toxins18010054` |
| R2 | `10.3389/fphar.2021.768210` |
| R3 | `10.1371/journal.pntd.0007899` |
| R4 | `10.1371/journal.pntd.0009659` |

---

## Reproducing the raw layer

```bash
python -m venomgap.cli fetch-sequences     # UniProt/ToxProt -> data/raw/uniprot/
python -m venomgap.cli fetch-occurrences   # GBIF            -> data/raw/gbif/
python -m venomgap.cli build-crossreact    # identity table  -> data/derived/
```

With `VENOMGAP_OFFLINE=1` set, no fetcher touches the network and any missing cache raises
`MissingCacheError`. CI runs with it set, so the committed derived tables are provably sufficient.
