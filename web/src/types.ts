export interface District {
  district_id: string;
  district: string;
  state: string;
  burden_state: string;
  lat: number;
  lon: number;
  coverage: number;
  deficit: number;
  uncertainty: number;
  status: "estimated" | "unknown";
  burden_weight: number;
  nearest_proteome_km: number;
  dominant_species: string[];
  limiting_families: string[];
  flags: string[];
}

export interface AtlasPayload {
  districts: District[];
  national_coverage: number;
  national_deficit: number;
  ell_km: number;
  ell_plateau_km: [number, number];
  unknown_cutoff_km: number;
  estimated: number;
  unknown: number;
  interpretation: string;
  disclaimer: string;
}

export interface Pin {
  pop_id: string;
  species: string;
  locality: string;
  state: string;
  lat: number;
  lon: number;
  weight: number;
}

export interface Solution {
  k: number;
  sites: string[];
  weights: number[];
  national_coverage: number;
  coverage_gain_vs_baseline: number;
  districts_moved_out_of_high_deficit: number;
  solver: string;
  pins: Pin[];
}

export interface SitingPayload {
  baseline_national_coverage: number;
  solutions: Solution[];
  coverage_vs_k: { k: number; greedy: number; greedy_local: number; weight_optimised: number }[];
  forced_k_curve: {
    k: number; uniform: number; weighted: number; sites: string[]; pins: Pin[];
  }[];
  turnover_k: number;
  turnover_basis: string;
  weight_reoptimisation_note: string;
  global_best_k_uniform: number;
  optimality_gap: { optimality_gap: number; k: number; subsets_evaluated: number };
  submodularity_note: string;
}

export interface SpeciesDetail {
  species: string;
  coverage: number;
  deficit: number;
  bite_share: number;
  uncertainty: number;
  nearest_proteome_km: number;
  status: string;
  composition: Record<string, number>;
  per_family_neutralised: Record<string, number>;
  limiting_families: string[];
  sequence_imputed: boolean;
}

export interface Provenance {
  pop_id: string;
  species: string;
  locality: string;
  state: string;
  doi: string;
  table: string;
  method: string;
  accessed: string;
  flags: string[];
  note: string;
}

export interface DistrictDetail extends District {
  per_species: SpeciesDetail[];
  provenance: Provenance[];
  flag_explanations: Record<string, string>;
}

export interface Population {
  pop_id: string;
  species: string;
  locality: string;
  state: string;
  country: string;
  lat: number;
  lon: number;
  dominant_family: string;
  composition: Record<string, number>;
  doi: string;
  table: string;
  flags: string[];
  holdout: boolean;
}

export interface MixtureResult {
  pop_ids: string[];
  weights: number[];
  national_coverage: number;
  national_deficit: number;
  baseline_national_coverage: number;
  change_vs_baseline: number;
  mixture_size: number;
}
