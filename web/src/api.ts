import type { AtlasPayload, DistrictDetail, MixtureResult, Population, SitingPayload } from "./types";

async function get<T>(path: string): Promise<T> {
  const response = await fetch(path);
  if (!response.ok) {
    throw new Error(`${path} returned ${response.status}: ${await response.text()}`);
  }
  return (await response.json()) as T;
}

export const fetchAtlas = () => get<AtlasPayload>("/api/atlas");
export const fetchSiting = () => get<SitingPayload>("/api/siting");
export const fetchDistrict = (id: string) =>
  get<DistrictDetail>(`/api/district/${encodeURIComponent(id)}`);
export const fetchPopulations = () =>
  get<{ populations: Population[]; families: string[] }>("/api/populations");

export async function evaluateMixture(
  popIds: string[],
  weights?: number[],
): Promise<MixtureResult> {
  const response = await fetch("/api/mixture", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ pop_ids: popIds, weights }),
  });
  if (!response.ok) {
    throw new Error(`mixture evaluation failed: ${await response.text()}`);
  }
  return (await response.json()) as MixtureResult;
}
