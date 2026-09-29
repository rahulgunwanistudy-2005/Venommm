import { useEffect, useMemo, useState } from "react";
import {
  CartesianGrid, Legend, Line, LineChart, ReferenceLine, ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { fetchAtlas, fetchDistrict, fetchSiting } from "./api";
import { DeficitMap, type Metric } from "./DeficitMap";
import { DistrictDetailPanel } from "./DistrictDetail";
import type { AtlasPayload, DistrictDetail, SitingPayload } from "./types";

export function Atlas() {
  const [atlas, setAtlas] = useState<AtlasPayload | null>(null);
  const [siting, setSiting] = useState<SitingPayload | null>(null);
  const [k, setK] = useState(0);
  const [metric, setMetric] = useState<Metric>("deficit");
  const [selected, setSelected] = useState<string | null>(null);
  const [detail, setDetail] = useState<DistrictDetail | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    Promise.all([fetchAtlas(), fetchSiting()])
      .then(([a, s]) => { setAtlas(a); setSiting(s); })
      .catch((e: Error) => setError(e.message));
  }, []);

  useEffect(() => {
    if (!selected) { setDetail(null); return; }
    setDetail(null);
    fetchDistrict(selected).then(setDetail).catch((e: Error) => setError(e.message));
  }, [selected]);

  const solution = useMemo(
    () => siting?.solutions.find((s) => s.k === k) ?? null,
    [siting, k],
  );
  // Pins show where k centres would go if exactly k are built, which is what the slider asks.
  const forcedAtK = useMemo(
    () => siting?.forced_k_curve?.find((row) => row.k === k) ?? null,
    [siting, k],
  );

  if (error) return <div className="panel"><strong>Could not load results.</strong><p className="note">{error}</p></div>;
  if (!atlas || !siting) return <div className="panel">Loading district results…</div>;

  const pins = forcedAtK?.pins ?? solution?.pins ?? [];
  const coverage = forcedAtK ? forcedAtK.weighted : atlas.national_coverage;

  // Both curves matter and they answer different questions: "exactly k added, vial split evenly"
  // is where the dilution penalty shows, and "weights re-optimised" is what removes it.
  // Guard rather than assume: an older results file predates forced_k_curve, and a chart that
  // throws takes the whole screen with it.
  const curve = (siting.forced_k_curve ?? []).map((row) => ({
    k: row.k, uniform: row.uniform, weighted: row.weighted,
  }));

  return (
    <>
      <div className="grid three" style={{ marginBottom: 14 }}>
        <div className="panel tile">
          <div className="label">Burden-weighted national coverage</div>
          <div className="value">{(coverage * 100).toFixed(1)}%</div>
          <div className="sub">
            current Big Four immunogen: {(siting.baseline_national_coverage * 100).toFixed(1)}%
            {k > 0 && forcedAtK &&
              ` · ${forcedAtK.weighted >= siting.baseline_national_coverage ? "+" : ""}` +
              `${((forcedAtK.weighted - siting.baseline_national_coverage) * 100).toFixed(1)} points`}
          </div>
        </div>
        <div className="panel tile">
          <div className="label">Districts modelled</div>
          <div className="value">{atlas.estimated}</div>
          <div className="sub">{atlas.unknown} reported unknown, rendered grey</div>
        </div>
        <div className="panel tile">
          <div className="label">Fitted spatial length scale</div>
          <div className="value">{atlas.ell_km.toFixed(0)} km</div>
          <div className="sub">
            plateau {atlas.ell_plateau_km[0].toFixed(0)}–{atlas.ell_plateau_km[1].toFixed(0)} km
          </div>
        </div>
      </div>

      <div className="controls panel">
        <label>
          New collection sites k = <strong>{k}</strong>
          <input type="range" min={0} max={6} step={1} value={k}
                 onChange={(e) => setK(Number(e.target.value))} style={{ width: 200 }} />
        </label>
        <label>
          Show
          <select value={metric} onChange={(e) => setMetric(e.target.value as Metric)}>
            <option value="deficit">deficit</option>
            <option value="uncertainty">uncertainty</option>
            <option value="burden">burden</option>
          </select>
        </label>
        {forcedAtK && k > 0 && (
          <span className="note">
            {forcedAtK.pins.length} site{forcedAtK.pins.length === 1 ? "" : "s"} ·
            {" "}uniform vial weights {(forcedAtK.uniform * 100).toFixed(1)}% ·
            {" "}re-optimised {(forcedAtK.weighted * 100).toFixed(1)}%
            {solution && ` · ${solution.districts_moved_out_of_high_deficit} districts out of high deficit`}
          </span>
        )}
      </div>

      <div className="grid two" style={{ marginTop: 14 }}>
        <div className="panel">
          <DeficitMap districts={atlas.districts} metric={metric} pins={pins}
                      selected={selected} onSelect={setSelected} />
          <p className="note" style={{ marginTop: 10 }}>
            Click a district for its per-family breakdown, provenance and flags.
          </p>
        </div>
        <div>
          <div className="panel" style={{ marginBottom: 14 }}>
            <h3 style={{ margin: "0 0 4px", fontSize: 14 }}>National coverage vs k</h3>
            <p className="note" style={{ marginTop: 0 }}>
              Turnover at k = {siting.turnover_k}: past it, coverage falls as the next venoms
              dilute the vial. Re-optimising the mixture weights removes the dip. Measured
              optimality gap {siting.optimality_gap.optimality_gap.toFixed(4)} against exact
              enumeration of {siting.optimality_gap.subsets_evaluated.toLocaleString()} subsets.
            </p>
            <ResponsiveContainer width="100%" height={190}>
              <LineChart data={curve} margin={{ top: 6, right: 12, left: 0, bottom: 4 }}>
                <CartesianGrid stroke="var(--border)" vertical={false} />
                <XAxis dataKey="k" stroke="var(--text-muted)" fontSize={11}
                       label={{ value: "new collection sites k", position: "insideBottom", offset: -2, fontSize: 11 }} />
                <YAxis stroke="var(--text-muted)" fontSize={11} domain={["auto", "auto"]}
                       tickFormatter={(v: number) => v.toFixed(2)} />
                <Legend fontSize={11} />
                <Tooltip formatter={(v: number) => v.toFixed(4)}
                         contentStyle={{ background: "var(--surface-1)", border: "1px solid var(--border)", fontSize: 12 }} />
                <ReferenceLine x={siting.turnover_k} stroke="var(--critical)" strokeDasharray="4 3" />
                <ReferenceLine x={k} stroke="var(--series-3)" />
                <Line type="monotone" dataKey="uniform" name="uniform vial weights"
                    stroke="var(--series-1)" strokeWidth={2} dot={{ r: 3 }} />
              <Line type="monotone" dataKey="weighted" name="re-optimised weights"
                    stroke="var(--series-3)" strokeWidth={2} dot={{ r: 3 }} />
              </LineChart>
            </ResponsiveContainer>
            <p className="note">{siting.weight_reoptimisation_note}</p>
            <p className="note">{siting.submodularity_note}</p>
          </div>
          <div className="panel">
            <h3 style={{ margin: "0 0 6px", fontSize: 14 }}>What the length scale means</h3>
            <p className="note" style={{ margin: 0 }}>{atlas.interpretation}</p>
          </div>
        </div>
      </div>

      {selected && (
        <div style={{ marginTop: 14 }}>
          {detail ? <DistrictDetailPanel detail={detail} /> : <div className="panel">Loading district…</div>}
        </div>
      )}
    </>
  );
}
