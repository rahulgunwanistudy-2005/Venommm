import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  Bar, BarChart, CartesianGrid, Cell, ReferenceLine, ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { evaluateMixture, fetchPopulations } from "./api";
import type { MixtureResult, Population } from "./types";

/**
 * Screen 3. Check venoms in and out of the immunising mixture and watch national coverage move --
 * including watching it go *down* when the mixture gets too large. That is the interactive proof
 * of R4, and the mechanism is a fixed antibody budget, not a fitted penalty on mixture size.
 */
export function MixtureExplorer() {
  const [populations, setPopulations] = useState<Population[]>([]);
  const [selected, setSelected] = useState<string[]>([]);
  const [result, setResult] = useState<MixtureResult | null>(null);
  const [history, setHistory] = useState<{ size: number; coverage: number; label: string }[]>([]);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [filter, setFilter] = useState("");
  const pending = useRef(0);

  useEffect(() => {
    fetchPopulations()
      .then(({ populations: rows }) => {
        setPopulations(rows);
        // Start from one venom per Big Four species so the first move is meaningful.
        const seed: string[] = [];
        for (const species of ["Naja naja", "Bungarus caeruleus", "Daboia russelii", "Echis carinatus"]) {
          const match = rows.find((r) => r.species === species && r.country === "India");
          if (match) seed.push(match.pop_id);
        }
        setSelected(seed);
      })
      .catch((e: Error) => setError(e.message));
  }, []);

  const evaluate = useCallback((popIds: string[]) => {
    if (popIds.length === 0) { setResult(null); return; }
    const ticket = ++pending.current;
    setBusy(true);
    evaluateMixture(popIds)
      .then((r) => {
        if (ticket !== pending.current) return;  // a newer request has superseded this one
        setResult(r);
        setHistory((h) => [
          ...h.slice(-19),
          { size: r.mixture_size, coverage: r.national_coverage, label: `|S| = ${r.mixture_size}` },
        ]);
      })
      .catch((e: Error) => setError(e.message))
      .finally(() => { if (ticket === pending.current) setBusy(false); });
  }, []);

  useEffect(() => { evaluate(selected); }, [selected, evaluate]);

  const toggle = (popId: string) =>
    setSelected((current) =>
      current.includes(popId) ? current.filter((p) => p !== popId) : [...current, popId],
    );

  const visible = useMemo(() => {
    const needle = filter.trim().toLowerCase();
    if (!needle) return populations;
    return populations.filter((p) =>
      `${p.species} ${p.locality} ${p.state} ${p.country} ${p.dominant_family}`
        .toLowerCase().includes(needle),
    );
  }, [populations, filter]);

  const delta = result ? result.change_vs_baseline : 0;
  const lastTwo = history.slice(-2);
  const justFell = lastTwo.length === 2 && lastTwo[1].coverage < lastTwo[0].coverage - 1e-9
    && lastTwo[1].size > lastTwo[0].size;

  if (error) return <div className="panel"><strong>Mixture explorer unavailable.</strong><p className="note">{error}</p></div>;

  return (
    <div className="grid two">
      <div className="panel">
        <h3 style={{ margin: "0 0 4px", fontSize: 15 }}>Immunising mixture</h3>
        <p className="note" style={{ marginTop: 0 }}>
          Add venoms and watch burden-weighted national coverage. Past the optimum it falls: a vial
          holds a fixed mass of antibody, so raising antibody against one more venom is antibody not
          raised against the toxins a given bite actually delivers.
        </p>
        <input value={filter} onChange={(e) => setFilter(e.target.value)}
               placeholder="Filter by species, state or dominant family"
               style={{ width: "100%", padding: "7px 9px", marginBottom: 8,
                        background: "var(--surface-2)", color: "var(--text-primary)",
                        border: "1px solid var(--border)", borderRadius: 6, font: "inherit" }} />
        <div className="pop-list">
          {visible.map((p) => (
            <div className="pop-row" key={p.pop_id}>
              <input type="checkbox" id={p.pop_id} checked={selected.includes(p.pop_id)}
                     onChange={() => toggle(p.pop_id)} />
              <label htmlFor={p.pop_id}>
                <em>{p.species}</em> — {p.locality}
                <div className="meta">
                  {p.state}, {p.country} · dominant {p.dominant_family} · {p.doi}
                  {p.holdout && <span className="flag warn" style={{ marginLeft: 6 }}>holdout target</span>}
                  {p.flags.map((f) => <span key={f} className="flag">{f.replace(/_/g, " ")}</span>)}
                </div>
              </label>
            </div>
          ))}
        </div>
      </div>

      <div>
        <div className="panel tile" style={{ marginBottom: 14 }}>
          <div className="label">Burden-weighted national coverage</div>
          <div className="value" style={{ color: delta < -1e-9 ? "var(--critical)" : undefined }}>
            {result ? `${(result.national_coverage * 100).toFixed(1)}%` : "—"}
            {busy && <span className="sub" style={{ marginLeft: 8 }}>updating…</span>}
          </div>
          <div className="sub">
            {result
              ? `${result.mixture_size} venom${result.mixture_size === 1 ? "" : "s"} · ` +
                `${delta >= 0 ? "+" : ""}${(delta * 100).toFixed(1)} points vs the current Big Four immunogen`
              : "select at least one venom"}
          </div>
          {justFell && (
            <div className="callout" style={{ borderLeftColor: "var(--critical)" }}>
              Coverage just <strong>fell</strong> as a venom was added. The model was never told this
              happens — it follows from a finite antibody budget, and it is what an experimental
              study reported in 2021.
            </div>
          )}
        </div>

        <div className="panel">
          <h3 style={{ margin: "0 0 6px", fontSize: 14 }}>Your moves so far</h3>
          <ResponsiveContainer width="100%" height={240}>
            <BarChart data={history.map((h, i) => ({ ...h, i }))}
                      margin={{ top: 8, right: 12, left: 0, bottom: 4 }}>
              <CartesianGrid stroke="var(--border)" vertical={false} />
              <XAxis dataKey="label" stroke="var(--text-muted)" fontSize={11} interval={0}
                     angle={-30} textAnchor="end" height={48} />
              <YAxis stroke="var(--text-muted)" fontSize={11} domain={[0, 1]}
                     tickFormatter={(v: number) => v.toFixed(1)}
                     label={{ value: "national coverage (fraction)", angle: -90,
                              position: "insideLeft", fontSize: 11 }} />
              <Tooltip formatter={(v: number) => v.toFixed(4)}
                       contentStyle={{ background: "var(--surface-1)", border: "1px solid var(--border)", fontSize: 12 }} />
              {result && (
                <ReferenceLine y={result.baseline_national_coverage} stroke="var(--text-muted)"
                               strokeDasharray="4 3"
                               label={{ value: "current Big Four", fontSize: 10, position: "right" }} />
              )}
              <Bar dataKey="coverage" radius={[4, 4, 0, 0]}>
                {history.map((h, i) => (
                  <Cell key={i}
                        fill={i > 0 && h.coverage < history[i - 1].coverage - 1e-9 && h.size > history[i - 1].size
                          ? "var(--series-2)" : "var(--series-1)"} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
          <p className="note">
            Orange bars are moves where adding a venom lowered coverage.{" "}
            <button onClick={() => setHistory([])}
                    style={{ background: "none", border: 0, color: "var(--series-1)",
                             cursor: "pointer", font: "inherit", padding: 0 }}>
              Clear history
            </button>
          </p>
        </div>
      </div>
    </div>
  );
}
