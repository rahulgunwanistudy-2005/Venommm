import type { DistrictDetail } from "./types";

/** Screen 2: why this district is where it is, in plain words, with every flag explained. */
export function DistrictDetailPanel({ detail }: { detail: DistrictDetail }) {
  const flags = Array.from(new Set([...detail.flags, ...detail.per_species.flatMap(
    (s) => (s.sequence_imputed ? ["sequence_imputed"] : []),
  )]));

  return (
    <div className="panel">
      <header style={{ display: "flex", flexWrap: "wrap", gap: 16, alignItems: "baseline" }}>
        <h2 style={{ margin: 0, fontSize: 17 }}>{detail.district}, {detail.state}</h2>
        <span className="note">
          coverage {(detail.coverage * 100).toFixed(1)}% · deficit {(detail.deficit * 100).toFixed(1)}%
          {" · "}uncertainty {detail.uncertainty.toFixed(2)}
          {" · "}status {detail.status}
        </span>
      </header>

      {detail.status === "unknown" && (
        <div className="callout">
          The nearest published venom proteome is {detail.nearest_proteome_km.toFixed(0)} km away,
          beyond the distance this model is willing to extrapolate. This district is reported as
          <strong> unknown</strong>, not as covered. Unknown is a finding, not a blank.
        </div>
      )}

      <div className="grid two" style={{ marginTop: 12 }}>
        <div>
          <h3 style={{ fontSize: 14, margin: "0 0 6px" }}>Per species, and what limits each</h3>
          <table>
            <thead>
              <tr>
                <th>Species</th>
                <th className="num">Coverage</th>
                <th className="num">Bite share</th>
                <th>Limiting families (worst first)</th>
                <th className="num">Nearest proteome</th>
              </tr>
            </thead>
            <tbody>
              {detail.per_species.map((row) => (
                <tr key={row.species}>
                  <td><em>{row.species}</em></td>
                  <td className="num">
                    <div className="bar-track" style={{ width: 70, display: "inline-block", marginRight: 6 }}>
                      <div className="bar-fill"
                           style={{ width: `${Math.max(2, row.coverage * 100)}%`, background: "var(--series-1)" }} />
                    </div>
                    {(row.coverage * 100).toFixed(0)}%
                  </td>
                  <td className="num">{row.bite_share.toFixed(2)}</td>
                  <td>{row.limiting_families.length ? row.limiting_families.join(" · ") : "none"}</td>
                  <td className="num">{row.nearest_proteome_km.toFixed(0)} km</td>
                </tr>
              ))}
            </tbody>
          </table>

          {detail.per_species[0] && (
            <>
              <h3 style={{ fontSize: 14, margin: "16px 0 6px" }}>
                Per-family neutralised fraction — <em>{detail.per_species[0].species}</em>
              </h3>
              <table>
                <thead>
                  <tr><th>Family</th><th className="num">Abundance</th><th className="num">Neutralised</th></tr>
                </thead>
                <tbody>
                  {Object.entries(detail.per_species[0].composition)
                    .filter(([, v]) => v > 0.005)
                    .sort((a, b) => b[1] - a[1])
                    .map(([family, abundance]) => {
                      const n = detail.per_species[0].per_family_neutralised[family] ?? 1;
                      const limiting = detail.per_species[0].limiting_families.includes(family);
                      return (
                        <tr key={family}>
                          <td>{family}{limiting && <span className="flag warn" style={{ marginLeft: 6 }}>limiting</span>}</td>
                          <td className="num">{(abundance * 100).toFixed(1)}%</td>
                          <td className="num">
                            <div className="bar-track" style={{ width: 80, display: "inline-block", marginRight: 6 }}>
                              <div className="bar-fill" style={{
                                width: `${Math.max(2, n * 100)}%`,
                                background: limiting ? "var(--series-2)" : "var(--series-1)",
                              }} />
                            </div>
                            {(n * 100).toFixed(0)}%
                          </td>
                        </tr>
                      );
                    })}
                </tbody>
              </table>
            </>
          )}
        </div>

        <div>
          <h3 style={{ fontSize: 14, margin: "0 0 6px" }}>What is assumed here</h3>
          {flags.length === 0 && <p className="note">No caveats flagged for this district.</p>}
          {flags.map((flag) => (
            <p key={flag} className="note" style={{ margin: "0 0 8px" }}>
              <span className="flag warn">{flag.replace(/_/g, " ")}</span>
              <br />
              {detail.flag_explanations[flag] ?? "No explanation recorded for this flag."}
            </p>
          ))}

          <h3 style={{ fontSize: 14, margin: "16px 0 6px" }}>
            Proteomes informing this estimate ({detail.provenance.length})
          </h3>
          <table>
            <thead><tr><th>Population</th><th>Source</th></tr></thead>
            <tbody>
              {detail.provenance.map((p) => (
                <tr key={p.pop_id}>
                  <td>
                    <em>{p.species}</em><br />
                    <span className="note">{p.locality}, {p.state}</span>
                  </td>
                  <td>
                    <span className="note">{p.doi}</span><br />
                    <span className="note">{p.table} · {p.method} · accessed {p.accessed}</span>
                    {p.flags.map((f) => <span key={f} className="flag">{f.replace(/_/g, " ")}</span>)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}
