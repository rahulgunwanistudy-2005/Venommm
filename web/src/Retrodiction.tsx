import { useEffect, useState } from "react";

interface Criterion {
  criterion_id: string;
  statement: string;
  computed: number | null;
  threshold: number | null;
  comparison: string;
  passed: boolean;
  detail: string;
}

interface Payload {
  targets_passed: number;
  targets_total: number;
  criteria: Criterion[];
  preregistration_commit: string;
  notes: string;
  diagnostics: {
    theta_scale_mismatch: { theta: number; fraction_below_theta: number; consequence: string };
    data_sparsity: { consequence: string };
  };
}

/** The blinded retrodiction, reported exactly as it came out. */
export function Retrodiction() {
  const [payload, setPayload] = useState<Payload | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetch("/api/retrodiction")
      .then((r) => (r.ok ? r.json() : Promise.reject(new Error(`${r.status}`))))
      .then(setPayload)
      .catch((e: Error) => setError(e.message));
  }, []);

  if (error) return <div className="panel">Retrodiction results unavailable: {error}</div>;
  if (!payload) return <div className="panel">Loading…</div>;

  return (
    <>
      <div className="panel" style={{ marginBottom: 14 }}>
        <h2 style={{ margin: "0 0 4px", fontSize: 17 }}>
          Blinded retrodiction: {payload.targets_passed} of {payload.targets_total} targets passed
        </h2>
        <p className="note" style={{ marginTop: 0 }}>
          Four published antivenom-failure findings were held out. Parameters were fitted without
          any of their antivenom results, and a guard in the calibration code raises rather than
          warns if a holdout row reaches a fit. Pre-registration commit{" "}
          <code>{payload.preregistration_commit.slice(0, 12)}</code>.
        </p>
        <div className="callout">{payload.notes}</div>
      </div>

      <div className="panel" style={{ marginBottom: 14 }}>
        <table>
          <thead>
            <tr>
              <th>Criterion</th><th>Pre-registered statement</th>
              <th className="num">Computed</th><th className="num">Threshold</th>
              <th>Outcome</th>
            </tr>
          </thead>
          <tbody>
            {payload.criteria.map((c) => (
              <tr key={c.criterion_id}>
                <td><strong>{c.criterion_id}</strong></td>
                <td style={{ maxWidth: 380 }}>
                  {c.statement}
                  <div className="note" style={{ marginTop: 3 }}>{c.detail}</div>
                </td>
                <td className="num">{c.computed === null ? "—" : c.computed.toFixed(4)}</td>
                <td className="num">{c.threshold === null ? "—" : c.threshold.toFixed(4)}</td>
                <td className={c.passed ? "pass" : "fail"}>{c.passed ? "PASS" : "FAIL"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="panel">
        <h3 style={{ margin: "0 0 6px", fontSize: 15 }}>Why the misses missed</h3>
        <p className="note">{payload.diagnostics.theta_scale_mismatch.consequence}</p>
        <p className="note">{payload.diagnostics.data_sparsity.consequence}</p>
      </div>
    </>
  );
}
