import { useMemo } from "react";
import type { District } from "./types";

/**
 * District choropleth drawn as inline SVG from district centroids.
 *
 * Deliberately not a tile-based map: the demo has to run with no network at all, and a tile server
 * is a network dependency. The information content is the same -- one mark per district, positioned
 * geographically, coloured by the selected metric -- and it renders identically offline.
 *
 * The `unknown` class is grey and is drawn in its own pass, so a data-free district can never pick
 * up a colour from the ramp.
 */

const RAMP = [
  "var(--ramp-0)", "var(--ramp-1)", "var(--ramp-2)",
  "var(--ramp-3)", "var(--ramp-4)", "var(--ramp-5)", "var(--ramp-6)",
];

export type Metric = "deficit" | "uncertainty" | "burden";

const METRIC_LABEL: Record<Metric, string> = {
  deficit: "predicted deficit (fraction of weighted venom dose unneutralised)",
  uncertainty: "estimate uncertainty (0 = tight, 1 = no usable nearby proteome)",
  burden: "share of national snakebite burden (relative, rescaled to the maximum district)",
};

function rampColour(value: number): string {
  const index = Math.min(RAMP.length - 1, Math.max(0, Math.floor(value * RAMP.length)));
  return RAMP[index];
}

interface Props {
  districts: District[];
  metric: Metric;
  pins: { lat: number; lon: number; pop_id: string; species: string }[];
  selected: string | null;
  onSelect: (districtId: string) => void;
}

export function DeficitMap({ districts, metric, pins, selected, onSelect }: Props) {
  const { width, height, project, maxBurden } = useMemo(() => {
    const lons = districts.map((d) => d.lon);
    const lats = districts.map((d) => d.lat);
    const minLon = Math.min(...lons) - 0.6;
    const maxLon = Math.max(...lons) + 0.6;
    const minLat = Math.min(...lats) - 0.6;
    const maxLat = Math.max(...lats) + 0.6;
    const w = 760;
    const h = Math.round((w * (maxLat - minLat)) / (maxLon - minLon));
    return {
      width: w,
      height: h,
      maxBurden: Math.max(...districts.map((d) => d.burden_weight), 1e-9),
      project: (lat: number, lon: number): [number, number] => [
        ((lon - minLon) / (maxLon - minLon)) * w,
        h - ((lat - minLat) / (maxLat - minLat)) * h,
      ],
    };
  }, [districts]);

  const valueOf = (d: District): number => {
    if (metric === "deficit") return d.deficit;
    if (metric === "uncertainty") return d.uncertainty;
    return d.burden_weight / maxBurden;
  };

  const unknown = districts.filter((d) => d.status === "unknown");
  const known = districts.filter((d) => d.status !== "unknown");

  return (
    <div className="svg-wrap">
      <svg viewBox={`0 0 ${width} ${height}`} width="100%" role="img"
           aria-label={`District map coloured by ${METRIC_LABEL[metric]}`}>
        {unknown.map((d) => {
          const [x, y] = project(d.lat, d.lon);
          return (
            <rect key={d.district_id} x={x - 3.4} y={y - 3.4} width={6.8} height={6.8} rx={1.4}
                  fill="var(--unknown)" onClick={() => onSelect(d.district_id)}
                  style={{ cursor: "pointer" }}>
              <title>{`${d.district}, ${d.state} — unknown (nearest proteome beyond the extrapolation limit)`}</title>
            </rect>
          );
        })}
        {known.map((d) => {
          const [x, y] = project(d.lat, d.lon);
          const isSelected = d.district_id === selected;
          return (
            <rect key={d.district_id} x={x - 3.4} y={y - 3.4} width={6.8} height={6.8} rx={1.4}
                  fill={rampColour(valueOf(d))}
                  stroke={isSelected ? "var(--critical)" : "none"} strokeWidth={isSelected ? 2 : 0}
                  onClick={() => onSelect(d.district_id)} style={{ cursor: "pointer" }}>
              <title>
                {`${d.district}, ${d.state}\ndeficit ${d.deficit.toFixed(3)} · uncertainty ${d.uncertainty.toFixed(2)}` +
                 `\nnearest proteome ${d.nearest_proteome_km.toFixed(0)} km` +
                 (d.dominant_species.length ? `\nspecies: ${d.dominant_species.join(", ")}` : "")}
              </title>
            </rect>
          );
        })}
        {(() => {
          const [x, y] = project(12.6819, 80.0);
          return (
            <g key="immunogen-source">
              <circle cx={x} cy={y} r={7} fill="none" stroke="var(--critical)" strokeWidth={2} />
              <circle cx={x} cy={y} r={2.4} fill="var(--critical)" />
              <title>Immunogen source: Irula Cooperative, Chengalpattu, Tamil Nadu</title>
            </g>
          );
        })()}
        {pins.map((p) => {
          const [x, y] = project(p.lat, p.lon);
          return (
            <g key={p.pop_id}>
              <path d={`M ${x} ${y} l -6 -13 a 6.8 6.8 0 1 1 12 0 z`} fill="var(--series-3)"
                    stroke="var(--surface-1)" strokeWidth={1.4} />
              <title>{`Proposed collection site: ${p.species} — ${p.pop_id}`}</title>
            </g>
          );
        })}
      </svg>
      <div className="legend" style={{ marginTop: 8 }}>
        <span>{METRIC_LABEL[metric]}</span>
        {RAMP.map((c, i) => (
          <span key={c}>
            <span className="swatch" style={{ background: c }} />
            {(i / RAMP.length).toFixed(2)}
          </span>
        ))}
        <span>
          <span className="swatch" style={{ background: "var(--unknown)" }} />
          unknown — no nearby proteome
        </span>
        <span>
          <span className="swatch" style={{ background: "var(--series-3)" }} />
          proposed collection site
        </span>
      </div>
    </div>
  );
}
