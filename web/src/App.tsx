import { useState } from "react";
import { Atlas } from "./Atlas";
import { MixtureExplorer } from "./MixtureExplorer";
import { Retrodiction } from "./Retrodiction";
import "./styles.css";

const TABS = [
  { id: "atlas", label: "Atlas" },
  { id: "mixture", label: "Mixture explorer" },
  { id: "retrodiction", label: "Retrodiction" },
] as const;

type TabId = (typeof TABS)[number]["id"];

export default function App() {
  const [tab, setTab] = useState<TabId>("atlas");
  return (
    <div className="shell">
      <header className="app">
        <h1>VenomGap</h1>
        <p>
          India&rsquo;s antivenom is an estimator trained on venom from one collection locality.
          This is what that costs, district by district, and where to sample next.
        </p>
      </header>

      <nav className="tabs">
        {TABS.map((t) => (
          <button key={t.id} aria-current={tab === t.id} onClick={() => setTab(t.id)}>
            {t.label}
          </button>
        ))}
      </nav>

      {tab === "atlas" && <Atlas />}
      {tab === "mixture" && <MixtureExplorer />}
      {tab === "retrodiction" && <Retrodiction />}

      <footer className="disclaimer">
        Research model. Not clinical guidance. Predictions are computational and require
        experimental validation. The model predicts a risk ordering, not clinical neutralisation
        percentages.
      </footer>
    </div>
  );
}
