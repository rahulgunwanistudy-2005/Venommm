"""FastAPI app serving the atlas. Reads precomputed results; recomputes only what the user changes.

The mixture explorer needs live evaluation -- that is the point of it -- so the model is held in
memory and scored on demand. Everything else is served from `results/*.json` so the web app and the
paper cannot disagree.
"""

from __future__ import annotations

import logging
from typing import Any

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from venomgap.api.routes import router
from venomgap.config import DISCLAIMER, RESULTS_DIR

logger = logging.getLogger(__name__)


def create_app() -> FastAPI:
    app = FastAPI(
        title="VenomGap",
        version="0.1.0",
        description=(
            "Stoichiometric antivenom coverage modelling for India. " + DISCLAIMER
        ),
    )
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
        allow_methods=["GET", "POST"],
        allow_headers=["*"],
    )
    app.include_router(router)

    @app.get("/health")
    def health() -> dict[str, Any]:
        missing = [
            name
            for name in (
                "spatial.json", "optimisation.json", "retrodiction.json", "sensitivity.json",
            )
            if not (RESULTS_DIR / name).exists()
        ]
        return {
            "ok": not missing,
            "missing_results": missing,
            "disclaimer": DISCLAIMER,
        }

    return app


app = create_app()
