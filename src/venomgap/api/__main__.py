"""Run the API: `python -m venomgap.api`.

Binds the port given in `PORT`, falling back to `DEFAULT_API_PORT`. Reading the environment
rather than hard-coding a CLI flag is what lets a supervisor place this server on a free port
instead of colliding with whatever else is already listening. 8000 is a crowded default on a
development machine, so it is deliberately not the fallback here.
"""

from __future__ import annotations

import os

import uvicorn

DEFAULT_API_PORT = 8010
DEFAULT_API_HOST = "127.0.0.1"


def main() -> None:
    port = int(os.environ.get("PORT", DEFAULT_API_PORT))
    host = os.environ.get("HOST", DEFAULT_API_HOST)
    uvicorn.run("venomgap.api.app:app", host=host, port=port, log_level="info")


if __name__ == "__main__":
    main()
