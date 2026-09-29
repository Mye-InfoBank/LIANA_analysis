#!/usr/bin/env python3

"""
AI Snapshot utilities for the LIANA Results Explorer.

Snapshots contain a machine-readable representation of the
current Shiny analysis state. They are stored temporarily inside
the running container and can be shared with external LLMs.

Snapshots:
- expire automatically after SNAPSHOT_TTL_DAYS
- disappear when the Docker container is rebuilt/replaced
- are limited in number to avoid uncontrolled storage growth
"""

import json
import time
import uuid
from pathlib import Path
from datetime import datetime, timezone


# ============================================================
# CONFIGURATION
# ============================================================

# Because the app runs with WORKDIR /app/shiny,
# Absolute path to:
# shiny/ai_snapshots
SNAPSHOT_DIR = (
    Path(__file__).resolve().parent.parent
    / "ai_snapshots"
)


# Links remain usable for up to x days.
SNAPSHOT_TTL_DAYS = 7

# Additional safety cap.
MAX_SNAPSHOTS = 100


SNAPSHOT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# CLEANUP
# ============================================================

def cleanup_old_snapshots():
    """
    Remove expired snapshots and enforce MAX_SNAPSHOTS.
    """

    now = time.time()

    max_age_seconds = (
        SNAPSHOT_TTL_DAYS
        * 24
        * 60
        * 60
    )

    snapshot_files = list(
        SNAPSHOT_DIR.glob("*.json")
    )

    # --------------------------------------------------------
    # Delete snapshots older than the TTL
    # --------------------------------------------------------

    for path in snapshot_files:

        try:

            age_seconds = (
                now
                - path.stat().st_mtime
            )

            if age_seconds > max_age_seconds:
                path.unlink()

        except Exception:
            # Cleanup problems should never crash the app.
            pass

    # --------------------------------------------------------
    # Enforce maximum number of snapshots
    # --------------------------------------------------------

    snapshot_files = list(
        SNAPSHOT_DIR.glob("*.json")
    )

    if len(snapshot_files) > MAX_SNAPSHOTS:

        snapshot_files.sort(
            key=lambda path:
                path.stat().st_mtime
        )

        number_to_delete = (
            len(snapshot_files)
            - MAX_SNAPSHOTS
        )

        for path in snapshot_files[
            :number_to_delete
        ]:

            try:
                path.unlink()
            except Exception:
                pass


# ============================================================
# PLOTLY SERIALIZATION
# ============================================================

def figure_to_dict(fig):
    """
    Convert a Plotly figure into a JSON-safe Python dictionary.

    This stores the actual information behind the visible plot:
    axes, labels, traces, values, hover data, etc.
    """

    if fig is None:
        return None

    try:
        return json.loads(
            fig.to_json()
        )

    except Exception:
        return None


# ============================================================
# SAVE
# ============================================================

def save_ai_snapshot(
    payload: dict
) -> str:
    """
    Save one AI snapshot.

    Returns
    -------
    str
        Unique snapshot ID.
    """

    cleanup_old_snapshots()

    snapshot_id = (
        uuid.uuid4().hex
    )

    payload["snapshot"] = {
        "id": snapshot_id,

        "created_at": datetime.now(
            timezone.utc
        ).isoformat(),

        "expires_after_days":
            SNAPSHOT_TTL_DAYS
    }

    output_path = (
        SNAPSHOT_DIR
        / f"{snapshot_id}.json"
    )

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        # Compact JSON saves disk space.
        json.dump(
            payload,
            file,
            ensure_ascii=False,
            separators=(",", ":"),
            default=str
        )

    return snapshot_id


# ============================================================
# LOAD
# ============================================================

def load_ai_snapshot(
    snapshot_id: str
):
    """
    Load a snapshot.

    Returns None if:
    - the ID is invalid
    - the snapshot does not exist
    - the snapshot has expired
    """

    if (
        not snapshot_id
        or not snapshot_id.isalnum()
    ):
        return None

    path = (
        SNAPSHOT_DIR
        / f"{snapshot_id}.json"
    )

    if not path.exists():
        return None

    age_seconds = (
        time.time()
        - path.stat().st_mtime
    )

    max_age_seconds = (
        SNAPSHOT_TTL_DAYS
        * 24
        * 60
        * 60
    )

    if age_seconds > max_age_seconds:

        try:
            path.unlink()
        except Exception:
            pass

        return None

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as file:

        return json.load(file)