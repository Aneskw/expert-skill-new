#!/usr/bin/env python3
"""Recover the current prepared dataset from this run's generated artifacts.

This is intentionally limited to the current run.  The v3.1 standard bundle is
the authority, while experiments.sqlite supplies the already-generated legacy
compatibility rows.  The standard merge command is run after this script to
regenerate all compatibility views through the normal tool path.
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path


RUN = Path(__file__).resolve().parents[1]
SUBMITTED = RUN / "submitted.json"
DATABASE = RUN / "experiments.sqlite"
DATASET = RUN / "dataset.json"


standard = json.loads(SUBMITTED.read_text())
reaction_order = [reaction["reaction_id"] for reaction in standard["reactions"]]

with sqlite3.connect(DATABASE) as connection:
    rows = connection.execute(
        "SELECT record_id, content FROM experiments"
    ).fetchall()

legacy_by_id = {record_id: json.loads(content) for record_id, content in rows}
missing = [record_id for record_id in reaction_order if record_id not in legacy_by_id]
extra = sorted(set(legacy_by_id) - set(reaction_order))
if missing or extra:
    raise RuntimeError(f"Compatibility-record mismatch: missing={missing}, extra={extra}")

dataset = {
    "tool_version": "3.1",
    "schema_version": "1.0",
    "created_at": "2026-10-18T13:58:25.393533+00:00",
    "documents": standard["documents"],
    "candidates": [],
    "experiments": [legacy_by_id[record_id] for record_id in reaction_order],
    "run_log": [
        {
            "event": "dataset_recovered_from_current_run_artifacts",
            "reason": "A missing jq executable left dataset.json empty during a targeted source_context deduplication attempt.",
            "authority": "runs/current/submitted.json",
            "compatibility_source": "runs/current/experiments.sqlite",
        }
    ],
    "standard": standard,
}

temporary = DATASET.with_suffix(".json.tmp")
temporary.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
temporary.replace(DATASET)
print(json.dumps({"recovered_reactions": len(reaction_order), "documents": len(standard["documents"])}))
