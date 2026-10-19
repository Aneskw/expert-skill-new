#!/usr/bin/env python3
"""Persist the deterministic advisory audit as a raw workflow artifact."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools

dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
reaction_tools.write(ROOT / "runs/current/raw_model/audit.json",
                     reaction_tools.dataset_audit(dataset))
