#!/usr/bin/env python3
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools

dataset = reaction_tools.read(ROOT / "runs/current/dataset.json")
reaction_tools.write(ROOT / "runs/current/raw_model/audit.json", reaction_tools.dataset_audit(dataset))
