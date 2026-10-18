#!/usr/bin/env python3
import json
from pathlib import Path

run = Path(__file__).resolve().parents[1]
checks = {
    "table1_entry25": {
        "file": "records/doc_1a0071146ff7cf43__table1_entry25.md",
        "tokens": ["0.31 mmol", "0.6 equiv", "2 equiv", "r.t.", "1 h", "75%", "isolated yield", "merged rather than double-counted"],
    },
    "scheme3_4a": {
        "file": "records/doc_1a0071146ff7cf43__scheme3_4a.md",
        "tokens": ["1.5 equiv", "1.8 equiv", "DMSO", "0–5 °C", "rt", "0.5 h", "87%", "1:0.68", "Temperature conflict"],
    },
    "scheme4b": {
        "file": "records/doc_1a0071146ff7cf43__scheme4b_isolated_4a_cyclization.md",
        "tokens": ["0.17 mmol", "2 equiv", "acetonitrile", "r.t.", "1 h", "78%", "Yield measurement basis is unreported"],
    },
}
for row in checks.values():
    text = (run / row["file"]).read_text()
    row["missing_tokens"] = [token for token in row["tokens"] if token not in text]
(run / "raw_model/verification_summary.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2) + "\n")
if any(row["missing_tokens"] for row in checks.values()):
    raise SystemExit(1)
print(json.dumps({k: v["missing_tokens"] for k, v in checks.items()}, ensure_ascii=False))
