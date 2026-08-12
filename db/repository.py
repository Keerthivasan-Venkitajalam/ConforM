"""Scientific memory repository.

Backend selection:
  CONFORM_DB_URL unset (default)   -> SQLite at artifacts/conform_memory.db
  CONFORM_DB_URL=postgresql://...  -> PostgreSQL via psycopg (if installed)

The scientific pipeline must never become unusable because PostgreSQL is not
running (master prompt rule #22), so SQLite is the default and is fully
sufficient for the closed loop, duplicate detection, and reporting.
"""
from __future__ import annotations

import json
import os
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_PATH = Path(__file__).with_name("schema.sql")
DEFAULT_SQLITE_PATH = Path("artifacts/conform_memory.db")


class Repository:
    def __init__(self, db_path: Path | None = None):
        url = os.environ.get("CONFORM_DB_URL", "")
        if url.startswith("postgres"):
            raise NotImplementedError(
                "PostgreSQL backend is declared in db/schema.sql but the psycopg "
                "driver path is not exercised in this build; unset CONFORM_DB_URL "
                "to use the SQLite backend."
            )
        self.path = Path(db_path or DEFAULT_SQLITE_PATH)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.backend = "sqlite"
        self._init_schema()

    def _init_schema(self):
        self.conn.executescript(SCHEMA_PATH.read_text())
        self.conn.commit()

    # ---- experiments -------------------------------------------------
    def create_experiment(self, experiment_id: str, target: str, mode: str, config: dict):
        self.conn.execute(
            "INSERT OR REPLACE INTO experiment (experiment_id, created_at, target, mode, config_json, status)"
            " VALUES (?, ?, ?, ?, ?, 'running')",
            (experiment_id, datetime.now(timezone.utc).isoformat(), target, mode, json.dumps(config, default=str)),
        )
        self.conn.commit()

    def finish_experiment(self, experiment_id: str, manifest: dict, status: str = "completed"):
        self.conn.execute(
            "UPDATE experiment SET manifest_json = ?, status = ? WHERE experiment_id = ?",
            (json.dumps(manifest, default=str), status, experiment_id),
        )
        self.conn.commit()

    # ---- agent steps -------------------------------------------------
    def log_step(self, experiment_id: str, iteration: int, action: str, input_hash: str,
                 tool: str | None = None, params: dict | None = None, metrics: dict | None = None,
                 artifacts: list | None = None, interpretation: str | None = None,
                 failure: str | None = None, next_action: str | None = None,
                 runtime_seconds: float | None = None):
        self.conn.execute(
            "INSERT INTO agent_step (experiment_id, iteration, timestamp, action, input_hash, tool,"
            " params_json, metrics_json, output_artifacts_json, scientific_interpretation, failure,"
            " next_action, runtime_seconds) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (experiment_id, iteration, datetime.now(timezone.utc).isoformat(), action, input_hash, tool,
             json.dumps(params or {}, default=str), json.dumps(metrics or {}, default=str),
             json.dumps(artifacts or [], default=str), interpretation, failure, next_action, runtime_seconds),
        )
        self.conn.commit()

    def has_completed(self, experiment_id: str, input_hash: str) -> bool:
        """Duplicate-experiment detection: was this exact action already run successfully?"""
        row = self.conn.execute(
            "SELECT 1 FROM agent_step WHERE experiment_id = ? AND input_hash = ? AND failure IS NULL LIMIT 1",
            (experiment_id, input_hash),
        ).fetchone()
        return row is not None

    def steps(self, experiment_id: str) -> list[dict]:
        rows = self.conn.execute(
            "SELECT * FROM agent_step WHERE experiment_id = ? ORDER BY id", (experiment_id,)
        ).fetchall()
        return [dict(r) for r in rows]

    # ---- domain records ----------------------------------------------
    def save_pockets(self, experiment_id: str, pockets: list):
        self.conn.executemany(
            "INSERT INTO pocket_candidate (experiment_id, state_pdb_id, pocket_index, volume,"
            " druggability, residues_json, ground_truth_overlap) VALUES (?,?,?,?,?,?,?)",
            [(experiment_id, p.state_pdb_id, p.pocket_index, p.volume, p.druggability,
              json.dumps(p.residues), p.ground_truth_overlap) for p in pockets],
        )
        self.conn.commit()

    def save_docking(self, experiment_id: str, records: list[dict]):
        self.conn.executemany(
            "INSERT INTO docking_result (experiment_id, ligand_name, smiles, receptor_pdb_id,"
            " pocket_index, engine, best_affinity_kcal, poses_json, discovery_score, origin)"
            " VALUES (?,?,?,?,?,?,?,?,?,?)",
            [(experiment_id, r["ligand_name"], r.get("smiles"), r.get("receptor_pdb_id"),
              r.get("pocket_index"), r.get("engine", "vina"), r.get("best_affinity_kcal"),
              json.dumps(r.get("poses", [])), r.get("discovery_score"), r.get("origin", "library"))
             for r in records],
        )
        self.conn.commit()

    def docking_results(self, experiment_id: str) -> list[dict]:
        rows = self.conn.execute(
            "SELECT * FROM docking_result WHERE experiment_id = ? ORDER BY best_affinity_kcal",
            (experiment_id,),
        ).fetchall()
        return [dict(r) for r in rows]

    def list_experiments(self) -> list[dict]:
        rows = self.conn.execute("SELECT * FROM experiment ORDER BY created_at DESC").fetchall()
        return [dict(r) for r in rows]

    def close(self):
        self.conn.close()
