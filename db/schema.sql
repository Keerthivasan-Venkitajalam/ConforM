-- Scientific memory schema.
-- Written for SQLite (the local default). PostgreSQL notes inline: swap
-- TEXT->JSONB for the json columns and add `CREATE EXTENSION vector;` plus a
-- vector(N) column on conformational_state if pgvector similarity search is
-- wanted. The pipeline is fully usable on SQLite alone.

CREATE TABLE IF NOT EXISTS experiment (
    experiment_id   TEXT PRIMARY KEY,
    created_at      TEXT NOT NULL,
    target          TEXT NOT NULL,
    mode            TEXT NOT NULL,
    config_json     TEXT NOT NULL,
    manifest_json   TEXT,
    status          TEXT NOT NULL DEFAULT 'running'
);

-- One row per agent action actually executed (success or failure).
CREATE TABLE IF NOT EXISTS agent_step (
    id                        INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id             TEXT NOT NULL REFERENCES experiment(experiment_id),
    iteration                 INTEGER NOT NULL,
    timestamp                 TEXT NOT NULL,
    action                    TEXT NOT NULL,
    input_hash                TEXT NOT NULL,
    tool                      TEXT,
    params_json               TEXT,
    metrics_json              TEXT,
    output_artifacts_json     TEXT,
    scientific_interpretation TEXT,
    failure                   TEXT,
    next_action               TEXT,
    runtime_seconds           REAL
);

CREATE INDEX IF NOT EXISTS idx_agent_step_hash ON agent_step(experiment_id, input_hash);

CREATE TABLE IF NOT EXISTS pocket_candidate (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id         TEXT NOT NULL REFERENCES experiment(experiment_id),
    state_pdb_id          TEXT NOT NULL,
    pocket_index          INTEGER NOT NULL,
    volume                REAL,
    druggability          REAL,
    residues_json         TEXT,
    ground_truth_overlap  REAL
);

CREATE TABLE IF NOT EXISTS docking_result (
    id                 INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id      TEXT NOT NULL REFERENCES experiment(experiment_id),
    ligand_name        TEXT NOT NULL,
    smiles             TEXT,
    receptor_pdb_id    TEXT,
    pocket_index       INTEGER,
    engine             TEXT NOT NULL,
    best_affinity_kcal REAL,
    poses_json         TEXT,
    discovery_score    REAL,
    origin             TEXT
);
