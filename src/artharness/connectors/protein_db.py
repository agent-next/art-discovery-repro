"""Metagenomic protein database connector — local SQLite subset + MMseqs2 wrapper.

Paper: sessions queried a metagenomic protein database; the census ran over
~1.9B protein clusters, and homology searches ran either against cluster
representatives or a 365M-sequence reference set (UniProt+nr) clustered at
90% identity.

The SQLite schema is a NOT-IN-PAPER local reproduction of the columns the
pipeline needs: cluster membership, contig coordinates, biosample provenance
and Pfam annotations. Genomic coordinates are folded into ``proteins`` (one
protein = one locus on one contig), so a separate coordinates table would only
duplicate the 1:1 key.

:class:`MMseqs2Search` builds ``mmseqs`` command lines only — it never executes
unless ``execute=True`` is passed.
"""

from __future__ import annotations

import shlex
import sqlite3
import subprocess
from collections.abc import Iterable
from pathlib import Path

SCHEMA_SQL = """\
CREATE TABLE IF NOT EXISTS biosamples (
    biosample_id TEXT PRIMARY KEY,
    source       TEXT NOT NULL,   -- e.g. "SRA", "MGnify", "IMG/M" (NOT-IN-PAPER labels)
    metadata     TEXT             -- free-form JSON blob
);
CREATE TABLE IF NOT EXISTS contigs (
    contig_id    TEXT PRIMARY KEY,
    source_db    TEXT NOT NULL,
    biosample_id TEXT REFERENCES biosamples(biosample_id),
    length       INTEGER
);
CREATE TABLE IF NOT EXISTS clusters (
    cluster_id               TEXT PRIMARY KEY,
    representative_protein_id TEXT NOT NULL,
    size                     INTEGER NOT NULL DEFAULT 1
);
CREATE TABLE IF NOT EXISTS proteins (
    protein_id TEXT PRIMARY KEY,
    cluster_id TEXT REFERENCES clusters(cluster_id),
    contig_id  TEXT REFERENCES contigs(contig_id),
    start      INTEGER,
    "end"      INTEGER,
    strand     TEXT CHECK (strand IN ('+', '-') OR strand IS NULL),
    length     INTEGER,
    sequence   TEXT               -- optional: subset DBs may store ids only
);
CREATE TABLE IF NOT EXISTS pfam (
    protein_id TEXT NOT NULL REFERENCES proteins(protein_id),
    pfam_id    TEXT NOT NULL,
    bitscore   REAL,
    coverage   REAL,
    PRIMARY KEY (protein_id, pfam_id)
);
CREATE INDEX IF NOT EXISTS idx_proteins_contig ON proteins(contig_id, start);
CREATE INDEX IF NOT EXISTS idx_proteins_cluster ON proteins(cluster_id);
CREATE INDEX IF NOT EXISTS idx_pfam_pfam ON pfam(pfam_id);
"""


class ProteinDB:
    """SQLite-backed subset of the metagenomic protein database.

    ``path=":memory:"`` gives a transient DB for tests; otherwise the file is
    created if missing. Insert helpers commit per call; the ``load_*`` loaders
    use ``executemany`` and commit once.
    """

    def __init__(self, path: str | Path):
        self.path = str(path)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA foreign_keys = ON")
        self.conn.executescript(SCHEMA_SQL)

    def close(self) -> None:
        self.conn.close()

    # -- loaders -----------------------------------------------------------

    def insert_biosample(self, biosample_id: str, source: str,
                         metadata: str | None = None) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO biosamples (biosample_id, source, metadata)"
            " VALUES (?, ?, ?)",
            (biosample_id, source, metadata),
        )
        self.conn.commit()

    def insert_contig(self, contig_id: str, source_db: str,
                      biosample_id: str | None = None,
                      length: int | None = None) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO contigs (contig_id, source_db, biosample_id, length)"
            " VALUES (?, ?, ?, ?)",
            (contig_id, source_db, biosample_id, length),
        )
        self.conn.commit()

    def insert_cluster(self, cluster_id: str, representative_protein_id: str,
                       size: int = 1) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO clusters"
            " (cluster_id, representative_protein_id, size) VALUES (?, ?, ?)",
            (cluster_id, representative_protein_id, size),
        )
        self.conn.commit()

    def insert_protein(self, protein_id: str, cluster_id: str | None = None,
                       contig_id: str | None = None, start: int | None = None,
                       end: int | None = None, strand: str | None = None,
                       length: int | None = None, sequence: str | None = None) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO proteins"
            " (protein_id, cluster_id, contig_id, start, \"end\", strand, length, sequence)"
            " VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (protein_id, cluster_id, contig_id, start, end, strand, length, sequence),
        )
        self.conn.commit()

    def insert_pfam(self, protein_id: str, pfam_id: str,
                    bitscore: float | None = None,
                    coverage: float | None = None) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO pfam (protein_id, pfam_id, bitscore, coverage)"
            " VALUES (?, ?, ?, ?)",
            (protein_id, pfam_id, bitscore, coverage),
        )
        self.conn.commit()

    def load_proteins(self, rows: Iterable[dict]) -> None:
        """Bulk-load protein dicts (keys matching :meth:`insert_protein`)."""
        self.conn.executemany(
            "INSERT OR REPLACE INTO proteins"
            " (protein_id, cluster_id, contig_id, start, \"end\", strand, length, sequence)"
            " VALUES (:protein_id, :cluster_id, :contig_id, :start, :end,"
            " :strand, :length, :sequence)",
            rows,
        )
        self.conn.commit()

    # -- queries -----------------------------------------------------------

    def _rows(self, sql: str, params: tuple) -> list[dict]:
        return [dict(r) for r in self.conn.execute(sql, params)]

    def neighbors(self, protein_id: str, flank_bp: int) -> list[dict]:
        """Proteins on the same contig overlapping ``[start-flank, end+flank]``.

        NOT-IN-PAPER: window semantics (interval overlap, target excluded).
        """
        row = self.conn.execute(
            'SELECT contig_id, start, "end" FROM proteins WHERE protein_id = ?',
            (protein_id,),
        ).fetchone()
        if row is None:
            raise KeyError(f"unknown protein_id: {protein_id}")
        lo, hi = row["start"] - flank_bp, row["end"] + flank_bp
        return self._rows(
            'SELECT * FROM proteins WHERE contig_id = ? AND protein_id != ?'
            ' AND start <= ? AND "end" >= ? ORDER BY start',
            (row["contig_id"], protein_id, hi, lo),
        )

    def proteins_beside(self, class_: str, flank_bp: int = 10_000) -> list[dict]:
        """Locus neighborhoods: proteins within ``flank_bp`` of any protein
        annotated with Pfam id ``class_`` (anchors included).

        flank default per paper Methods 'Sampling of RT neighborhoods' p.29-30:
        "10 kb of flanking sequence on each side". NOT-IN-PAPER: only the
        semantics — interval overlap and anchoring on Pfam-annotated proteins
        rather than RT-cluster loci.
        """
        return self._rows(
            "SELECT DISTINCT p.* FROM proteins p"
            " JOIN proteins anchor ON anchor.contig_id = p.contig_id"
            " JOIN pfam f ON f.protein_id = anchor.protein_id"
            ' WHERE f.pfam_id = ? AND p.start <= anchor."end" + ?'
            ' AND p."end" >= anchor.start - ?'
            " ORDER BY p.contig_id, p.start",
            (class_, flank_bp, flank_bp),
        )

    def cluster_members(self, cluster_id: str) -> list[dict]:
        """All proteins assigned to ``cluster_id``."""
        return self._rows(
            "SELECT * FROM proteins WHERE cluster_id = ? ORDER BY protein_id",
            (cluster_id,),
        )


class MMseqs2Search:
    """Command-line builder for MMseqs2 homology searches (offline-first).

    ``mode`` selects the target database:

    - ``"representatives"`` — the cluster-representative subset DB.
    - ``"reference365m"``  — the 365M-sequence UniProt+nr reference clustered
      at 90% identity.  # paper: 365M reference set (UniProt+nr) at 90% identity

    ``db_paths`` maps mode -> target DB path; :meth:`build_command` resolves it
    when no explicit ``target_db`` is given.
    """

    MODES = ("representatives", "reference365m")

    def __init__(self, mode: str = "representatives",
                 db_paths: dict[str, str | Path] | None = None,
                 mmseqs_bin: str = "mmseqs",
                 min_seq_id: float | None = None):
        if mode not in self.MODES:
            raise ValueError(f"mode must be one of {self.MODES}, got {mode!r}")
        self.mode = mode
        self.db_paths = {k: str(v) for k, v in (db_paths or {}).items()}
        self.mmseqs_bin = mmseqs_bin
        self.min_seq_id = min_seq_id

    def resolve_target(self) -> str:
        """Target DB path for the configured mode."""
        try:
            return self.db_paths[self.mode]
        except KeyError:
            raise KeyError(
                f"no db_paths entry for mode {self.mode!r}; pass target_db explicitly"
            ) from None

    def build_command(self, query_fasta: str | Path,
                      target_db: str | Path | None = None,
                      result_db: str | Path = "resultDB",
                      tmp_dir: str | Path = "tmp",
                      easy: bool = True,
                      extra: Iterable[str] = ()) -> list[str]:
        """Build an ``mmseqs easy-search`` (default) or ``mmseqs search`` line.

        NOT-IN-PAPER: binary name, temp dir and output layout are local choices.
        The paper's only quoted SEARCH threshold is E <= 0.001; its "90%" describes
        how the 365M reference SET was clustered, not a search filter. The old
        default --min-seq-id 0.9 silently discarded every hit below 90% identity —
        i.e. exactly the remote homologs these searches exist to find (S1 finding
        C9). min_seq_id is now opt-in (None = omit the flag; MMseqs2's own
        defaults apply; filter by e-value downstream).
        """
        target = str(target_db) if target_db is not None else self.resolve_target()
        sub = "easy-search" if easy else "search"
        cmd = [self.mmseqs_bin, sub, str(query_fasta), target,
               str(result_db), str(tmp_dir)]
        if self.min_seq_id is not None:
            cmd += ["--min-seq-id", str(self.min_seq_id)]
        cmd.extend(extra)
        return cmd

    def run(self, *args, execute: bool = False,
            **kwargs) -> subprocess.CompletedProcess:
        """Build the command and — only with ``execute=True`` — run it."""
        cmd = self.build_command(*args, **kwargs)
        if not execute:
            raise RuntimeError(
                "MMseqs2Search is offline-first; pass execute=True to run: "
                + shlex.join(cmd)
            )
        return subprocess.run(cmd, check=True)
