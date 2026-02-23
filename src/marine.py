#!/usr/bin/env python3
"""
BlackRoad Marine Biology Lab — Bioinformatics & Marine Science Module
DNA/protein analysis, sequence alignment, ecological stats, field data logging.

Usage:
    python marine.py gc --seq ATGCATGCATGC
    python marine.py revcomp --seq ATGCATGC
    python marine.py orfs --seq ATGAAATGA
    python marine.py translate --seq ATGAAATGA
    python marine.py weight --protein MKVL
    python marine.py align-global --s1 AGTACGCA --s2 TATGC
    python marine.py align-local  --s1 AGTACGCA --s2 TATGC
    python marine.py fasta --file sequences.fasta
    python marine.py blast --query ATGCATGC --db-file db.fasta --top 5
    python marine.py log-sample --station S001 --lat 36.5 --lon -121.8 --depth 10 --temp 15.2
    python marine.py diversity --station S001
"""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
from dataclasses import dataclass, asdict, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DB_PATH = Path.home() / ".blackroad" / "marine_biology.db"

# ── Genetic code (standard) ─────────────────────────────────────────────────

CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Monoisotopic residue masses (Da)
AA_MASS: Dict[str, float] = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
}
WATER_MASS = 18.01056


# ── Dataclasses ────────────────────────────────────────────────────────────────

@dataclass
class Sequence:
    id:          str
    description: str
    sequence:    str
    seq_type:    str = "dna"   # dna | rna | protein

    def __len__(self) -> int:
        return len(self.sequence)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class FieldSample:
    id:          int
    station_id:  str
    latitude:    float
    longitude:   float
    depth_m:     float
    temp_c:      float
    salinity:    float
    ph:          float
    do_mgl:      float          # dissolved oxygen mg/L
    turbidity:   float
    species:     str            # JSON list of species observed
    count:       int
    notes:       str
    sampled_at:  str

    def to_dict(self) -> dict:
        d = asdict(self)
        try:
            d["species_list"] = json.loads(self.species)
        except Exception:
            d["species_list"] = []
        return d


# ── Sequence utilities ─────────────────────────────────────────────────────────

def gc_content(seq: str) -> float:
    """Return GC content as a fraction [0, 1]."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(gc / len(seq), 4)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g",
            "N": "N", "n": "n"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def translate_codon(codon: str) -> str:
    """Translate a single RNA/DNA codon to one-letter amino acid."""
    return CODON_TABLE.get(codon.upper().replace("U", "T"), "X")


def translate_sequence(seq: str) -> str:
    """Translate a DNA/RNA sequence to protein (stops at first *)."""
    seq  = seq.upper().replace("U", "T")
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i+3], "X")
        if aa == "*":
            break
        prot.append(aa)
    return "".join(prot)


def find_orfs(seq: str, min_length: int = 100) -> List[dict]:
    """
    Find all open reading frames in all 6 frames.
    Returns ORFs with start, stop, frame, strand, protein.
    """
    seq  = seq.upper().replace("U", "T")
    orfs = []

    def _scan(s: str, strand: str) -> None:
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                codon = s[i:i+3]
                if codon == "ATG":
                    start = i
                    prot  = []
                    j     = i
                    while j < len(s) - 2:
                        c = s[j:j+3]
                        aa = CODON_TABLE.get(c, "X")
                        if aa == "*":
                            orf_len = j + 3 - start
                            if orf_len >= min_length:
                                orfs.append({
                                    "start":   start,
                                    "stop":    j + 3,
                                    "frame":   frame + 1,
                                    "strand":  strand,
                                    "length":  orf_len,
                                    "protein": "".join(prot),
                                })
                            break
                        prot.append(aa)
                        j += 3
                i += 3

    _scan(seq, "+")
    _scan(reverse_complement(seq), "-")
    return sorted(orfs, key=lambda x: -x["length"])


def protein_weight(sequence: str) -> float:
    """Calculate monoisotopic molecular weight (Da) of a protein sequence."""
    sequence = sequence.upper().replace("*", "")
    mass = WATER_MASS
    for aa in sequence:
        mass += AA_MASS.get(aa, 111.1)  # use average if unknown
    return round(mass, 4)


def fasta_parser(text: str) -> List[Sequence]:
    """Parse multi-FASTA text into a list of Sequence objects."""
    seqs  = []
    cur_id   = None
    cur_desc = ""
    cur_seq  = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                seqs.append(Sequence(cur_id, cur_desc, "".join(cur_seq)))
            parts    = line[1:].split(None, 1)
            cur_id   = parts[0]
            cur_desc = parts[1] if len(parts) > 1 else ""
            cur_seq  = []
        else:
            cur_seq.append(line)
    if cur_id is not None:
        seqs.append(Sequence(cur_id, cur_desc, "".join(cur_seq)))
    return seqs


# ── Sequence Alignment ─────────────────────────────────────────────────────────

def _score(a: str, b: str, match: int = 2, mismatch: int = -1) -> int:
    return match if a == b else mismatch


def needleman_wunsch(seq1: str, seq2: str, gap: int = -2) -> dict:
    """Needleman-Wunsch global alignment. Returns aligned sequences and score."""
    m, n = len(seq1), len(seq2)
    dp   = [[0.0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            dp[i][j] = max(
                dp[i-1][j-1] + _score(seq1[i-1], seq2[j-1]),
                dp[i-1][j]   + gap,
                dp[i][j-1]   + gap,
            )

    # Traceback
    a1, a2 = [], []
    i, j   = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + _score(seq1[i-1], seq2[j-1]):
            a1.append(seq1[i-1]); a2.append(seq2[j-1]); i -= 1; j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            a1.append(seq1[i-1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(seq2[j-1]); j -= 1

    al1  = "".join(reversed(a1))
    al2  = "".join(reversed(a2))
    ident = sum(1 for a, b in zip(al1, al2) if a == b and a != "-")
    return {
        "algorithm":  "needleman_wunsch",
        "score":      dp[m][n],
        "aligned1":   al1,
        "aligned2":   al2,
        "identity":   round(ident / max(len(al1), 1) * 100, 2),
        "length":     len(al1),
    }


def smith_waterman(seq1: str, seq2: str, gap: int = -2) -> dict:
    """Smith-Waterman local alignment."""
    m, n    = len(seq1), len(seq2)
    dp      = [[0.0] * (n + 1) for _ in range(m + 1)]
    best    = 0.0
    best_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            val = max(
                0,
                dp[i-1][j-1] + _score(seq1[i-1], seq2[j-1]),
                dp[i-1][j]   + gap,
                dp[i][j-1]   + gap,
            )
            dp[i][j] = val
            if val > best:
                best     = val
                best_pos = (i, j)

    # Traceback from best
    a1, a2 = [], []
    i, j   = best_pos
    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i-1][j-1] + _score(seq1[i-1], seq2[j-1]):
            a1.append(seq1[i-1]); a2.append(seq2[j-1]); i -= 1; j -= 1
        elif dp[i-1][j] >= dp[i][j-1]:
            a1.append(seq1[i-1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(seq2[j-1]); j -= 1

    al1   = "".join(reversed(a1))
    al2   = "".join(reversed(a2))
    ident = sum(1 for a, b in zip(al1, al2) if a == b and a != "-")
    return {
        "algorithm": "smith_waterman",
        "score":     best,
        "aligned1":  al1,
        "aligned2":  al2,
        "identity":  round(ident / max(len(al1), 1) * 100, 2),
        "length":    len(al1),
    }


def blast_mock(query: str, db_seqs: List[Sequence], top_n: int = 5) -> List[dict]:
    """
    Simplified BLAST-like local alignment of query vs database sequences.
    Uses Smith-Waterman internally; returns top_n hits sorted by score.
    """
    results = []
    for seq in db_seqs:
        aln = smith_waterman(query.upper(), seq.sequence.upper())
        results.append({
            "subject_id":   seq.id,
            "description":  seq.description,
            "score":        aln["score"],
            "identity_pct": aln["identity"],
            "aligned_len":  aln["length"],
            "query_aln":    aln["aligned1"],
            "subject_aln":  aln["aligned2"],
            "e_value":      round(math.exp(-aln["score"] / 10), 6),
        })
    results.sort(key=lambda x: -x["score"])
    return results[:top_n]


# ── Database ───────────────────────────────────────────────────────────────────

def get_conn() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    _init_db(conn)
    return conn


def _init_db(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS field_samples (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            station_id  TEXT NOT NULL,
            latitude    REAL NOT NULL DEFAULT 0,
            longitude   REAL NOT NULL DEFAULT 0,
            depth_m     REAL NOT NULL DEFAULT 0,
            temp_c      REAL NOT NULL DEFAULT 0,
            salinity    REAL NOT NULL DEFAULT 35,
            ph          REAL NOT NULL DEFAULT 8.1,
            do_mgl      REAL NOT NULL DEFAULT 6,
            turbidity   REAL NOT NULL DEFAULT 0,
            species     TEXT NOT NULL DEFAULT '[]',
            count       INTEGER NOT NULL DEFAULT 0,
            notes       TEXT NOT NULL DEFAULT '',
            sampled_at  TEXT NOT NULL DEFAULT (datetime('now'))
        );
        CREATE INDEX IF NOT EXISTS idx_fs_station ON field_samples(station_id);
        CREATE INDEX IF NOT EXISTS idx_fs_ts      ON field_samples(sampled_at);

        CREATE TABLE IF NOT EXISTS sequence_db (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            seq_id      TEXT NOT NULL UNIQUE,
            description TEXT NOT NULL DEFAULT '',
            sequence    TEXT NOT NULL,
            seq_type    TEXT NOT NULL DEFAULT 'dna',
            gc_content  REAL,
            length      INTEGER,
            created_at  TEXT NOT NULL DEFAULT (datetime('now'))
        );
    """)
    conn.commit()


def log_sample(
    station_id: str,
    latitude: float,
    longitude: float,
    depth_m: float,
    temp_c: float,
    salinity: float = 35.0,
    ph: float = 8.1,
    do_mgl: float = 6.0,
    turbidity: float = 0.0,
    species: Optional[List[str]] = None,
    count: int = 0,
    notes: str = "",
) -> FieldSample:
    sp_json = json.dumps(species or [])
    ts      = datetime.now().isoformat()
    with get_conn() as conn:
        cur = conn.execute(
            "INSERT INTO field_samples(station_id,latitude,longitude,depth_m,temp_c,"
            "salinity,ph,do_mgl,turbidity,species,count,notes,sampled_at)"
            " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (station_id, latitude, longitude, depth_m, temp_c,
             salinity, ph, do_mgl, turbidity, sp_json, count, notes, ts),
        )
        conn.commit()
        row_id = cur.lastrowid
    return FieldSample(row_id, station_id, latitude, longitude, depth_m, temp_c,
                       salinity, ph, do_mgl, turbidity, sp_json, count, notes, ts)


def station_diversity(station_id: str, days: int = 90) -> dict:
    """Shannon diversity index and species richness for a station."""
    from collections import Counter
    import math as _math

    since = (datetime.now() - __import__("datetime").timedelta(days=days)).isoformat()
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT species, count FROM field_samples WHERE station_id=? AND sampled_at>=?",
            (station_id, since),
        ).fetchall()

    counter: Counter = Counter()
    for row in rows:
        try:
            sp_list = json.loads(row["species"])
            cnt     = row["count"] or 1
            for sp in sp_list:
                counter[sp] += cnt
        except Exception:
            pass

    total = sum(counter.values())
    if total == 0:
        return {"station_id": station_id, "richness": 0, "shannon_h": 0.0, "species": {}}

    shannon = -sum(
        (c / total) * _math.log(c / total)
        for c in counter.values() if c > 0
    )
    return {
        "station_id": station_id,
        "days":       days,
        "richness":   len(counter),
        "total_obs":  total,
        "shannon_h":  round(shannon, 4),
        "evenness":   round(shannon / _math.log(len(counter)), 4) if len(counter) > 1 else 1.0,
        "species":    dict(counter.most_common()),
    }


# ── CLI ────────────────────────────────────────────────────────────────────────

def _print(obj):
    print(json.dumps(obj, indent=2, default=str))


def main():
    parser = argparse.ArgumentParser(description="BlackRoad Marine Biology Lab")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("gc", help="GC content of a DNA sequence")
    p.add_argument("--seq", required=True)

    p = sub.add_parser("revcomp", help="Reverse complement")
    p.add_argument("--seq", required=True)

    p = sub.add_parser("orfs", help="Find open reading frames")
    p.add_argument("--seq",    required=True)
    p.add_argument("--minlen", type=int, default=30, dest="min_length")

    p = sub.add_parser("translate", help="Translate DNA → protein")
    p.add_argument("--seq", required=True)

    p = sub.add_parser("codon", help="Translate a single codon")
    p.add_argument("--codon", required=True)

    p = sub.add_parser("weight", help="Protein molecular weight")
    p.add_argument("--protein", required=True)

    p = sub.add_parser("align-global", help="Needleman-Wunsch global alignment")
    p.add_argument("--s1", required=True)
    p.add_argument("--s2", required=True)
    p.add_argument("--gap", type=int, default=-2)

    p = sub.add_parser("align-local", help="Smith-Waterman local alignment")
    p.add_argument("--s1", required=True)
    p.add_argument("--s2", required=True)
    p.add_argument("--gap", type=int, default=-2)

    p = sub.add_parser("fasta", help="Parse a FASTA file")
    p.add_argument("--file", required=True)

    p = sub.add_parser("blast", help="BLAST-like local alignment vs FASTA DB")
    p.add_argument("--query",   required=True)
    p.add_argument("--db-file", required=True, dest="db_file")
    p.add_argument("--top",     type=int, default=5)

    p = sub.add_parser("log-sample", help="Log a field sample")
    p.add_argument("--station",   required=True, dest="station_id")
    p.add_argument("--lat",       type=float, required=True, dest="latitude")
    p.add_argument("--lon",       type=float, required=True, dest="longitude")
    p.add_argument("--depth",     type=float, default=0, dest="depth_m")
    p.add_argument("--temp",      type=float, default=20, dest="temp_c")
    p.add_argument("--salinity",  type=float, default=35)
    p.add_argument("--ph",        type=float, default=8.1)
    p.add_argument("--do",        type=float, default=6.0, dest="do_mgl")
    p.add_argument("--turbidity", type=float, default=0)
    p.add_argument("--species",   nargs="*", default=[])
    p.add_argument("--count",     type=int, default=0)
    p.add_argument("--notes",     default="")

    p = sub.add_parser("diversity", help="Species diversity stats for a station")
    p.add_argument("--station", required=True, dest="station_id")
    p.add_argument("--days",    type=int, default=90)

    args = parser.parse_args()

    if args.cmd == "gc":
        gc = gc_content(args.seq)
        _print({"sequence": args.seq[:40], "gc_content": gc, "at_content": round(1 - gc, 4)})

    elif args.cmd == "revcomp":
        _print({"original": args.seq, "reverse_complement": reverse_complement(args.seq)})

    elif args.cmd == "orfs":
        orfs = find_orfs(args.seq, args.min_length)
        _print({"n_orfs": len(orfs), "orfs": orfs})

    elif args.cmd == "translate":
        _print({"dna": args.seq, "protein": translate_sequence(args.seq)})

    elif args.cmd == "codon":
        _print({"codon": args.codon.upper(), "amino_acid": translate_codon(args.codon)})

    elif args.cmd == "weight":
        _print({"protein": args.protein, "mass_da": protein_weight(args.protein),
                "length": len(args.protein.replace("*",""))})

    elif args.cmd == "align-global":
        _print(needleman_wunsch(args.s1, args.s2, args.gap))

    elif args.cmd == "align-local":
        _print(smith_waterman(args.s1, args.s2, args.gap))

    elif args.cmd == "fasta":
        text = Path(args.file).read_text()
        seqs = fasta_parser(text)
        _print([s.to_dict() for s in seqs])

    elif args.cmd == "blast":
        query   = args.query
        db_text = Path(args.db_file).read_text()
        db_seqs = fasta_parser(db_text)
        hits    = blast_mock(query, db_seqs, args.top)
        _print({"query": query, "n_hits": len(hits), "hits": hits})

    elif args.cmd == "log-sample":
        s = log_sample(
            args.station_id, args.latitude, args.longitude, args.depth_m,
            args.temp_c, args.salinity, args.ph, args.do_mgl,
            args.turbidity, args.species, args.count, args.notes,
        )
        _print({"status": "logged", "id": s.id, "station": s.station_id})

    elif args.cmd == "diversity":
        _print(station_diversity(args.station_id, args.days))


if __name__ == "__main__":
    main()
