# BlackRoad Marine Biology Lab

[![CI](https://github.com/BlackRoad-OS/blackroad-marine-biology-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/BlackRoad-OS/blackroad-marine-biology-lab/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-proprietary-red.svg)](LICENSE)
[![BlackRoad OS](https://img.shields.io/badge/BlackRoad-OS-black.svg)](https://blackroad.io)

> Bioinformatics: DNA analysis, sequence alignment (NW/SW), BLAST, field sampling

Part of the **BlackRoad OS** health & science platform — production-grade implementations with SQLite persistence, pytest coverage, and CI/CD.

## Features

### Bioinformatics
- `gc_content(seq)` — GC fraction of a DNA sequence
- `reverse_complement(seq)` — reverse complement
- `find_orfs(seq, min_length=100)` — all 6-frame ORFs with protein translation
- `translate_codon(codon)` / `translate_sequence(seq)` — standard genetic code
- `protein_weight(sequence)` — monoisotopic molecular weight (Da)
- `needleman_wunsch(s1, s2)` — global alignment with traceback
- `smith_waterman(s1, s2)` — local alignment with identity %
- `fasta_parser(text)` — multi-FASTA text → `Sequence` objects
- `blast_mock(query, db_seqs, top_n=5)` — ranked hits with E-value

### Field Data
- `log_sample(station_id, lat, lon, depth, temp, ...)` — SQLite field sample logging
- `station_diversity(station_id)` — Shannon H, evenness, species richness

## Quick Start

```bash
# Bioinformatics
python src/marine.py gc --seq ATGCATGCATGC
python src/marine.py revcomp --seq ATGCATGC
python src/marine.py orfs --seq ATGAAAGGGAAATAA --minlen 9
python src/marine.py translate --seq ATGAAATAA
python src/marine.py weight --protein MKVLSPADKTNVK
python src/marine.py align-global --s1 AGTACGCA --s2 TATGC
python src/marine.py align-local  --s1 AGTACGCA --s2 TATGC
python src/marine.py blast --query ATGCATGC --db-file sequences.fasta

# Field data
python src/marine.py log-sample --station S001 --lat 36.5 --lon -121.8 --depth 10 --temp 15.2 --species Kelp "Sea Otter" --count 5
python src/marine.py diversity --station S001
```

## Installation

```bash
# No dependencies required — pure Python stdlib + sqlite3
python src/marine.py --help
```

## Testing

```bash
pip install pytest pytest-cov
pytest tests/ -v --cov=src
```

## Data Storage

All data is stored locally in `~/.blackroad/marine-biology-lab.db` (SQLite). Zero external dependencies.

## License

Proprietary — © BlackRoad OS, Inc. All rights reserved.
