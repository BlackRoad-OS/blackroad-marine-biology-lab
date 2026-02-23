"""Tests for BlackRoad Marine Biology Lab."""
import os, tempfile, pytest
os.environ["HOME"] = tempfile.mkdtemp()

from marine import (
    gc_content, reverse_complement, translate_codon, translate_sequence,
    find_orfs, protein_weight, fasta_parser, needleman_wunsch, smith_waterman,
    blast_mock, log_sample, station_diversity, Sequence,
)

def test_gc_content():
    assert gc_content("ATGC") == 0.5
    assert gc_content("AAAA") == 0.0
    assert gc_content("GGCC") == 1.0

def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AATTCC") == "GGAATT"

def test_translate_codon():
    assert translate_codon("ATG") == "M"
    assert translate_codon("TAA") == "*"
    assert translate_codon("GGG") == "G"

def test_translate_sequence():
    # ATG (M) + AAA (K) + TGA (stop)
    prot = translate_sequence("ATGAAATGA")
    assert prot == "MK"

def test_find_orfs():
    seq = "ATG" + "AAA" * 20 + "TAA"  # 66 nt ORF
    orfs = find_orfs(seq, min_length=10)
    assert len(orfs) >= 1
    assert orfs[0]["protein"].startswith("M")

def test_protein_weight():
    w = protein_weight("MK")
    assert w > 0
    assert w > protein_weight("M")

def test_fasta_parser():
    fasta = ">seq1 desc one\nATGCATGC\n>seq2\nCCGGTTAA\n"
    seqs  = fasta_parser(fasta)
    assert len(seqs) == 2
    assert seqs[0].id == "seq1"
    assert seqs[0].sequence == "ATGCATGC"

def test_needleman_wunsch_identical():
    r = needleman_wunsch("ACGT", "ACGT")
    assert r["identity"] == 100.0

def test_needleman_wunsch_divergent():
    r = needleman_wunsch("AAAA", "TTTT")
    assert r["identity"] < 50

def test_smith_waterman_local():
    r = smith_waterman("AGTACGCA", "TATGC")
    assert r["score"] >= 0
    assert r["algorithm"] == "smith_waterman"

def test_blast_mock_top_hit():
    query = "ATGCATGCATGC"
    db = [
        Sequence("exact", "exact match", "ATGCATGCATGC"),
        Sequence("diff",  "different",   "TTTTTTTTTTTT"),
    ]
    hits = blast_mock(query, db, top_n=2)
    assert hits[0]["subject_id"] == "exact"
    assert hits[0]["score"] >= hits[1]["score"]

def test_log_sample():
    s = log_sample("STN001", 36.5, -121.8, depth_m=10, temp_c=15.2,
                   species=["Kelp", "Sea Otter"], count=3)
    assert s.id > 0
    assert s.station_id == "STN001"

def test_diversity_empty():
    d = station_diversity("NONEXISTENT", days=1)
    assert d["richness"] == 0
