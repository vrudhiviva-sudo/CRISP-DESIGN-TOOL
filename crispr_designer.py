"""Core sequence-analysis functions for the CRISPR gRNA Designer.

This teaching implementation scans both strands for SpCas9 NGG PAMs.
Coordinates returned to the user are 1-based, inclusive, on the input sequence.
"""
import re
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def validate_dna(sequence: str) -> str:
    """Remove whitespace and FASTA header lines; validate an A/C/G/T sequence."""
    if not isinstance(sequence, str):
        raise ValueError("Sequence must be text.")
    lines = [line.strip() for line in sequence.splitlines() if line.strip()]
    lines = [line for line in lines if not line.startswith(">")]
    cleaned = re.sub(r"\s+", "", "".join(lines)).upper()
    if not cleaned:
        raise ValueError("Please provide a DNA sequence.")
    invalid = sorted(set(cleaned) - set("ACGT"))
    if invalid:
        raise ValueError("Invalid DNA character(s): " + ", ".join(invalid) + ". Use A, C, G and T only.")
    return cleaned


def parse_fasta(text: str) -> str:
    """Parse the first FASTA record using Biopython."""
    from io import StringIO
    from Bio import SeqIO
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    if not records:
        raise ValueError("No FASTA record found. Check that the file has a header beginning with '>'.")
    return validate_dna(str(records[0].seq))


def _has_poly_t(guide: str) -> bool:
    return "TTTT" in guide


def analyze_sequence(sequence: str, guide_length: int = 20, min_gc: float = 40,
                     max_gc: float = 60, exclude_poly_t: bool = True,
                     pam_pattern: str = "NGG"):
    """Return candidate guide records and sequence statistics."""
    seq = validate_dna(sequence)
    if pam_pattern != "NGG":
        raise ValueError("This prototype currently supports the SpCas9 NGG PAM only.")
    if not (0 <= min_gc <= max_gc <= 100):
        raise ValueError("GC limits must satisfy 0 ≤ minimum ≤ maximum ≤ 100.")
    if guide_length < 1:
        raise ValueError("Guide length must be positive.")

    n = len(seq)
    candidates = []
    pam_sites = 0

    # Plus strand: 5' - [guide] [NGG PAM] - 3'
    for pam_start in range(guide_length, n - 2):
        pam = seq[pam_start:pam_start + 3]
        if pam[1:] == "GG":
            pam_sites += 1
            guide = seq[pam_start - guide_length:pam_start]
            gc = gc_fraction(guide) * 100
            if min_gc <= gc <= max_gc and not (exclude_poly_t and _has_poly_t(guide)):
                candidates.append({
                    "Guide (5'→3')": guide,
                    "PAM": pam,
                    "Strand": "+",
                    "Start": pam_start - guide_length + 1,
                    "End": pam_start,
                    "GC %": round(gc, 1),
                    "GC distance from 50": round(abs(gc - 50), 1),
                })

    # Minus strand: reverse-complement sequence and scan as if it were plus.
    # Convert reverse-complement coordinates back to the submitted sequence.
    rc = str(Seq(seq).reverse_complement())
    for pam_start in range(guide_length, n - 2):
        pam = rc[pam_start:pam_start + 3]
        if pam[1:] == "GG":
            pam_sites += 1
            guide = rc[pam_start - guide_length:pam_start]
            gc = gc_fraction(guide) * 100
            if min_gc <= gc <= max_gc and not (exclude_poly_t and _has_poly_t(guide)):
                start = n - pam_start + 1
                end = n - pam_start + guide_length
                candidates.append({
                    "Guide (5'→3')": guide,
                    "PAM": pam,
                    "Strand": "-",
                    "Start": start,
                    "End": end,
                    "GC %": round(gc, 1),
                    "GC distance from 50": round(abs(gc - 50), 1),
                })

    candidates.sort(key=lambda x: (x["GC distance from 50"], x["Start"], x["Strand"]))
    stats = {
        "length": n,
        "sequence_gc": gc_fraction(seq) * 100,
        "pam_sites": pam_sites,
    }
    return candidates, stats
