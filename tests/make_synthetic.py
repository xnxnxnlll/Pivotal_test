"""Generate a tiny synthetic reference + reads with KNOWN coverage so the
alignment and coverage modules can be unit-tested without external data.
"""
from __future__ import annotations
import random


def random_genome(length: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(seq: str, err: float, seed: int = 1) -> str:
    """Apply substitutions at rate `err` to mimic ONT mismatch noise."""
    rng = random.Random(seed)
    out = []
    for b in seq:
        if rng.random() < err:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        else:
            out.append(b)
    return "".join(out)


def write_fasta(path: str, records: list[tuple[str, str]]):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def make_case(refdir: str, readpath: str):
    """One 3000 bp reference 'virusX', one 1500 bp 'virusY' (no reads map to it).
    Reads tile positions 0-1000 (depth 2 over 0-500) of virusX with 5% error.
    Returns expected coverage facts for assertions.
    """
    gx = random_genome(3000, seed=42)
    gy = random_genome(1500, seed=7)
    write_fasta(f"{refdir}/panel.fasta", [("virusX", gx), ("virusY", gy)])

    reads = []
    # two reads covering 0-500 (depth 2), one read covering 500-1000 (depth 1)
    reads.append(("r1", mutate(gx[0:500], 0.05, seed=11)))
    reads.append(("r2", mutate(gx[0:500], 0.05, seed=12)))
    reads.append(("r3", mutate(gx[500:1000], 0.05, seed=13)))
    write_fasta(readpath, reads)

    return {
        "ref_lengths": {"virusX": 3000, "virusY": 1500},
        "virusX_covered_bp": 1000,      # positions 0-1000
        "virusX_breadth": 1000 / 3000,
        "virusX_max_depth": 2,
        "virusY_covered_bp": 0,
    }


if __name__ == "__main__":
    import tempfile, os
    d = tempfile.mkdtemp()
    facts = make_case(d, os.path.join(d, "reads.fasta"))
    print("wrote synthetic data to", d)
    print(facts)
