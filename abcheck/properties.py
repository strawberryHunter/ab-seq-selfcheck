"""Physicochemical properties — pI, GRAVY, and net charge at a given pH.

All of these are **sequence-based estimates**. Measured pI depends on structure,
glycosylation and buffer composition, so use these values for relative comparison only.

Computed for the whole domain and for each CDR. The pI of an excised fragment such as a
CDR is computed treating it as a standalone peptide — i.e. assuming a free N-terminal amine
and C-terminal carboxyl at its ends. Real CDRs are embedded in the domain and have no free
termini, so read CDR pI only as an indicator of "this region's charge composition" (see DECISIONS.md).
"""

from __future__ import annotations

from dataclasses import dataclass

# --------------------------------------------------------------------------
# pKa — EMBOSS (iep) values
# --------------------------------------------------------------------------

#: Side-chain pKa values, as used in EMBOSS's Epk.dat.
SIDE_CHAIN_PKA: dict[str, float] = {
    "D": 3.9,
    "E": 4.1,
    "C": 8.5,
    "Y": 10.1,
    "H": 6.5,
    "K": 10.8,
    "R": 12.5,
}

#: Positively charged side chains (+1 when protonated)
POSITIVE_RESIDUES = frozenset("KRH")

#: Negatively charged side chains (-1 when deprotonated)
NEGATIVE_RESIDUES = frozenset("DECY")

PKA_N_TERMINUS = 8.6
PKA_C_TERMINUS = 3.6

#: Kyte-Doolittle hydropathy index
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

#: Default pH for reporting net charge, close to physiological conditions.
DEFAULT_PH = 7.4

#: pI search range and convergence criteria
_PI_LOW, _PI_HIGH = 0.0, 14.0
_PI_TOLERANCE = 1e-6
_PI_MAX_ITERATIONS = 200


# --------------------------------------------------------------------------
# Calculations
# --------------------------------------------------------------------------

def net_charge(seq: str, ph: float = DEFAULT_PH, *, include_termini: bool = True) -> float:
    """Net charge at the given pH.

    Sums the average charge of each ionizable group via Henderson-Hasselbalch.
    Positive groups: 1 / (1 + 10^(pH - pKa))
    Negative groups: -1 / (1 + 10^(pKa - pH))
    """
    seq = seq.upper()
    charge = 0.0

    if include_termini and seq:
        charge += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERMINUS))
        charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - ph))

    for residue in seq:
        pka = SIDE_CHAIN_PKA.get(residue)
        if pka is None:
            continue
        if residue in POSITIVE_RESIDUES:
            charge += 1.0 / (1.0 + 10 ** (ph - pka))
        elif residue in NEGATIVE_RESIDUES:
            charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def isoelectric_point(seq: str, *, include_termini: bool = True) -> float | None:
    """The pH at which net charge is zero, found by bisection.

    Returns None if there are no ionizable groups (e.g. 'AAAA' with termini excluded).
    """
    seq = seq.upper()
    if not seq:
        return None

    has_group = include_termini or any(r in SIDE_CHAIN_PKA for r in seq)
    if not has_group:
        return None

    low, high = _PI_LOW, _PI_HIGH
    charge_low = net_charge(seq, low, include_termini=include_termini)
    charge_high = net_charge(seq, high, include_termini=include_termini)

    # Same sign across the whole range means pI is outside it — return the boundary.
    if charge_low * charge_high > 0:
        return low if abs(charge_low) < abs(charge_high) else high

    for _ in range(_PI_MAX_ITERATIONS):
        mid = (low + high) / 2.0
        charge_mid = net_charge(seq, mid, include_termini=include_termini)
        if abs(charge_mid) < _PI_TOLERANCE or (high - low) < _PI_TOLERANCE:
            return mid
        # charge decreases monotonically with pH
        if charge_mid > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def gravy(seq: str) -> float | None:
    """GRAVY — the arithmetic mean of Kyte-Doolittle hydropathy.

    Positive means overall hydrophobic, negative hydrophilic. Characters not in the table are skipped.
    """
    seq = seq.upper()
    # Every Kyte-Doolittle value has one decimal place, so sum them as integer tenths and
    # divide once. Summing the floats directly gave different last bits on different
    # Pythons (sum() became compensated in 3.12), which flipped GRAVY at the 3rd decimal
    # (-0.257 on 3.11, -0.258 on 3.13 for the same sequence). Integer sums are exact on
    # every version and in the browser engine, and int / int is correctly rounded.
    tenths = [round(KYTE_DOOLITTLE[r] * 10) for r in seq if r in KYTE_DOOLITTLE]
    if not tenths:
        return None
    return sum(tenths) / (10 * len(tenths))


# --------------------------------------------------------------------------
# Bundled results
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Properties:
    """Properties of a single region."""

    label: str
    sequence: str
    length: int
    pi: float | None
    gravy: float | None
    net_charge_at_ph: float
    ph: float

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "sequence": self.sequence,
            "length": self.length,
            "pI": None if self.pi is None else round(self.pi, 2),
            "GRAVY": None if self.gravy is None else round(self.gravy, 3),
            "net_charge": round(self.net_charge_at_ph, 2),
            "pH": self.ph,
        }


def compute(label: str, seq: str, ph: float = DEFAULT_PH) -> Properties:
    """Compute properties for a single region."""
    seq = seq.upper()
    return Properties(
        label=label,
        sequence=seq,
        length=len(seq),
        pi=isoelectric_point(seq),
        gravy=gravy(seq),
        net_charge_at_ph=net_charge(seq, ph),
        ph=ph,
    )


def compute_all(sequence: str, cdrs: dict, ph: float = DEFAULT_PH) -> dict[str, Properties]:
    """Compute properties for the whole domain and each CDR.

    cdrs just needs the same shape as abcheck.annotate.Annotation.cdrs
    ({name: object with start/end/seq}). With no CDRs, only the domain is returned.
    """
    result: dict[str, Properties] = {"도메인 전체": compute("도메인 전체", sequence, ph)}
    for name, cdr in sorted(cdrs.items(), key=lambda kv: kv[1].start):
        result[name] = compute(name, cdr.seq, ph)
    return result
