/** Complement of each IUPAC nucleotide code; U pairs as T. */
const IUPAC_COMPLEMENT: Record<string, string> = {
  A: "T", C: "G", G: "C", T: "A", U: "A",
  R: "Y", Y: "R", K: "M", M: "K", B: "V", V: "B", D: "H", H: "D",
  S: "S", W: "W", N: "N",
};

/**
 * AIRR `rev_comp`: true when the row's coordinates and alignments refer to
 * the reverse complement of `sequence`. Accepts T/TRUE/1 in any case.
 */
export function isReverseComplement(value: string | undefined): boolean {
  return /^(?:T|TRUE|1)$/i.test(String(value ?? "").trim());
}

/** Upper-case reverse complement; IUPAC codes are complemented, other symbols become N. */
export function reverseComplement(sequence: string): string {
  return [...sequence.toUpperCase()].reverse().map((base) => IUPAC_COMPLEMENT[base] ?? "N").join("");
}

/** `sequence` as the aligner aligned it: upper-case, reverse-complemented when `rev_comp` is true. */
export function orientedSequence(sequence: string | undefined, revComp: string | undefined): string {
  const raw = sequence ?? "";
  return isReverseComplement(revComp) ? reverseComplement(raw) : raw.toUpperCase();
}
