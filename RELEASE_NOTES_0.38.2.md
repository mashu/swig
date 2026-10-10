# Swig 0.38.2

## Sensitive-D

- Adds the opt-in AER-R-only `sensitive_d` calling profile in Swig Web, saved/exported pipeline configuration, the general CLI, and `swig-cli --vdj`.
- Sensitive-D is R-optimized with one deliberate change: the consecutive exact D-candidate floor is 4 nt instead of 5 nt. The D-presence cost, evidence-conditioned relaxation, ambiguity reporting, V/J scores, and every other threshold are unchanged.
- On the supplied 20k lower-SHM and 20k IgG-like simulations it adds 8 and 30 truth-D-containing reported sets while adding 5 and 6 no-D false calls, respectively.
- The supplied mutation-interrupted 14-nt D regression is now a permanent test: R-optimized omits it, while Sensitive-D recovers the expected `IGHD5-32` sequence class.

## V/J boundary correction

- AER-R now keeps independent V/J alignments immutable while evaluating D candidates. Temporary clip-V/retain-D and retain-V/clip-D allocations are compared as complete partitions; only the winning partition is committed.
- A symmetric allocation must beat both its matching allocation and the global partition strictly. A narrowly validated one-base V tie rescue is permitted only for an exact V continuation with at least ten unaligned 3′ V-reference bases. There is no automatic one-base J tie preference because it worsened aggregate J-start error.
- Both V end and J start can materialize a five-base, score-zero terminal block in a high-SHM context, but only inside the interval bounded by the actual selected segments. Rejected D hypotheses cannot change either endpoint.
- On the IgG-like simulation, R-optimized V-end MAE improves from 0.7987 to 0.7866 nt and J-start MAE from 0.8351 to 0.8205 nt. On lower-SHM, V-end MAE improves from 0.4842 to 0.4836 nt and the complete J endpoint distribution is unchanged.

## Interfaces, compatibility, and speed

- Browser and CLI version strings are 0.38.2. Both surfaces show `Sensitive-D`, enforce its AER-R requirement, retain it in saved/exported configuration, and reuse it for post-analysis query inference.
- Profile identifier 3 is added to the WASM ABI; identifiers 0–2 and all default behavior remain unchanged.
- The default remains RIAT-MP with truth-optimized calling. Ordinary AER, standard SwiftIG, legacy AER-R profiles, and R-optimized selection remain opt-in/unchanged except for the shared AER-R endpoint correction.
- Full one-worker timing showed a 1–8% CPU-time increase versus v0.37.6 AER-R (4–13% versus v0.37.6 plain AER), depending on dataset/profile. See [`BENCHMARK_SENSITIVE_D_0.38.2.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_SENSITIVE_D_0.38.2.md) for accuracy, trimming, fixture, and timing details.
