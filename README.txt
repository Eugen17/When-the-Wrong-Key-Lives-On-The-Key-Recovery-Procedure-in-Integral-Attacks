ARKST — alpha-Admissible Key-Difference Finder (README, short)
==============================================================

Contents
--------
• alpha_admissible_{cipher_name}.c  — C source to search α-admissible subspaces
• alpha_admissible_{cipher_name}    — compiled binary (produced from the C file)
• elimination_GIFT.py               — post-process GIFT output (intersection with actual key-add positions)
• echelon_form.sage                 — convert a basis to row-echelon form

Purpose
-------
Find a basis of α-admissible key differences for PRESENT, GIFT, and RECTANGLE—i.e., key
directions over the last r rounds that do not change whether a chosen integral/divisibility
property (mask α) holds. These form the key subspace you can exclude from key guessing.

Build
-----
# Linux/macOS (clang or gcc)
cc -O3 -march=native -DNDEBUG alpha_admissible_{cipher_name}.c -o alpha_admissible_{cipher_name}

Help
----
./alpha_admissible_present -h

Usage
-----
./alpha_admissible_{cipher_name} [-r ROUNDS] [-Z SAMPLES] [-a ALPHA_HEX] [-s SEED] [--no-zero] [-h]

Defaults
--------
• PRESENT:   r=3, Z=50, α=0x1
• GIFT:      r=3, Z=50, α=0x1
• RECTANGLE: r=3, Z=50, α=0x100

Example
-------
./alpha_admissible_present -r 3 -Z 50 -a 0x1 -s 42

Runtime
-------
1–2 rounds: instant; 3 rounds: typically ≤ 10 minutes (depends on CPU and Z).

Output
------
A basis printed as tuples per round, e.g.:
(0x…, 0x…, 0x…)
Non-unit vectors indicate non-trivial cross-round compensations.

Notes
-----
• Z is the sample size for affine propagation; 50 is usually enough. Use a fixed -s (seed) for reproducibility.
• GIFT: the round key is added only to half the state. After running the finder, intersect the resulting
  subspace with the span of unit vectors at the bit positions where the round key is actually added.
  Use elimination_GIFT.py (paste the C output into `input_str`) to perform this intersection.
• Row-echelon form: use echelon_form.sage by putting the basis (from the C tool or elimination_GIFT.py)
  into the DATA string variable; it prints an echelonized basis.
• Number of key bits to guess =
    (#rounds × size_of_round_key) − (dimension of the resulting α-admissible subspace),
  where the dimension is the count of basis vectors after any post-processing/intersection.
• The tool naturally detects non-trivial key relations: a wrong bit can be compensated by another
  (even in a later round); such dependencies appear as non-unit vectors in the echelonized basis.
