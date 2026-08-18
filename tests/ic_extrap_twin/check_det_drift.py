#!/usr/bin/env python3
"""
The deterministic IC fill exists TWICE — once in this tree, once in
port_kokkos_int — and the twin gate requires the two to agree bit for bit. That
is not a property of the algorithm; it is a property of the tolerance, the
iteration caps and the summation order all being identical. So the two copies
are meant to be verbatim, and this asserts it in a second rather than at the
end of a job on a compute node.

If you deliberately change the fill, change it in BOTH trees in the same commit
and re-run the twin gate. If the copies ever have to differ, this check is the
place to record why.

Exit 0 when the two blocks are identical, 1 when they are not.
"""
import difflib
import sys

C_SRC  = "/home/a/a270088/port2/fesom2_port_zstar/src/fesom_phc.c"
KK_SRC = "/home/a/a270088/port_kokkos_int/src/fesom_phc.cpp"
BEGIN  = "/*--- extrap_nod3D_det"
END    = "/*--- Top-level entry"


def block(path):
    s = open(path).read()
    try:
        i, j = s.index(BEGIN), s.index(END)
    except ValueError:
        raise SystemExit(f"{path}: could not find the det block markers")
    if j < i:
        raise SystemExit(f"{path}: markers out of order")
    return s[i:j].rstrip().split("\n")


def main():
    c, kk = block(C_SRC), block(KK_SRC)
    if c == kk:
        print(f"PASS: the two deterministic-fill copies are identical ({len(c)} lines)")
        return 0
    print(f"FAIL: the deterministic fill has drifted between the trees "
          f"(C {len(c)} lines, Kokkos {len(kk)})")
    for line in list(difflib.unified_diff(c, kk, "C", "Kokkos", lineterm=""))[:60]:
        print(line)
    return 1


if __name__ == "__main__":
    sys.exit(main())
