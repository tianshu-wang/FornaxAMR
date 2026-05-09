#!/usr/bin/env python3
"""
Compare two HDF5 dump files from magnetized_sedov runs (e.g. 1-rank vs 7-rank).
Blocks are matched by (level, coordinate) since MPI may reorder them.
All fields are compared bitwise (exact double equality via numpy view).
"""

import sys
import numpy as np
import h5py


CELL_FIELDS = ["density", "eint", "v1", "v2", "v3", "ye", "B1", "B2", "B3"]
FACE_FIELDS = ["Bf1", "Bf2", "Bf3"]
ALL_FIELDS  = CELL_FIELDS + FACE_FIELDS


def load_file(path):
    with h5py.File(path, "r") as f:
        coords = f["coordinate"][:]          # (nblocks, 3)
        levels = f["level"][:]               # (nblocks,)
        data = {name: f[name][:] for name in ALL_FIELDS}
    return coords, levels, data


def make_block_key(coord, level):
    # Round to avoid floating-point jitter in coordinate storage
    return (int(level), round(float(coord[0]), 10),
                         round(float(coord[1]), 10),
                         round(float(coord[2]), 10))


def compare_files(path1, path2):
    coords1, levels1, data1 = load_file(path1)
    coords2, levels2, data2 = load_file(path2)

    n1 = len(levels1)
    n2 = len(levels2)
    print(f"  File A: {n1} blocks")
    print(f"  File B: {n2} blocks")

    if n1 != n2:
        print(f"  ERROR: block count mismatch ({n1} vs {n2})")
        return False

    # Build key -> index map for file 2
    idx2 = {}
    for b in range(n2):
        key = make_block_key(coords2[b], levels2[b])
        if key in idx2:
            print(f"  ERROR: duplicate block key in file B: {key}")
            return False
        idx2[key] = b

    first_diff_step = None   # (field, block_key, b1_idx, b2_idx, cell_pos)
    total_blocks_diff = 0

    for b1 in range(n1):
        key = make_block_key(coords1[b1], levels1[b1])
        b2 = idx2.get(key)
        if b2 is None:
            print(f"  ERROR: block {key} from file A not found in file B")
            return False

        block_has_diff = False
        for field in ALL_FIELDS:
            arr1 = data1[field][b1]
            arr2 = data2[field][b2]
            # Bitwise comparison via uint64 view
            u1 = arr1.view(np.uint64)
            u2 = arr2.view(np.uint64)
            diff_mask = (u1 != u2)
            if diff_mask.any():
                if not block_has_diff:
                    block_has_diff = True
                    total_blocks_diff += 1
                if first_diff_step is None:
                    idx = np.argwhere(diff_mask)[0]
                    first_diff_step = (field, key, b1, b2, tuple(idx))

    return first_diff_step, total_blocks_diff, n1


def main():
    if len(sys.argv) != 3:
        print("Usage: compare_dumps.py <file_a.h5> <file_b.h5>")
        sys.exit(1)

    path1, path2 = sys.argv[1], sys.argv[2]
    print(f"\nComparing:\n  A: {path1}\n  B: {path2}")

    result = compare_files(path1, path2)
    if result is False:
        print("  => COMPARISON FAILED (structure mismatch)")
        return

    first_diff, total_blocks_diff, nblocks = result

    if first_diff is None:
        print(f"  => BITWISE IDENTICAL (all {nblocks} blocks match in all fields)")
    else:
        field, key, b1_idx, b2_idx, cell_idx = first_diff
        print(f"  => NOT identical: {total_blocks_diff}/{nblocks} blocks differ")
        print(f"\n  First difference detected:")
        print(f"    Field     : {field}")
        print(f"    Block key : level={key[0]}, coord=({key[1]:.6g}, {key[2]:.6g}, {key[3]:.6g})")
        print(f"    Block idx : A[{b1_idx}]  B[{b2_idx}]")
        print(f"    Cell index: {cell_idx}")
        # Print the actual values
        with h5py.File(path1, "r") as fa, h5py.File(path2, "r") as fb:
            va = fa[field][b1_idx][tuple(cell_idx)]
            vb = fb[field][b2_idx][tuple(cell_idx)]
        print(f"    Value in A: {va!r}")
        print(f"    Value in B: {vb!r}")
        print(f"    Difference: {abs(va - vb):.6e}  (rel: {abs(va-vb)/max(abs(va),abs(vb),1e-300):.6e})")


if __name__ == "__main__":
    main()
