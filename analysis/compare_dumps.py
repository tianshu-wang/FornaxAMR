#!/usr/bin/env python3
"""
Compare two HDF5 dump files from magnetized_sedov runs (e.g. 1-rank vs 7-rank).
Blocks are matched by (level, coordinate) since MPI may reorder them.
All grouped dump variables are compared bitwise.
"""

import sys
import numpy as np
import h5py


def decode_variable_name(value):
    if isinstance(value, bytes):
        return value.decode("utf-8").rstrip("\x00")
    return str(value).rstrip("\x00")


def decode_variable_names(attr):
    values = np.asarray(attr)
    if values.shape == ():
        text = decode_variable_name(values.item())
        return [part.strip() for part in text.split(",") if part.strip()]
    return [decode_variable_name(value) for value in values.tolist()]


def build_variable_map(h5):
    variable_map = {}
    variable_order = []
    for dataset_name, obj in h5.items():
        if not isinstance(obj, h5py.Dataset) or "variable_names" not in obj.attrs:
            continue
        names = decode_variable_names(obj.attrs["variable_names"])
        if obj.ndim < 2 or obj.shape[1] != len(names):
            raise ValueError(f"dataset {dataset_name!r} has inconsistent variable_names")
        for var_idx, name in enumerate(names):
            if name in variable_map:
                raise ValueError(f"duplicate dump variable name {name!r}")
            variable_map[name] = (dataset_name, var_idx)
            variable_order.append(name)
    if not variable_map:
        raise ValueError("dump file has no grouped datasets with variable_names attributes")
    return variable_map, variable_order


def unsigned_view(arr):
    if arr.dtype.itemsize == 4:
        return arr.view(np.uint32)
    if arr.dtype.itemsize == 8:
        return arr.view(np.uint64)
    raise TypeError(f"unsupported dump dtype for bitwise comparison: {arr.dtype}")


def read_variable_block(h5, variable_map, field, block_idx):
    dataset_name, var_idx = variable_map[field]
    return h5[dataset_name][block_idx, var_idx]


def read_variable_value(h5, variable_map, field, block_idx, cell_idx):
    dataset_name, var_idx = variable_map[field]
    return h5[dataset_name][(block_idx, var_idx) + tuple(cell_idx)]


def make_block_key(coord, level):
    # Round to avoid floating-point jitter in coordinate storage
    return (int(level), round(float(coord[0]), 10),
                         round(float(coord[1]), 10),
                         round(float(coord[2]), 10))


def compare_files(path1, path2):
    with h5py.File(path1, "r") as f1, h5py.File(path2, "r") as f2:
        coords1 = f1["coordinate"][:]          # (nblocks, 3)
        levels1 = f1["level"][:]               # (nblocks,)
        coords2 = f2["coordinate"][:]
        levels2 = f2["level"][:]
        variables1, order1 = build_variable_map(f1)
        variables2, _ = build_variable_map(f2)

        n1 = len(levels1)
        n2 = len(levels2)
        print(f"  File A: {n1} blocks")
        print(f"  File B: {n2} blocks")

        if n1 != n2:
            print(f"  ERROR: block count mismatch ({n1} vs {n2})")
            return False

        fields1 = set(variables1)
        fields2 = set(variables2)
        if fields1 != fields2:
            only_a = sorted(fields1 - fields2)
            only_b = sorted(fields2 - fields1)
            if only_a:
                print(f"  ERROR: variables only in file A: {only_a}")
            if only_b:
                print(f"  ERROR: variables only in file B: {only_b}")
            return False
        fields = [field for field in order1 if field in fields2]

        for field in fields:
            dset1_name, _ = variables1[field]
            dset2_name, _ = variables2[field]
            dset1 = f1[dset1_name]
            dset2 = f2[dset2_name]
            if dset1.dtype != dset2.dtype or dset1.shape[2:] != dset2.shape[2:]:
                print(f"  ERROR: variable {field} has incompatible storage")
                print(f"    A: {dset1_name}{dset1.shape} {dset1.dtype}")
                print(f"    B: {dset2_name}{dset2.shape} {dset2.dtype}")
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
            for field in fields:
                arr1 = read_variable_block(f1, variables1, field, b1)
                arr2 = read_variable_block(f2, variables2, field, b2)
                if arr1.shape != arr2.shape:
                    print(f"  ERROR: variable {field} block shape mismatch: {arr1.shape} vs {arr2.shape}")
                    return False
                diff_mask = unsigned_view(arr1) != unsigned_view(arr2)
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
            variables_a, _ = build_variable_map(fa)
            variables_b, _ = build_variable_map(fb)
            va = read_variable_value(fa, variables_a, field, b1_idx, cell_idx)
            vb = read_variable_value(fb, variables_b, field, b2_idx, cell_idx)
        print(f"    Value in A: {va!r}")
        print(f"    Value in B: {vb!r}")
        print(f"    Difference: {abs(va - vb):.6e}  (rel: {abs(va-vb)/max(abs(va),abs(vb),1e-300):.6e})")


if __name__ == "__main__":
    main()
