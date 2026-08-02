#!/usr/bin/env python3
"""Shared HDF5 extraction for the one-dimensional radiation benchmarks."""
import h5py
import numpy as np

BS_DEFAULT = 16
M_ID, M_ACTIVE = 0, 3
M_XMIN, M_DX = slice(4, 7), slice(10, 13)


def load_line(path, ngroup):
    with h5py.File(path, "r") as h5:
        bs = int(h5.attrs.get("block_size", BS_DEFAULT))
        scale = float(h5.attrs.get("rad_scale", 1.0))
        time = float(h5.attrs["time"])
        meta = h5["MetaData"][:]
        data = h5["Data"]
        eos = h5["Eos"]
        rows = []
        for bid, m in enumerate(meta):
            if m[M_ID] < 0 or m[M_ACTIVE] < 0.5:
                continue
            x = m[M_XMIN][0] + (np.arange(bs) + 0.5) * m[M_DX][0]
            jc = kc = bs // 2
            rho = np.asarray(data[bid, 0]).reshape(bs, bs, bs)[:, jc, kc]
            vel = np.asarray([
                np.asarray(data[bid, v]).reshape(bs, bs, bs)[:, jc, kc]
                for v in (1, 2, 3)
            ])
            temp = np.asarray(eos[bid, 1]).reshape(bs, bs, bs)[:, jc, kc]
            eg = []
            fg = []
            for g in range(ngroup):
                eg.append(np.asarray(data[bid, 6 + 4*g]).reshape(bs, bs, bs)[:, jc, kc] * scale)
                fg.append(np.asarray(data[bid, 7 + 4*g]).reshape(bs, bs, bs)[:, jc, kc] * scale)
            rows.append((x, rho, vel, temp, np.asarray(eg), np.asarray(fg)))
    order = np.argsort(np.concatenate([r[0] for r in rows]))
    return dict(time=time, x=np.concatenate([r[0] for r in rows])[order],
        rho=np.concatenate([r[1] for r in rows])[order],
        vel=np.concatenate([r[2] for r in rows], axis=1)[:, order],
        temp=np.concatenate([r[3] for r in rows])[order],
        E=np.concatenate([r[4] for r in rows], axis=1)[:, order],
        F=np.concatenate([r[5] for r in rows], axis=1)[:, order])
