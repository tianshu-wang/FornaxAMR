#!/usr/bin/env python3
"""Visualize the FSA icosahedral angular grid as an orthographic SVG."""

import argparse
import math


FACES = (
    (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
    (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
    (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
    (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(a, s):
    return (s * a[0], s * a[1], s * a[2])


def normalize(a):
    mag = math.sqrt(dot(a, a))
    if mag <= 0.0:
        raise ValueError("zero-length vector")
    return (a[0] / mag, a[1] / mag, a[2] / mag)


def clamp_unit(x):
    return max(-1.0, min(1.0, x))


def icosahedron_vertices():
    phi = 0.5 * (1.0 + math.sqrt(5.0))
    return [
        normalize(v) for v in (
            (-1.0,  phi,  0.0), ( 1.0,  phi,  0.0),
            (-1.0, -phi,  0.0), ( 1.0, -phi,  0.0),
            ( 0.0, -1.0,  phi), ( 0.0,  1.0,  phi),
            ( 0.0, -1.0, -phi), ( 0.0,  1.0, -phi),
            ( phi,  0.0, -1.0), ( phi,  0.0,  1.0),
            (-phi,  0.0, -1.0), (-phi,  0.0,  1.0),
        )
    ]


def flat_face_point(ico, face, level, i, j):
    w0 = level - i - j
    w1 = i
    w2 = j
    return (
        (w0 * ico[face[0]][0] + w1 * ico[face[1]][0] + w2 * ico[face[2]][0]) / level,
        (w0 * ico[face[0]][1] + w1 * ico[face[1]][1] + w2 * ico[face[2]][1]) / level,
        (w0 * ico[face[0]][2] + w1 * ico[face[1]][2] + w2 * ico[face[2]][2]) / level,
    )


def find_or_add_vertex(vertices, cells, p):
    p = normalize(p)
    for idx, q in enumerate(vertices):
        if sum((p[d] - q[d]) ** 2 for d in range(3)) < 1.0e-24:
            return idx
    vertices.append(p)
    cells.append([])
    return len(vertices) - 1


def triangle_center(ico, face, level, coords):
    pts = [normalize(flat_face_point(ico, face, level, i, j)) for i, j in coords]
    e01 = sub(pts[0], pts[1])
    e02 = sub(pts[0], pts[2])
    center = cross(e01, e02)
    if dot(center, add(add(pts[0], pts[1]), pts[2])) < 0.0:
        center = scale(center, -1.0)
    return normalize(center)


def build_grid(level):
    ico = icosahedron_vertices()
    vertices = []
    triangles = []
    cells = []

    for face in FACES:
        local = {}
        for i in range(level + 1):
            for j in range(level + 1 - i):
                p = flat_face_point(ico, face, level, i, j)
                local[(i, j)] = find_or_add_vertex(vertices, cells, p)

        for i in range(level):
            for j in range(level - i):
                ids = (local[(i, j)], local[(i + 1, j)], local[(i, j + 1)])
                center = triangle_center(ico, face, level, ((i, j), (i + 1, j), (i, j + 1)))
                triangles.append((ids, center))
                tri_id = len(triangles) - 1
                for vertex_id in ids:
                    cells[vertex_id].append(tri_id)

                if i + j < level - 1:
                    ids = (local[(i + 1, j)], local[(i + 1, j + 1)], local[(i, j + 1)])
                    center = triangle_center(
                        ico, face, level, ((i + 1, j), (i + 1, j + 1), (i, j + 1))
                    )
                    triangles.append((ids, center))
                    tri_id = len(triangles) - 1
                    for vertex_id in ids:
                        cells[vertex_id].append(tri_id)

    expected_vertices = 10 * level * level + 2
    expected_triangles = 20 * level * level
    if len(vertices) != expected_vertices or len(triangles) != expected_triangles:
        raise RuntimeError(
            f"unexpected grid size: {len(vertices)} cells, {len(triangles)} triangles"
        )
    return vertices, triangles, cells


def build_arcs(vertices, triangles):
    edge_map = {}
    arcs = []

    for tri_id, (ids, _center) in enumerate(triangles):
        for c0, c1 in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            key = tuple(sorted((c0, c1)))
            edge_map.setdefault(key, []).append(tri_id)

    for (c0, c1), tri_ids in sorted(edge_map.items()):
        if len(tri_ids) != 2:
            raise RuntimeError(f"edge {(c0, c1)} has {len(tri_ids)} triangle neighbors")
        p0 = triangles[tri_ids[0]][1]
        p1 = triangles[tri_ids[1]][1]
        midpoint = normalize(add(p0, p1))
        delta = sub(vertices[c1], vertices[c0])
        vec = sub(delta, scale(midpoint, dot(delta, midpoint)))
        vec = normalize(vec)
        if dot(vec, delta) < 0.0:
            vec = scale(vec, -1.0)
        arcs.append({
            "cells": (c0, c1),
            "triangles": tuple(tri_ids),
            "angle": math.acos(clamp_unit(dot(p0, p1))),
            "midpoint": midpoint,
            "vec": vec,
        })

    expected_arcs = 30 * len(triangles) // 20
    if len(arcs) != expected_arcs:
        raise RuntimeError(f"unexpected arc count: {len(arcs)} != {expected_arcs}")
    return arcs


def tangent_basis(center):
    ref = (1.0, 0.0, 0.0) if abs(center[2]) > 0.9 else (0.0, 0.0, 1.0)
    e1 = normalize(cross(ref, center))
    e2 = cross(center, e1)
    return e1, e2


def ordered_cell_triangles(center, cell_triangles, triangles):
    e1, e2 = tangent_basis(center)
    return sorted(
        cell_triangles,
        key=lambda tri_id: math.atan2(dot(triangles[tri_id][1], e2), dot(triangles[tri_id][1], e1)),
    )


def slerp(a, b, t):
    cos_omega = max(-1.0, min(1.0, dot(a, b)))
    omega = math.acos(cos_omega)
    if omega < 1.0e-14:
        return a
    sin_omega = math.sin(omega)
    wa = math.sin((1.0 - t) * omega) / sin_omega
    wb = math.sin(t * omega) / sin_omega
    return normalize((
        wa * a[0] + wb * b[0],
        wa * a[1] + wb * b[1],
        wa * a[2] + wb * b[2],
    ))


def camera_basis():
    view = normalize((0.18, -0.38, 1.0))
    right = normalize(cross((0.0, 1.0, 0.0), view))
    up = cross(view, right)
    return view, right, up


def project(p, view, right, up, center, radius):
    return (
        center[0] + radius * dot(p, right),
        center[1] - radius * dot(p, up),
        dot(p, view),
    )


def path_data(points):
    if len(points) < 2:
        return ""
    return "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in points)


def edge_segments(a, b, want_front, view, right, up, center, radius):
    projected = []
    segments = []
    for step in range(33):
        p = slerp(a, b, step / 32.0)
        x, y, z = project(p, view, right, up, center, radius)
        is_front = z >= 0.0
        if is_front == want_front:
            projected.append((x, y))
        else:
            if len(projected) >= 2:
                segments.append(projected)
            projected = []
    if len(projected) >= 2:
        segments.append(projected)
    return segments


def write_svg(path, level):
    vertices, triangles, cells = build_grid(level)
    arcs = build_arcs(vertices, triangles)
    view, right, up = camera_basis()
    size = 900
    radius = 360
    center = (size / 2, size / 2)
    dual_points = [center_dir for _, center_dir in triangles]

    with open(path, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" viewBox="0 0 {size} {size}">\n')
        fh.write('<defs><marker id="arc-arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M 0 0 L 6 3 L 0 6 z" fill="#117a65"/></marker></defs>\n')
        fh.write('<rect width="100%" height="100%" fill="#ffffff"/>\n')
        fh.write(f'<circle cx="{center[0]}" cy="{center[1]}" r="{radius}" fill="#eee7f5" stroke="#dfd0ee" stroke-width="6"/>\n')

        for want_front, color, width, opacity in (
            (False, "#d7c5ea", 5, 0.45),
            (True, "#8e5cc2", 6, 0.92),
        ):
            for arc in arcs:
                i, j = arc["triangles"]
                a = dual_points[i]
                b = dual_points[j]
                for segment in edge_segments(a, b, want_front, view, right, up, center, radius):
                    fh.write(
                        f'<path d="{path_data(segment)}" fill="none" stroke="{color}" '
                        f'stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round" '
                        f'opacity="{opacity}"/>\n'
                    )

        arrow_half_angle = 0.035
        for arc in arcs:
            midpoint = arc["midpoint"]
            vec = arc["vec"]
            _xm, _ym, zm = project(midpoint, view, right, up, center, radius)
            if zm < 0.0:
                continue
            p_start = normalize(add(scale(midpoint, math.cos(arrow_half_angle)),
                                    scale(vec, -math.sin(arrow_half_angle))))
            p_end = normalize(add(scale(midpoint, math.cos(arrow_half_angle)),
                                  scale(vec, math.sin(arrow_half_angle))))
            x0, y0, z0 = project(p_start, view, right, up, center, radius)
            x1, y1, z1 = project(p_end, view, right, up, center, radius)
            if z0 < -0.02 or z1 < -0.02:
                continue
            fh.write(
                f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
                f'stroke="#117a65" stroke-width="2.6" stroke-linecap="round" '
                f'opacity="0.82" marker-end="url(#arc-arrow)"/>\n'
            )

        for p in dual_points:
            x, y, z = project(p, view, right, up, center, radius)
            if z >= 0.0:
                fh.write(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.6" fill="#c8addf" opacity="0.7"/>\n')

        fh.write(
            f'<text x="32" y="{size - 34}" font-family="Arial, sans-serif" '
            f'font-size="24" fill="#5e417d">N_ANGLE_LEV={level}, Nang={len(vertices)}</text>\n'
        )
        fh.write('</svg>\n')

    total_area = 4.0 * math.pi
    degree = [len(cell) for cell in cells]
    pent_hex = sum(
        1 for arc in arcs
        if sorted((degree[arc["cells"][0]], degree[arc["cells"][1]])) == [5, 6]
    )
    hex_hex = sum(
        1 for arc in arcs
        if degree[arc["cells"][0]] == 6 and degree[arc["cells"][1]] == 6
    )
    pent_pent = sum(
        1 for arc in arcs
        if degree[arc["cells"][0]] == 5 and degree[arc["cells"][1]] == 5
    )
    print(f"wrote {path}")
    print(
        f"N_ANGLE_LEV={level} Nang={len(vertices)} Narc={len(arcs)} "
        f"pent-pent={pent_pent} pent-hex={pent_hex} hex-hex={hex_hex} "
        f"expected_total_solid_angle={total_area:.16e}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--level", type=int, default=1, help="positive FSA angular level")
    parser.add_argument("--out", required=True, help="output SVG path")
    args = parser.parse_args()
    if args.level < 1:
        raise SystemExit("--level must be positive")
    write_svg(args.out, args.level)


if __name__ == "__main__":
    main()
