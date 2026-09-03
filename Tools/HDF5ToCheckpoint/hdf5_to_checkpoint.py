#!/usr/bin/env python3
"""Convert a Chombo-format HDF5 file (e.g. GRTresna's InitialDataFinal.3d.hdf5,
or a GRChombo checkpoint) into a GRTeclyn/AMReX checkpoint that can be loaded
with `amr.restart`.

GRTeclyn does not implement its own checkpoint I/O: it uses AMReX's stock
Amr/AmrLevel/StateData/VisMF machinery unchanged, and that on-disk format
isn't documented anywhere as a stable spec. This script reconstructs it from
first principles by reading AMReX's own serialisation code (Geometry/CoordSys/
RealBox/Box stream operators in amrex/Src/Base, VisMF's Header format in
amrex/Src/Base/AMReX_VisMF.cpp) rather than depending on a real checkpoint to
copy the format from. Everything it needs -- box layout, domain, cell size,
periodicity, ref_ratio, per-level dt -- comes straight out of the HDF5 file's
own per-level groups/attributes; the physical domain origin (prob_lo) isn't
stored in the Chombo file, so it defaults to (0, 0, 0) (true of every
GRChombo/GRTeclyn example so far) and can be overridden with --prob-lo.

The one thing that's a genuine simulation choice, not physics data -- the
number of ghost cells (amr.evolution.num_ghosts) -- isn't read from anywhere:
it's inferred directly from how much halo padding the HDF5 file itself stores
around each box, which is what must match anyway for GRTeclyn to interpret
the data correctly.

State variable components are transcribed positionally, straight from
whatever order the HDF5 file's own component_0, component_1, ... attributes
list -- no name matching, no hardcoded variable list. GRTeclyn addresses
state variables by position, not name, so this produces a correct checkpoint
only if the file's own component order already matches the order the target
executable's StateVariables.hpp expects. That's true for GRTresna/GRChombo
output in every case seen so far (component_0='chi' ... matching CCZ4's
order, with any matter fields appended after), but isn't independently
verified here -- if a source file ever used a different order, this would
transcribe it wrong with no error.

Usage
-----
    python3 hdf5_to_checkpoint.py InitialDataFinal.3d.hdf5 my_restart_chk

Then point your params file at the result:

    amr.restart = /path/to/my_restart_chk

Only cell-centered, single state descriptor data is supported (this matches
every current GRTeclyn example).
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

Box = tuple  # ((lo_i, lo_j, lo_k), (hi_i, hi_j, hi_k))

# AMReX's FABio::FAB_NATIVE format descriptor for an IEEE-754 double,
# amrex/Src/Base/AMReX_FPC.cpp: ieee_double = ("64 bits, sign at 0, 11
# exponent bits at 1, 52 mantissa bits at 12, bias 1023") is the same on
# every platform (IEEE-754 is universal); only the byte order differs, and
# FPC::NativeRealDescriptor() picks reverse_double_order on little-endian
# machines (x86_64, ARM64/Apple Silicon -- effectively everything today) and
# normal_double_order on big-endian ones.
_IEEE_DOUBLE = "8, (64 11 52 0 1 12 0 1023)"
_LITTLE_ENDIAN_ORDER = "8, (8 7 6 5 4 3 2 1)"
_BIG_ENDIAN_ORDER = "8, (1 2 3 4 5 6 7 8)"


def native_fab_format() -> str:
    order = _LITTLE_ENDIAN_ORDER if sys.byteorder == "little" else _BIG_ENDIAN_ORDER
    return f"FAB (({_IEEE_DOUBLE}),({order}))"


def box_ncells(box: Box) -> tuple[int, int, int]:
    lo, hi = box
    return (hi[0] - lo[0] + 1, hi[1] - lo[1] + 1, hi[2] - lo[2] + 1)


def box_text(box: Box) -> str:
    lo, hi = box
    return f"(({lo[0]},{lo[1]},{lo[2]}) ({hi[0]},{hi[1]},{hi[2]}) (0,0,0))"


def boxarray_text(boxes: list[Box]) -> str:
    return "\n".join([f"({len(boxes)} 0", *[box_text(b) for b in boxes], ")"])


def fmt_real(x: float) -> str:
    """Match C++'s default (non-scientific) ostream formatting with
    precision(17), which is what Amr::checkPoint() uses for everything in
    the top-level Header outside the VisMF min/max blocks."""
    return f"{x:.17g}"


def format_geom_entry(domain_box: Box, dx: float, prob_lo: tuple[float, float, float],
                       periodic: tuple[int, int, int]) -> str:
    """Reproduce amrex::operator<<(ostream&, const Geometry&), which is
    CoordSys << RealBox << Domain() << 'P' << IntVect(periodicity):
      amrex/Src/Base/AMReX_Geometry.cpp:18-22
      amrex/Src/Base/AMReX_CoordSys.cpp:446-459 (offset == prob_lo, see
        amrex/Src/Base/AMReX_Geometry.cpp:154/158/188, SetOffset(rb.lo()))
      amrex/Src/Base/AMReX_RealBox.cpp:31-40
    e.g. "(0 (0,0,0)(40,40,40) 1)\\n(RealBox 0 640 0 640 0 640 )((0,0,0) (15,15,15) (0,0,0))P(0,0,0)"
    """
    nx, ny, nz = box_ncells(domain_box)
    prob_hi = (prob_lo[0] + dx * nx, prob_lo[1] + dx * ny, prob_lo[2] + dx * nz)
    coordsys = (
        f"(0 ({fmt_real(prob_lo[0])},{fmt_real(prob_lo[1])},{fmt_real(prob_lo[2])})"
        f"({fmt_real(dx)},{fmt_real(dx)},{fmt_real(dx)}) 1)"
    )
    realbox = (
        f"(RealBox {fmt_real(prob_lo[0])} {fmt_real(prob_hi[0])} "
        f"{fmt_real(prob_lo[1])} {fmt_real(prob_hi[1])} "
        f"{fmt_real(prob_lo[2])} {fmt_real(prob_hi[2])} )"
    )
    per = f"P({int(periodic[0])},{int(periodic[1])},{int(periodic[2])})"
    return f"{coordsys}\n{realbox}{box_text(domain_box)}{per}"


def infer_ghost_width(box: Box, stored_cells: int, ncomp: int) -> int:
    nx, ny, nz = box_ncells(box)
    per_comp = stored_cells // ncomp
    if per_comp * ncomp != stored_cells:
        raise ValueError("stored cell count is not a multiple of num_components")
    for g in range(0, 10):
        if (nx + 2 * g) * (ny + 2 * g) * (nz + 2 * g) == per_comp:
            return g
    raise ValueError(
        f"could not infer a consistent ghost width for box {box} "
        f"with {per_comp} stored cells per component"
    )


def reghost(block: np.ndarray, box: Box, g_src: int, g_dst: int) -> np.ndarray:
    """Re-pad/crop a flat (Fortran-order) per-component block from ghost
    width g_src to g_dst. Only used when a box's stored halo width doesn't
    match the rest of the file (ghost width is otherwise inferred once and
    assumed uniform)."""
    nx, ny, nz = box_ncells(box)
    src = block.reshape((nx + 2 * g_src, ny + 2 * g_src, nz + 2 * g_src), order="F")
    if g_dst <= g_src:
        d = g_src - g_dst
        dst = src[d: d + nx + 2 * g_dst, d: d + ny + 2 * g_dst, d: d + nz + 2 * g_dst]
    else:
        d = g_dst - g_src
        dst = np.pad(src, d, mode="edge")
    return np.asfortranarray(dst).reshape(-1, order="F")


def read_hdf5_level(f: h5py.File, level: int, names: list[str], num_components: int):
    lvl = f[f"level_{level}"]
    raw_boxes = lvl["boxes"][:]
    boxes = [
        ((int(b["lo_i"]), int(b["lo_j"]), int(b["lo_k"])),
         (int(b["hi_i"]), int(b["hi_j"]), int(b["hi_k"])))
        for b in raw_boxes
    ]
    attrs = lvl.attrs
    domain = attrs["prob_domain"]
    domain_box = ((int(domain[0]), int(domain[1]), int(domain[2])),
                  (int(domain[3]), int(domain[4]), int(domain[5])))
    periodic = (int(attrs.get("is_periodic_0", 0)),
                int(attrs.get("is_periodic_1", 0)),
                int(attrs.get("is_periodic_2", 0)))
    return {
        "names": names,
        "num_components": num_components,
        "boxes": boxes,
        "offsets": lvl["data:offsets=0"][:],
        "data": lvl["data:datatype=0"][:],
        "domain_box": domain_box,
        "dx": float(attrs["dx"]),
        "dt": float(attrs.get("dt", 0.0)),
        "ref_ratio": int(attrs.get("ref_ratio", 2)),
        "periodic": periodic,
    }


def read_hdf5_source(path: Path):
    f = h5py.File(path, "r")
    if "Chombo_global" not in f:
        raise ValueError(f"{path} does not look like a Chombo HDF5 file (no Chombo_global group)")
    spacedim = int(f["Chombo_global"].attrs.get("SpaceDim", 3))
    if spacedim != 3:
        raise ValueError(f"only 3D data is supported, file has SpaceDim={spacedim}")

    num_levels = int(f.attrs["num_levels"])
    num_components = int(f.attrs["num_components"])
    names = [f.attrs[f"component_{i}"].decode() if isinstance(f.attrs[f"component_{i}"], bytes)
             else str(f.attrs[f"component_{i}"])
             for i in range(num_components)]

    levels = [read_hdf5_level(f, lvl, names, num_components) for lvl in range(num_levels)]
    return {"names": names, "num_components": num_components, "levels": levels}


def build_level_data(level: int, boxes: list[Box], src_level, ngrow: int,
                      fab_format: str, ncomp: int):
    """Assemble the VisMF _H text and _D binary payload for one level.

    Components are transcribed positionally, in whatever order the HDF5
    file itself stores them -- GRTeclyn addresses state variables by
    position, not name, so this only produces a correct checkpoint if the
    file's component order already matches the target executable's."""
    src_box_index = {b: i for i, b in enumerate(src_level["boxes"])}
    unmatched = [b for b in boxes if b not in src_box_index]
    if unmatched:
        raise SystemExit(
            f"error: level {level} boxes {unmatched} are not present in "
            f"the HDF5 file's level_{level} box list: {src_level['boxes']}"
        )

    out_bytes = bytearray()
    fod_lines = []
    box_min: list[list[float]] = []
    box_max: list[list[float]] = []
    offset = 0
    for box in boxes:
        nx, ny, nz = box_ncells(box)
        cells_dst = (nx + 2 * ngrow) * (ny + 2 * ngrow) * (nz + 2 * ngrow)
        grown_box = (
            (box[0][0] - ngrow, box[0][1] - ngrow, box[0][2] - ngrow),
            (box[1][0] + ngrow, box[1][1] + ngrow, box[1][2] + ngrow),
        )
        preamble = f"{fab_format}{box_text(grown_box)} {ncomp}\n".encode()

        src_bi = src_box_index[box]
        src_off0, src_off1 = src_level["offsets"][src_bi], src_level["offsets"][src_bi + 1]
        per_comp_src = (src_off1 - src_off0) // src_level["num_components"]
        g_src = infer_ghost_width(box, src_off1 - src_off0, src_level["num_components"])

        comp_min, comp_max, comp_blocks = [], [], []
        for ci in range(ncomp):
            block = np.asarray(
                src_level["data"][src_off0 + ci * per_comp_src: src_off0 + (ci + 1) * per_comp_src],
                dtype=np.float64,
            )
            if g_src != ngrow:
                block = reghost(block, box, g_src, ngrow)
            comp_min.append(float(block.min()))
            comp_max.append(float(block.max()))
            comp_blocks.append(block)

        box_min.append(comp_min)
        box_max.append(comp_max)
        payload = np.concatenate(comp_blocks).astype(np.float64, copy=False)
        assert payload.size == ncomp * cells_dst

        fod_lines.append(f"FabOnDisk: SD_0_New_MF_D_00000 {offset}")
        out_bytes.extend(preamble)
        out_bytes.extend(payload.tobytes())
        offset += len(preamble) + payload.nbytes

    def fmt_row(vals: list[float]) -> str:
        return "".join(f"{v:.17e}," for v in vals) + "\n"

    header_text = (
        f"1\n1\n{ncomp}\n{ngrow}\n{boxarray_text(boxes)}\n"
        f"{len(boxes)}\n" + "\n".join(fod_lines) + "\n\n"
        f"{len(boxes)},{ncomp}\n" + "".join(fmt_row(r) for r in box_min) +
        f"\n{len(boxes)},{ncomp}\n" + "".join(fmt_row(r) for r in box_max)
    )
    return header_text, bytes(out_bytes)


def build_top_header(num_levels: int, geom_entries: list[str],
                      ref_ratios: list[int], dt_levels: list[float]) -> str:
    """Reproduce the fixed preamble Amr::checkPoint() writes before the
    per-level AmrLevel/StateData blocks: amrex/Src/Amr/AMReX_Amr.cpp:1838-1859."""
    max_level = num_levels - 1
    finest_level = max_level
    header = "CheckPointVersion_1.0\n3\n0\n" + f"{max_level}\n{finest_level}\n"
    header += "".join(e + " " for e in geom_entries) + "\n"
    header += "".join(f"({r},{r},{r}) " for r in ref_ratios) + "\n"
    header += "".join(fmt_real(v) + " " for v in dt_levels) + "\n"  # dt_level
    header += "".join(fmt_real(v) + " " for v in dt_levels) + "\n"  # dt_min
    n_cycle = [1] + list(ref_ratios)
    header += "".join(f"{v} " for v in n_cycle) + "\n"
    header += "".join("0 " for _ in range(num_levels)) + "\n"  # level_steps
    header += "".join("0 " for _ in range(num_levels)) + "\n"  # level_count
    return header


def convert(hdf5_path: Path, output_dir: Path,
            prob_lo: tuple[float, float, float], force: bool) -> None:
    if output_dir.exists():
        if not force:
            raise SystemExit(
                f"error: output directory {output_dir} already exists "
                "(pass --force to overwrite)"
            )
        shutil.rmtree(output_dir)

    src = read_hdf5_source(hdf5_path)
    num_levels = len(src["levels"])
    ncomp = src["num_components"]

    # Ghost width is inferred once from level 0's first box and assumed
    # uniform everywhere -- it's a simulation choice (evolution.num_ghosts),
    # not physics data, so there's nowhere else to get it from except by
    # checking it matches what every box in the file actually stores.
    lvl0 = src["levels"][0]
    ngrow = infer_ghost_width(lvl0["boxes"][0], lvl0["offsets"][1] - lvl0["offsets"][0],
                              lvl0["num_components"])
    fab_format = native_fab_format()

    print(f"Converting {hdf5_path} -> {output_dir}")
    print(f"  levels: {num_levels}, ngrow: {ngrow}, variables ({ncomp}): {src['names']}")
    for lvl in range(num_levels):
        print(f"  level {lvl}: {len(src['levels'][lvl]['boxes'])} box(es)")

    output_dir.mkdir(parents=True)
    (output_dir / "FabArrayHeaders.txt").write_text(
        "".join(f"Level_{lvl}/SD_0_New_MF\n" for lvl in range(num_levels))
    )

    geom_entries = [
        format_geom_entry(src["levels"][lvl]["domain_box"], src["levels"][lvl]["dx"],
                          prob_lo, src["levels"][lvl]["periodic"])
        for lvl in range(num_levels)
    ]
    ref_ratios = [src["levels"][lvl]["ref_ratio"] for lvl in range(num_levels - 1)]
    dt_levels = [src["levels"][lvl]["dt"] for lvl in range(num_levels)]

    level_blocks = []
    for lvl in range(num_levels):
        boxes = src["levels"][lvl]["boxes"]
        state_header_text, data_bytes = build_level_data(
            lvl, boxes, src["levels"][lvl], ngrow, fab_format, ncomp,
        )
        level_dir = output_dir / f"Level_{lvl}"
        level_dir.mkdir()
        (level_dir / "SD_0_New_MF_H").write_text(state_header_text)
        (level_dir / "SD_0_New_MF_D_00000").write_bytes(data_bytes)

        ba_text = boxarray_text(boxes)
        domain_text = box_text(src["levels"][lvl]["domain_box"])
        dt = dt_levels[lvl]
        level_blocks.append(
            f"{lvl}\n{geom_entries[lvl]}\n{ba_text}1\n"
            f"{domain_text}\n{ba_text}{fmt_real(-dt)}\n{fmt_real(-dt)}\n0\n0\n1\n"
            f"Level_{lvl}/SD_0_New_MF\n"
        )

    header_text = build_top_header(num_levels, geom_entries, ref_ratios, dt_levels) + "".join(level_blocks)
    (output_dir / "Header").write_text(header_text)

    print(f"Done. Set amr.restart = {output_dir} in your params file.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hdf5_file", type=Path,
                         help="Chombo-format HDF5 file (GRTresna output or GRChombo checkpoint)")
    parser.add_argument("output_checkpoint", type=Path,
                         help="Where to write the new GRTeclyn-compatible checkpoint")
    parser.add_argument("--prob-lo", default="0,0,0",
                         help="Physical domain lower corner 'x,y,z' (default '0,0,0', true "
                              "of every GRChombo/GRTeclyn example seen so far). The Chombo "
                              "HDF5 file doesn't store this, only cell size and index-space "
                              "domain, so it has to come from you if it isn't the origin.")
    parser.add_argument("--force", action="store_true",
                         help="Overwrite output_checkpoint if it already exists")
    args = parser.parse_args()

    prob_lo = tuple(float(x) for x in args.prob_lo.split(","))
    if len(prob_lo) != 3:
        raise SystemExit(f"error: --prob-lo must have 3 components, got {args.prob_lo!r}")
    convert(args.hdf5_file, args.output_checkpoint, prob_lo, args.force)


if __name__ == "__main__":
    main()
