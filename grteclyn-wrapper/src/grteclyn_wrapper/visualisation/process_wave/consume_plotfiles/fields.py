from __future__ import annotations

import numpy as np


def _canonical_field_name(name: str) -> str:
    # Accept a few convenient aliases.
    aliases = {
        "Weyl": "Weyl4_Re",
        "Weyl4": "Weyl4_Re",
        "Weyl_Re": "Weyl4_Re",
        "Weyl_Im": "Weyl4_Im",
        "Weyl_Mag": "Weyl4_Mag",
        # Complex scalar: imaginary part lives in lump slot 0 (StateVariables.hpp).
        "phi2": "phi_lump0",
        "Pi2": "Pi_lump0",
    }
    return aliases.get(name, name)


def _field_type_for(ds, *field_names: str) -> str:
    """Prefer BoxLib fields, but support yt stream datasets for isolated tests."""
    field_list = list(getattr(ds, "field_list", []))
    derived_field_list = list(getattr(ds, "derived_field_list", []))
    available = field_list + derived_field_list
    for name in field_names:
        if ("boxlib", name) in available:
            return "boxlib"
    for name in field_names:
        if ("stream", name) in available:
            return "stream"
    for name in field_names:
        for ftype, fname in available:
            if fname == name:
                return ftype
    return "boxlib"


def _field_key(ds, field_name: str) -> tuple[str, str]:
    return (_field_type_for(ds, field_name), field_name)


def _register_derived_fields(ds, field: str) -> None:
    base_ftype = _field_type_for(ds, "phi", "Pi", "chi")
    if (base_ftype, field) in list(getattr(ds, "field_list", [])) + list(getattr(ds, "derived_field_list", [])):
        return

    available_names = {
        name for (_ftype, name) in list(getattr(ds, "field_list", []))
        + list(getattr(ds, "derived_field_list", []))
    }
    lump_pairs = [
        (f"phi_lump{k}", f"Pi_lump{k}")
        for k in range(5)
        if f"phi_lump{k}" in available_names and f"Pi_lump{k}" in available_names
    ]

    def _sum_lump_component(data, component: str):
        total = None
        for k in range(5):
            name = f"{component}_lump{k}"
            if name not in available_names:
                continue
            term = data[base_ftype, name]
            total = term if total is None else total + term
        return total

    # GW proxy fields
    if field == "GW_Plus":
        def _gw_plus(field, data):
            return data[base_ftype, "A11"] - data[base_ftype, "A22"]
        ds.add_field((base_ftype, "GW_Plus"), function=_gw_plus, sampling_type="cell", units="")
    elif field == "GW_Cross":
        def _gw_cross(field, data):
            return 2.0 * data[base_ftype, "A12"]
        ds.add_field((base_ftype, "GW_Cross"), function=_gw_cross, sampling_type="cell", units="")
    elif field == "weyl4":
        # A_ij GW strain proxy (same formula as central_timeseries extraction).
        def _weyl4_proxy(field, data):
            a11 = data[base_ftype, "A11"]
            a12 = data[base_ftype, "A12"]
            a22 = data[base_ftype, "A22"]
            return np.sqrt((a11 - a22) ** 2 + (2.0 * a12) ** 2)
        ds.add_field((base_ftype, "weyl4"), function=_weyl4_proxy, sampling_type="cell", units="")
    elif field == "Weyl4_Mag":
        def _weyl4_mag(field, data):
            re_v = data[base_ftype, "Weyl4_Re"]
            im_v = data[base_ftype, "Weyl4_Im"]
            return np.sqrt(re_v**2 + im_v**2)
        ds.add_field((base_ftype, "Weyl4_Mag"), function=_weyl4_mag, sampling_type="cell", units="")
    elif field == "chi_minus_1":
        def _chi_minus_1(field, data):
            return data[base_ftype, "chi"] - 1.0
        ds.add_field((base_ftype, "chi_minus_1"), function=_chi_minus_1, sampling_type="cell", units="")
    elif field == "phi_lump_sum":
        def _phi_lump_sum(field, data):
            total = _sum_lump_component(data, "phi")
            if total is not None:
                return total
            return data[base_ftype, "phi"]
        ds.add_field((base_ftype, "phi_lump_sum"), function=_phi_lump_sum, sampling_type="cell", units="")
    elif field == "Pi_lump_sum":
        def _pi_lump_sum(field, data):
            total = _sum_lump_component(data, "Pi")
            if total is not None:
                return total
            return data[base_ftype, "Pi"]
        ds.add_field((base_ftype, "Pi_lump_sum"), function=_pi_lump_sum, sampling_type="cell", units="")
    elif field == "scalar_activity":
        def _scalar_activity(field, data):
            phi = data[base_ftype, "phi"]
            pi = data[base_ftype, "Pi"]
            # U(1) complex scalar (boson star): phi/Pi = Re(Phi, Pi), phi_lump0/Pi_lump0
            # = Im components in the first lump slots — combine into one |field| norm.
            if (
                lump_pairs == [("phi_lump0", "Pi_lump0")]
                and "phi" in available_names
                and "Pi" in available_names
            ):
                phi2 = data[base_ftype, "phi_lump0"]
                pi2 = data[base_ftype, "Pi_lump0"]
                return np.sqrt(phi**2 + pi**2 + phi2**2 + pi2**2)
            if lump_pairs:
                total = None
                for phi_name, pi_name in lump_pairs:
                    term = np.sqrt(
                        data[base_ftype, phi_name] ** 2 + data[base_ftype, pi_name] ** 2
                    )
                    total = term if total is None else total + term
                return total
            return np.sqrt(phi**2 + pi**2)
        ds.add_field((base_ftype, "scalar_activity"), function=_scalar_activity, sampling_type="cell", units="")
    elif field == "lump_activity":
        has_combined_scalar = "phi" in available_names and "Pi" in available_names
        ref_field = next(
            (name for name in ("rho_req", "chi", "K", "Weyl4_Re") if name in available_names),
            None,
        )

        def _lump_activity(field, data):
            if lump_pairs:
                total = None
                for phi_name, pi_name in lump_pairs:
                    term = np.sqrt(
                        data[base_ftype, phi_name] ** 2 + data[base_ftype, pi_name] ** 2
                    )
                    total = term if total is None else total + term
                return total
            if has_combined_scalar:
                phi = data[base_ftype, "phi"]
                pi = data[base_ftype, "Pi"]
                return np.sqrt(phi**2 + pi**2)
            if ref_field is None:
                raise RuntimeError(
                    "lump_activity requires phi_lump*/Pi_lump*, phi/Pi, or a reference field"
                )
            return np.zeros_like(data[base_ftype, ref_field])

        ds.add_field((base_ftype, "lump_activity"), function=_lump_activity, sampling_type="cell", units="")
    elif field == "local_speed":
        def _local_speed(field, data):
            eps = 1.0e-12
            chi = data[base_ftype, "chi"]
            lapse = data[base_ftype, "lapse"]
            c1 = np.abs(data[base_ftype, "shift1"]) + lapse * np.sqrt(chi / np.maximum(data[base_ftype, "h11"], eps))
            c2 = np.abs(data[base_ftype, "shift2"]) + lapse * np.sqrt(chi / np.maximum(data[base_ftype, "h22"], eps))
            c3 = np.abs(data[base_ftype, "shift3"]) + lapse * np.sqrt(chi / np.maximum(data[base_ftype, "h33"], eps))
            return np.maximum(np.maximum(c1, c2), c3)
        ds.add_field((base_ftype, "local_speed"), function=_local_speed, sampling_type="cell", units="")

    elif field == "Weyl4_Re":
        # Ensure base fields are available if asked for explicitly?
        # Usually they are just read from disk.
        pass
