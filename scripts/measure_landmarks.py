#!/usr/bin/env python3
"""Calculate distances from supplied native-image landmarks and DICOM geometry.

Standard library only. Does not read DICOM, select anatomy, validate the source
of metadata, or establish diagnostic accuracy. Input coordinates are zero-based
pixel centers [column, row], including fractional positions between centers.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


def finite(value, name):
    if not math.isfinite(value):
        raise ValueError(f"{name}: non-finite result; verify coordinate magnitudes")
    return value


def finite_sum(values, name):
    try:
        return finite(math.fsum(values), name)
    except OverflowError as error:
        raise ValueError(f"{name}: numeric overflow; verify coordinate magnitudes") from error


def validate_decimals(decimals):
    if type(decimals) is not int or not 0 <= decimals <= 4:
        raise ValueError("decimals must be an integer from 0 to 4")


def vector(value, size, name):
    if not isinstance(value, (list, tuple)) or len(value) != size:
        raise ValueError(f"{name}: expected {size} finite numbers")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in value):
        raise ValueError(f"{name}: values must be numbers, not strings or booleans")
    try:
        result = tuple(float(x) for x in value)
    except OverflowError as error:
        raise ValueError(f"{name}: number outside floating-point range") from error
    if not all(math.isfinite(x) for x in result):
        raise ValueError(f"{name}: non-finite value")
    return result


def dot(a, b):
    return finite_sum((finite(x * y, "dot product term") for x, y in zip(a, b)),
                      "dot product")


def sub(a, b):
    return tuple(finite(x - y, "coordinate difference") for x, y in zip(a, b))


def norm(a):
    return finite(math.hypot(*a), "vector magnitude")


def cross(a, b):
    return tuple(finite_sum((finite(a[i] * b[j], "cross product term"),
                             -finite(a[j] * b[i], "cross product term")),
                            "cross product")
                 for i, j in ((1, 2), (2, 0), (0, 1)))


def calculate(record, decimals=1):
    validate_decimals(decimals)
    if not isinstance(record, dict):
        raise ValueError("measurement must be an object")
    if not isinstance(record.get("id"), str) or not record["id"].strip():
        raise ValueError("measurement requires a nonempty string id")
    spacing = vector(record.get("pixel_spacing_mm"), 2, "pixel_spacing_mm")
    if min(spacing) <= 0:
        raise ValueError("pixel_spacing_mm must be positive")
    ipp = vector(record.get("ipp_mm"), 3, "ipp_mm")
    iop = vector(record.get("iop"), 6, "iop")
    u, v = iop[:3], iop[3:]
    if (abs(norm(u) - 1) > 1e-4 or abs(norm(v) - 1) > 1e-4
            or abs(dot(u, v)) > 1e-4):
        raise ValueError("IOP vectors must be unit and orthogonal; verify metadata")
    rows, columns = record.get("rows"), record.get("columns")
    if rows is not None or columns is not None:
        if any(type(x) is not int or x < 1 for x in (rows, columns)):
            raise ValueError("rows and columns must both be positive integers")

    def xyz(point):
        c, r = vector(point, 2, "landmark [column,row]")
        if c < 0 or r < 0:
            raise ValueError("landmarks must be on the source grid, not negative")
        if rows is not None and (c > columns - 1 or r > rows - 1):
            raise ValueError("landmark outside source pixel-center bounds")
        column_mm = finite(c * spacing[1], "column offset")
        row_mm = finite(r * spacing[0], "row offset")
        return tuple(finite_sum((ipp[j],
                                 finite(column_mm * u[j], "column coordinate"),
                                 finite(row_mm * v[j], "row coordinate")),
                                "patient coordinate")
                     for j in range(3))

    def two_points(key):
        points = record.get(key)
        if not isinstance(points, list) or len(points) != 2:
            raise ValueError(f"{key}: expected exactly two [column,row] points")
        return list(map(xyz, points))

    kind = record.get("kind")
    result = {"id": record["id"], "kind": kind, "input_record": record}
    if kind == "perpendicular":
        a, b = two_points("line")
        p = xyz(record.get("point"))
        direction = sub(b, a)
        size = norm(direction)
        if size <= 1e-12:
            raise ValueError("reference line has coincident endpoints")
        # Normalizing first avoids squaring a large length or forming inf/inf.
        unit = tuple(finite(component / size, "unit direction") for component in direction)
        offset = sub(p, a)
        along = dot(offset, unit)
        value = norm(cross(offset, unit))
        parameter = finite(along / size, "projection parameter")
        foot = tuple(finite_sum((a[j], finite(along * unit[j], "foot offset")),
                                "perpendicular foot") for j in range(3))
        result.update({"landmarks_patient_coordinates_mm": [a, b, p],
                       "projection_parameter": parameter,
                       "perpendicular_foot_patient_mm": foot,
                       "foot_outside_reference_segment": not 0 <= parameter <= 1})
    elif kind in ("length", "patient_axis_projection"):
        a, b = two_points("points")
        difference = sub(b, a)
        result["landmarks_patient_coordinates_mm"] = [a, b]
        if kind == "length":
            value = norm(difference)
        else:
            axis = vector(record.get("axis"), 3, "axis")
            if abs(norm(axis) - 1) > 1e-4:
                raise ValueError("axis must be a unit vector in patient coordinates")
            signed = dot(difference, axis)
            result["signed_axis_difference_mm"] = signed
            value = abs(signed)
    else:
        raise ValueError(f"unsupported measurement kind: {kind!r}")
    result["computed_value_mm"] = finite(value, "measurement")
    result["display_mm"] = finite(round(value, decimals), "rounded measurement")
    result["calculation_only"] = True
    return result


def calculate_document(document, decimals=1):
    validate_decimals(decimals)
    if not isinstance(document, dict):
        raise ValueError("input must be an object")
    if document.get("coordinate_convention") != "zero_based_pixel_centers":
        raise ValueError("explicit coordinate_convention=zero_based_pixel_centers required")
    records = document.get("measurements")
    if not isinstance(records, list) or not records:
        raise ValueError("measurements must be a nonempty array")
    results, seen = [], set()
    for index, record in enumerate(records):
        try:
            result = calculate(record, decimals)
            study = record.get("study", "")
            if not isinstance(study, str):
                raise ValueError("study, when present, must be a string")
            identity = (study, record["id"])
            if identity in seen:
                raise ValueError("duplicate measurement id within study")
            seen.add(identity)
            results.append(result)
        except (ValueError, OverflowError) as error:
            raise ValueError(f"measurement {index + 1}: {error}") from error
    return {"schema_version": 1,
            "coordinate_convention": "zero_based_pixel_centers",
            "limitations": "Arithmetic from supplied geometry only; source identity and anatomical boundaries are not validated.",
            "measurements": results}


def write_new_output(path, payload):
    """Publish complete UTF-8 bytes without replacing an existing destination.

    A same-directory staging file is linked into place only after a successful
    write. A filesystem that cannot create hard links fails without a partial
    destination. Staging files are private and removed on ordinary exceptions.
    """
    encoded = payload.encode("utf-8")
    staging = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", prefix=".landmarks-",
                                         suffix=".tmp", dir=path.parent,
                                         delete=False) as stream:
            staging = Path(stream.name)
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(staging, path)
    finally:
        if staging is not None:
            staging.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, help="new file only; existing files are never replaced")
    parser.add_argument("--decimals", type=int, choices=range(5), default=1)
    args = parser.parse_args()
    try:
        raw = args.input.read_bytes()
        result = calculate_document(json.loads(raw), args.decimals)
        result["input_sha256"] = hashlib.sha256(raw).hexdigest()
        payload = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        if args.output:
            write_new_output(args.output, payload)
        else:
            print(payload, end="")
    except (OSError, ValueError, OverflowError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
