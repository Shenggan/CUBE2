#!/usr/bin/env python3
"""Generate CUDA compile-time layout constants from parameters.f90."""

import re
import sys
from pathlib import Path


REQUIRED = (
    "ng",
    "nnt",
    "nns",
    "ratio_cs",
    "np_nc",
    "ngb",
    "ncore",
    "body_centered_cubic",
    "image_buffer",
    "tile_buffer",
    "app",
)


def read_parameters(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.split("!", 1)[0]
        if "parameter" not in line.lower():
            continue
        for name in REQUIRED:
            match = re.search(rf"\b{name}\s*=\s*([^,\s]+)", line, re.IGNORECASE)
            if match:
                values[name] = match.group(1)
    missing = [name for name in REQUIRED if name not in values]
    if missing:
        raise ValueError(f"missing parameter(s): {', '.join(missing)}")
    return values


def resolve(value: str, values: dict[str, str]) -> str:
    seen = set()
    while value.lower() in values and value.lower() not in seen:
        seen.add(value.lower())
        value = values[value.lower()]
    return value


def cpp_real(value: str, values: dict[str, str]) -> str:
    value = resolve(value, values).replace("D", "e").replace("d", "e")
    if not re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?", value):
        raise ValueError(f"not a numeric real expression: {value}")
    return value


def cpp_int(value: str, values: dict[str, str]) -> str:
    value = resolve(value, values)
    if not re.fullmatch(r"[+-]?\d+", value):
        raise ValueError(f"not an integer expression: {value}")
    return value


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: gen_cuda_parameters.py parameters.f90 cuda_parameters_generated.h")

    values = {key.lower(): value for key, value in read_parameters(Path(sys.argv[1])).items()}
    body_centered = resolve(values["body_centered_cubic"], values).lower()
    if body_centered not in (".true.", ".false."):
        raise ValueError(f"not a logical constant: {body_centered}")

    header = f"""// Generated from parameters.f90; do not edit manually.
#pragma once

constexpr int CUBE_NG = {cpp_int(values['ng'], values)};
constexpr int CUBE_NNT = {cpp_int(values['nnt'], values)};
constexpr int CUBE_NNS = {cpp_int(values['nns'], values)};
constexpr int CUBE_RATIO_CS = {cpp_int(values['ratio_cs'], values)};
constexpr int CUBE_NP_NC = {cpp_int(values['np_nc'], values)};
constexpr int CUBE_NGB = {cpp_int(values['ngb'], values)};
constexpr int CUBE_NCORE = {cpp_int(values['ncore'], values)};
constexpr bool CUBE_BODY_CENTERED_CUBIC = {'true' if body_centered == '.true.' else 'false'};
constexpr double CUBE_IMAGE_BUFFER = {cpp_real(values['image_buffer'], values)};
constexpr double CUBE_TILE_BUFFER = {cpp_real(values['tile_buffer'], values)};
constexpr double CUBE_APP = {cpp_real(values['app'], values)};
"""
    Path(sys.argv[2]).write_text(header)


if __name__ == "__main__":
    main()
