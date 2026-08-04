from __future__ import annotations

import math
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version
from pathlib import Path

from geoclide.mathope import gamma_f32, gamma_f64

try:
    import tomllib
except ModuleNotFoundError:  # python 3.10
    tomllib = None

DIR_ROOT: Path = Path(__file__).resolve().parent.parent


def _get_version() -> str:
    """
    Read the geoclide version

    The version is read from pyproject.toml when running from a
    source checkout (python >= 3.11), else from the installed
    package metadata, with '0.0.0' as final fallback.
    """
    if tomllib is not None:
        try:
            with open(DIR_ROOT / "pyproject.toml", "rb") as f:
                return tomllib.load(f)["project"]["version"]
        except (FileNotFoundError, KeyError):
            pass
    try:
        return _installed_version("geoclide")
    except PackageNotFoundError:
        return "0.0.0"


VERSION: str = _get_version()

GAMMA2_F32 = gamma_f32(2)
GAMMA3_F32 = gamma_f32(3)
GAMMA5_F32 = gamma_f32(5)

GAMMA2_F64 = gamma_f64(2)
GAMMA3_F64 = gamma_f64(3)
GAMMA5_F64 = gamma_f64(5)

TWO_PI = math.pi * 2.0
