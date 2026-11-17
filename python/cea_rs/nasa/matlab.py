"""NASA MATLAB wrapper signatures; these wrappers are not ported."""
from __future__ import annotations

from types import SimpleNamespace
from typing import Optional, Sequence, SupportsFloat
import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
VectorLike = Sequence[SupportsFloat] | FloatArray
ScalarOrArray = SupportsFloat | VectorLike
__all__ = ["eq_solve", "rocket_solve", "shock_solve", "detonation_solve"]

def eq_solve(
    eq_type: int,
    reactants: list[str],
    T: Optional[float] = None,
    H: Optional[float] = None,
    S: Optional[float] = None,
    U: Optional[float] = None,
    P: Optional[float] = None,
    V: Optional[float] = None,
    T_reac: Optional[ScalarOrArray] = None,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[ScalarOrArray] = None,
    phi: Optional[ScalarOrArray] = None,
    r_eq: Optional[ScalarOrArray] = None,
    pct_fuel: Optional[ScalarOrArray] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
) -> SimpleNamespace:
    raise NotImplementedError('MATLAB eq_solve wrapper')

def rocket_solve(
    reactants: list[str],
    pc: SupportsFloat,
    pi_p: Optional[ScalarOrArray] = None,
    subar: Optional[ScalarOrArray] = None,
    supar: Optional[ScalarOrArray] = None,
    T_reac: Optional[ScalarOrArray] = None,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
    iac: bool = True,
    n_frz: Optional[int] = None,
    hc: Optional[SupportsFloat] = None,
    tc: Optional[SupportsFloat] = None,
    mdot: Optional[SupportsFloat] = None,
    ac_at: Optional[SupportsFloat] = None,
    tc_est: Optional[SupportsFloat] = None,
) -> SimpleNamespace:
    raise NotImplementedError('MATLAB rocket_solve wrapper')

def shock_solve(
    reactants: list[str],
    T0: SupportsFloat,
    p0: SupportsFloat,
    *,
    u1: Optional[SupportsFloat | VectorLike] = None,
    Mach1: Optional[SupportsFloat | VectorLike] = None,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
    reflected: bool = True,
    incident_frozen: bool = False,
    reflected_frozen: bool = False,
) -> SimpleNamespace:
    raise NotImplementedError('MATLAB shock_solve wrapper')

def detonation_solve(
    reactants: list[str],
    T1: SupportsFloat,
    p1: SupportsFloat,
    *,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
    frozen: bool = False,
) -> SimpleNamespace:
    raise NotImplementedError('MATLAB detonation_solve wrapper')

def __getattr__(name):
    if name.startswith("Matlab"):
        raise NotImplementedError(f"MATLAB {name}")
    raise AttributeError(name)
