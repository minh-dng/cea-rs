"""The separate NASA MATLAB convenience wrappers are outside this adapter."""

__all__ = ["eq_solve", "rocket_solve", "shock_solve", "detonation_solve"]


def __getattr__(name):
    if name in __all__ or name.startswith("Matlab"):
        raise NotImplementedError(f"cea_rs.nasa.matlab.{name}: MATLAB wrappers")
    raise AttributeError(name)
