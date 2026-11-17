"""Re-export the compatibility API without loading a Fortran library."""

from .. import *
from .. import __all__
from .. import lib_version as _version
from .. import lib_version_major as _version_major
from .. import lib_version_minor as _version_minor
from .. import lib_version_patch as _version_patch
