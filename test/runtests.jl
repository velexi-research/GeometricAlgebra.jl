"""
Unit tests for GeometricAlgebra.jl package.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
using Test

# External packages
using TestTools: jltest

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Run tests

jltest.run_tests(@__DIR__)
