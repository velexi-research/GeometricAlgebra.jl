"""
Utility functions and macros for unit tests.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

"""
    get_random_value(value_to_add=0)

Return a random value with absolute value in the range [`value_to_add`, `value_to_add` + 1).
If `value_to_add` is not provided, 0 will be used.
Value is either negative or positive with equal probability.
"""
macro get_random_value(value_to_add=0)
    return quote
        local random_number = rand() + $value_to_add
        (rand() > 0.5) ? random_number : -random_number
    end
end
