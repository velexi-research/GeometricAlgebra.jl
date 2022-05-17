#   Copyright 2020 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
operations.jl defines operations on GeometricAlgebra types
"""

# Additive operations
include("operations/add.jl")
include("operations/subtract.jl")

# Multiplicative operations
include("operations/multiply.jl")
include("operations/divide.jl")
include("operations/wedge.jl")
include("operations/contract_left.jl")
include("operations/dot.jl")

# Geometric operations
include("operations/dual.jl")
include("operations/project.jl")
include("operations/reject.jl")
