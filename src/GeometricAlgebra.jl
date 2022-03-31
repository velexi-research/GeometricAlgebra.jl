#   Copyright (c) 2020-2022 Velexi Corporation
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
The GeometricAlgebra.jl module defines geometric algebra types and functions.
"""
module GeometricAlgebra

# Types
include("AbstractMultivector.jl")
include("AbstractBlade.jl")
include("AbstractScalar.jl")
include("Zero.jl")
include("One.jl")
include("Scalar.jl")
include("Blade.jl")
include("Pseudoscalar.jl")
include("Multivector.jl")

# Methods
include("operations.jl")

end  # End GeometricAlgebra module
