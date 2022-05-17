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
Utility functions for unit tests.
"""

"""
    get_random_value(value_to_add=0)

Return a random value with absolute value in the range [`value_to_add`, `value_to_add` + 1).
If `value_to_add` is not provided, 0 will be used.
Value is either negative or positive with equal probability.
"""
function get_random_value(value_to_add=0)
    random_number = rand() + value_to_add
    return (rand() > 0.5) ? random_number : -random_number
end
