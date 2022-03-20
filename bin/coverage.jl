#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --color=yes --startup-file=no \
           --project=`dirname "${BASH_SOURCE[0]}"` "${BASH_SOURCE[0]}" "$@"
=#
"""
Generate coverage analysis.

------------------------------------------------------------------------------
Copyright (c) 2020-2021 Velexi Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
------------------------------------------------------------------------------
"""
# --- Imports

using ArgParse
using Coverage
using Logging
using Printf

# --- Main program

function main()

    # --- Preparations

    # Define command-line interface
    arg_table = ArgParseSettings()
    @add_arg_table! arg_table begin
        "--keep-cov-files", "-k"
        help = "retain *.cov files"
        action = :store_true
        "--pkg-dir", "-d"
        help = "package directory"
        default = "."
        "--verbose", "-v"
        help = "enable verbose mode"
        action = :store_true
    end

    # Parse command-line arguments
    args::Dict = parse_args(ARGS, arg_table)
    keep_cov_files::Bool = args["keep-cov-files"]
    pkg_dir::String = args["pkg-dir"]
    verbose::Bool = args["verbose"]

    # Set log level
    if !verbose
        disable_logging(Logging.Info)
    end

    # Construct paths to src and test directories
    src_dir = joinpath(pkg_dir, "src")
    test_dir = joinpath(pkg_dir, "test")

    # --- Analyze code coverage and display results

    coverage = analyze_coverage(src_dir::String, test_dir::String)
    display_results(coverage)

    return nothing
end

# --- Functions

"""
    analyze_coverage(src_dir::String, tmp_dir::String)

Analyze test coverage.
"""
function analyze_coverage(src_dir::String, test_dir::String)
    # Process '*.cov' files
    coverage = process_folder(src_dir)

    # Process '*.info' files
    coverage = merge_coverage_counts(
        coverage,
        filter!(
            let prefixes = (src_dir, "")
                c -> any(p -> startswith(c.filename, p), prefixes)
            end,
            LCOV.readfolder(test_dir),
        ),
    )

    return coverage
end

"""
    display_results(coverage::Array)

Display coverage results.
"""
function display_results(coverage::Array)

    # Line formats
    header_line_format = "%-35s %15s %10s %10s\n"
    results_line_format = "%-35s %15d %10d %10s\n"
    horizontal_rule = "-"^79

    # Print header line
    println(horizontal_rule)
    print_formatted(header_line_format, "File", "Lines of Code", "Missed", "Coverage")
    println(horizontal_rule)

    # Initialize line counters
    total_lines_of_code = 0
    total_covered_lines_of_code = 0

    # Print coverage for individual files
    for file_coverage in coverage
        filename = file_coverage.filename
        filename = filename[(findlast("src/", filename)[1] + 4):end]

        covered_lines_of_code, lines_of_code = get_summary(
            process_file(file_coverage.filename)
        )
        missed_lines_of_code = lines_of_code - covered_lines_of_code
        coverage_pct = 100 * covered_lines_of_code / lines_of_code
        coverage_pct_str = isnan(coverage_pct) ? "N/A" : @sprintf "%9.1f%%" coverage_pct

        print_formatted(
            results_line_format,
            filename,
            lines_of_code,
            missed_lines_of_code,
            coverage_pct_str,
        )

        # Increment line counteres
        total_lines_of_code += lines_of_code
        total_covered_lines_of_code += covered_lines_of_code
    end

    # Print coverage summary
    total_missed_lines_of_code = total_lines_of_code - total_covered_lines_of_code
    coverage_pct = 100 * total_covered_lines_of_code / total_lines_of_code
    coverage_pct_str = isnan(coverage_pct) ? "N/A" : @sprintf "%9.1f%%" coverage_pct

    println(horizontal_rule)
    print_formatted(
        results_line_format,
        "TOTAL",
        total_lines_of_code,
        total_missed_lines_of_code,
        coverage_pct_str,
    )

    # TODO: add count of tests passed, skipped, failed.
    # TODO: add test runtime

    return nothing
end

"""
    print_formatted(fmt, args...)

Print formatted text.
"""
print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))

# --- Run main program

main()
