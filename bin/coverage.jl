#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --color=yes --startup-file=no \
           --project=`dirname "${BASH_SOURCE[0]}"` "${BASH_SOURCE[0]}" "$@"
=#
"""
Generate coverage analysis.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
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
    coverage = merge_coverage_counts(coverage,
        filter!(
            let prefixes = (src_dir, "")
                c -> any(p -> startswith(c.filename, p), prefixes)
            end,
            LCOV.readfolder(test_dir)
        )
    )
end

"""
    display_results(coverage::Array)

Display coverage results.
"""
function display_results(coverage::Array)

    header_line_format = "%-35s %15s %10s %10s\n"
    results_line_format = "%-35s %15d %10d %9.1f%%\n"
    horizontal_rule = "-"^79

    # Print header line
    println(horizontal_rule)
    print_formatted(header_line_format,
                    "File", "Lines of Code", "Missed", "Coverage")
    println(horizontal_rule)

    # Print coverage for individual files
    for file_coverage in coverage
        covered_lines, total_lines =
            get_summary(process_file(file_coverage.filename))

        filename = file_coverage.filename[
            findlast("src/", file_coverage.filename)[1]+4:end]

        print_formatted(results_line_format, filename,
                        total_lines, total_lines - covered_lines,
                        covered_lines / total_lines * 100)
    end

    # Print coverage summary
    covered_lines, total_lines = get_summary(coverage)
    println(horizontal_rule)
    print_formatted(results_line_format, "TOTAL",
                    total_lines, total_lines - covered_lines,
                    covered_lines / total_lines * 100)

    # TODO: add count of tests passed, skipped, failed.
    # TODO: add test runtime
end

"""
    print_formatted(fmt, args...)

Print formatted text.
"""
print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))

# --- Run main program

main()
