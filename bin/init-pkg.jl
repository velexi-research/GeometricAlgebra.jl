#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --color=yes --startup-file=no \
           --project=`dirname "${BASH_SOURCE[0]}"` "${BASH_SOURCE[0]}" "$@"
=#
"""
Initialize Julia package.

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
using Logging
using PkgTemplates
using UUIDs

# --- Main program

function main()

    # --- Preparations

    # Define command-line interface
    arg_table = ArgParseSettings()
    @add_arg_table! arg_table begin
        "--overwrite", "-f"
            help = "overwrite pre-existing package files"
            action = :store_true
        "--julia-version", "-j"
            help = "minimum Julia version"
            arg_type = VersionNumber
            default = v"1.4.0"
        "--dest-dir", "-d"
            help = "directory where Julia package will reside"
            default = "."
        "--license", "-l"
            help = "package license"
            default = "ASL"
        "pkg_name"
            help = "package name"
            required = true
    end

    # Parse command-line arguments
    args::Dict = parse_args(ARGS, arg_table)
    pkg_name::String = args["pkg_name"]
    overwrite::Bool = args["overwrite"]
    julia_version::VersionNumber = args["julia-version"]
    dest_dir::String = args["dest-dir"]
    license::String = args["license"]

    # Construct name of temporary directory
    tmp_dir = "tmp.init-pkg." * string(uuid4())

    # --- Initialize Julia package

    # Create package
    create_pkg(pkg_name, tmp_dir, julia_version, license)

    # Move package files to destination directory
    move_pkg(pkg_name, tmp_dir, dest_dir, overwrite)

    # --- Clean up

    clean_up(tmp_dir)
end

# --- Helper Functions

"""
    create_pkg(pkg_name::String, tmp_dir::String,
                   julia_version::VersionNumber, license::String)

Create Julia package.
"""
function create_pkg(pkg_name::String, tmp_dir::String,
                    julia_version::VersionNumber, license::String)

    @info "Creating '$pkg_name' package"

    pkg_template = Template(dir=tmp_dir,
                            julia=julia_version,
                            plugins=[License(name=license),
                                     !SrcDir, !Tests, !Readme,
                                     !Git, !CompatHelper, !TagBot])

    pkg_template(pkg_name)
end

"""
    move_pkg(pkg_name::String, tmp_dir::String, dest_dir::String,
                 overwrite::Bool)

Move Julia package in `tmp_dir/pkg_name` to `dest_dir`.
"""
function move_pkg(pkg_name::String, tmp_dir::String, dest_dir::String,
                  overwrite::Bool=false)

    # Emit progress message
    message = string("Moving package to destination directory ",
                     dest_dir in (".", "..") ? "'$dest_dir'" : dest_dir)
    @info message

    # Preparations
    tmp_pkg_dir::String = joinpath(tmp_dir, pkg_name)
    pkg_contents::Vector = readdir(tmp_pkg_dir, sort=false)

    # Check for items that will be overwritten
    items_will_be_overwritten::Bool = false
    if !overwrite
        for item in pkg_contents
            dest_path = joinpath(dest_dir, item)
            if ispath(dest_path)
                items_will_be_overwritten = true

                # Emit error message
                message = string(
                    "$item already exists in destination directory ",
                    dest_dir in (".", "..") ? "'$dest_dir'" : dest_dir)
                @error message
            end
        end
    end

    # Move package contents
    if !items_will_be_overwritten
        for item in pkg_contents
            mv(joinpath(tmp_pkg_dir, item), joinpath(dest_dir, item),
               force=overwrite)
        end
    end
end

"""
    clean_up(tmp_dir::String)

Remove temporary directories.
"""
function clean_up(tmp_dir::String)
    @info "Cleaning up"
    rm(tmp_dir, force=true, recursive=true)
end

# --- Run main program

main()
