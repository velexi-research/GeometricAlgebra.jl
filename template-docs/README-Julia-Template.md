Velexi Template: Julia Package
==============================

___Authors___
Kevin T. Chu `<kevin@velexi.com>`

------------------------------------------------------------------------------

Contents
--------

1. [Overview][#1]

  1.1. [Software Dependencies][#1.1]

  1.2. [Directory Structure][#1.2]

  1.3. [Template Files][#1.3]

2. [Usage][#2]

  2.1. [Setting Up][#2.1]

  2.2. [Running Tests][#2.2]

  2.3. [Cleaning the Project Directory][#1.3]

------------------------------------------------------------------------------

## 1. Overview

This project template is intended to

* streamline the process of setting up a Julia package and

* provide a collection of tools to support Julia software development (e.g.,
  analysis of unit test coverage).

### 1.1. Software Dependencies

#### Base Requirements

* Julia (>=1.4)

#### Julia Packages ####

See `[deps]` section of the `bin/Project.toml` and `test/Project.toml` files.

### 1.2. Directory Structure

    README.md
    README-Template-Usage.md
    README.md.template
    Makefile
    bin/
    extras/
    src/
    test/

* `README.md`: symbolic link to `README-Template-Usage.md`

* `README-Template-Usage.md`: this file

* `README.md.template`: template README file for Julia package

* `Makefile`: Makefile defining a collection of useful commands to maintain
  software (e.g., `test`, `clean`)

* `bin`: directory containing command-line utilities

  * `init-pkg.jl`: utility to initialize Julia package

  * `coverage.jl`: utility to analyze test coverage

* `src`: directory for source code

  * `Example.jl.template` is an example Julia module.

* `test`: directory for test code

  * `runtests.jl` is an example test setup file.

  * `Example_tests.jl.template` contains example unit tests.

* `extras`: directory containing optional content (e.g., configuration files,
  example code snippets, etc.)

### 1.3. Template Files

Template files and directories are indicated by the 'template' suffix. These
files and directories are intended to simplify the set up of package. When
appropriate, they should be renamed (with the 'template' suffix removed).

------------------------------------------------------------------------------

## 2. Usage

### 2.1. Setting Up

* (OPTIONAL) Configure `direnv`.

  * `direnv`. Copy `extras/envrc.template` to the top-level package directory
    and rename it to `.envrc`.

* From the top-level package directory, use `make setup` to add the Julia
  packages required to set up the development environment.

  ```shell
  $ cd PKG_DIR
  $ make setup
  ```

* Initialize Julia package.

  * Use the `init-pkg.jl` utility to initialize the package. In the following
    command, replace `PKG_NAME` with the name of your package.

    ```shell
    $ bin/init-pkg.jl PKG_NAME
    ```

    Note that `init-pkg.jl` supports several command-line options. To see the
    full list of options, use `init-pkg.jl --help`.

  * Update template files in the `src` and `test` directories to be consistent
    with `PKG_NAME`.

    * Replace all references to the `Example` module to `PKG_NAME` in
      filenames and source code.

* Update documentation

  * Replace `README.md` with `README.md.template`.

  * Replace any template fields in the `LICENSE` file with the appropriate
    project-specific information.

### 2.2. Running Tests

* Use `make` to run the unit tests and display coverage information.

  ```shell
  $ make
  ```

  Note that `test` is the default Make target, so `make test` will achieve the
  same result.

### 2.3. Cleaning the Package Directory

* Use `make clean` to automatically remove files and directories generated
  during setup or testing (e.g., temporary directories, coverage files).

  ```shell
  $ make clean
  ```

------------------------------------------------------------------------------

[-----------------------------INTERNAL LINKS-----------------------------]: #

[#1]: #1-overview
[#1.1]: #11-software-dependencies
[#1.2]: #12-directory-structure
[#1.3]: #13-template-files

[#2]: #2-usage
[#2.1]: #21-setting-up
[#2.2]: #22-running-tests
[#2.3]: #23-cleaning-the-package-directory

[#3]: #3-references
