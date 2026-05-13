---
project: LowMachReact-Hex
summary: Fortran 2008 MPI finite-volume solver for constant-density low-Mach / incompressible flow with passive species, passive sensible enthalpy, Cantera thermodynamics, and future radiation coupling.
author: Khaled Mosharraf Mukut
project_github: https://github.com/kmmukut/LowMachReact-Hex
project_website: https://kmmukut.github.io/LowMachReact-Hex
src_dir: src
output_dir: docs
page_dir: doc_src
display: public
    protected
    private
source: true
incl_src: true
proc_internals: false
graph: true
call_graph: true
used_by_graph: true
inheritance_graph: true
show_proc_parent: true
search: false
sort: permission-alpha
md_extensions: markdown.extensions.extra
    markdown.extensions.codehilite
    markdown.extensions.toc
    markdown.extensions.admonition
extra_filetypes: cpp // c_cpp.CppLexer
doxygen: true
encoding: utf-8
print_creation_date: true
warn: true
---

LowMachReact-Hex is a Fortran 2008 MPI finite-volume solver for constant-density low-Mach / incompressible flow with passive species, passive sensible enthalpy, Cantera thermodynamics, and future radiation coupling.

For hand-written documentation pages, see the pages generated from `doc_src/`.
