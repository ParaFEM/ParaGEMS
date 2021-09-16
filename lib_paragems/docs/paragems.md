ersion: 2020.09.1
project: ParaGEMS
project_website: https://bitbucket.org/pieterboom/paragems/src/master/
date: January 26, 2021
author: Pieter Boom
author_description: PDRA in the Mechanics and Physics of Solids (MaPoS) research group at the University of Manchester Department of Mechanical, Aerospace and Civil Engineering (MACE)
email: pieter.boom@manchester.ac.uk
src_dir: ../src/
exclude_dir: ../src/libraries/
             ../src/modules/deprecated/
include: ../src/modules/common
output_dir: ./doc_ford
page_dir: ./doc_pages
docmark: >
docmark_alt: @@@
predocmark: ->
predocmark_alt: @@@@@@
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
license: bsd
extra_filetypes: inc !
extra_mods: petsc:https://www.mcs.anl.gov/petsc/
extra_vartypes: PETSC_DECIDE
                Mat
                Vec
                IS
                PetscErrorCode
proc_internals: true
sort: type-alpha

A parallel math library for discrete exterior calculus, along with miniApps for various applications in science and engineering.

