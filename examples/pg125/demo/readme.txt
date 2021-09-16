
  PG125 (adapted from P125) Demonstrator File Set

  This set of files is for use with the demonstrator application. Program pg125
  reads the "PG125 Input Files" and outputs the "PG125 Output Files".
  The "PG125 ParaView Files" are required to visualize the model using ParaView.

  To load the model into ParaView with pre-assigned settings, open the
  pg125_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.

  To load the model into ParaView with default settings, open the
  pg125_demo.ensi.case file.


  PG125 Input Files

  pg125_demo.dat                 ! control data
  pg125_demo.d                   ! geometry
  pg125_demo.bnd                 ! boundary conditions


  PG125 Output Files

  pg125_demo.res                 ! summary output
  pg125_demo.ensi.NDPRE-******   ! nodal pressure


  PG125 Paraview Files

  pg125_demo.pvsm                ! Paraview settings
  pg125_demo.ensi.case           ! Paraview case file
  pg125_demo.ensi.geo            ! geometry
  pg125_demo.ensi.MATID          ! material numbers
  pg125_demo.ensi.NDBND          ! restrained nodes
