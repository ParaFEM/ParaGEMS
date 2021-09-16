
  PG124 (adapted from P124) Demonstrator File Set

  This set of files is for use with the demonstrator application. Program pg124
  reads the "PG124 Input Files" and outputs the "PG124 Output Files".
  The "PG124 ParaView Files" are required to visualize the model using ParaView.

  To load the model into ParaView with pre-assigned settings, open the
  pg124_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the
  pg124_demo.ensi.case file.


  PG124 Input Files

  pg124_demo.dat                 ! control data
  pg124_demo.d                   ! geometry
  pg124_demo.bnd                 ! boundary conditions


  PG124 Output Files

  pg124_demo.res                 ! summary output
  pg124_demo.ensi.NDTTR-******   ! nodal temperature


  PG124 Paraview Files

  pg124_demo.pvsm                ! Paraview settings
  pg124_demo.ensi.case           ! Paraview case file
  pg124_demo.ensi.geo            ! geometry
  pg124_demo.ensi.MATID          ! material numbers
  pg124_demo.ensi.NDBND          ! restrained nodes
