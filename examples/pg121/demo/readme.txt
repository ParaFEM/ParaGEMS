
  pg121 Demonstrator File Set

  This set of files is for use with the demonstrator application. Program pg121
  reads the "pg121 Input Files" and outputs the "pg121 Output Files".
  The "pg121 ParaView Files" are required to visualize the model using ParaView.

  To load the model into ParaView with pre-assigned settings, open the
  pg121_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.

  To load the model into ParaView with default settings, open the
  pg121_demo.ensi.case file.


  pg121 Input Files

  pg121_demo.dat                 ! control data
  pg121_demo.d                   ! geometry
  pg121_demo.bnd                 ! boundary conditions
  pg121_demo.lds                 ! applied forces


  pg121 Output Files

  pg121_demo.res                 ! summary output
  pg121_demo.ensi.DISPL-000001   ! displacements

  
  pg121 Paraview Files

  pg121_demo.pvsm                ! Paraview settings
  pg121_demo.ensi.case           ! Paraview case file
  pg121_demo.ensi.geo            ! geometry
  pg121_demo.ensi.MATID          ! material numbers
  pg121_demo.ensi.NDBND          ! restrained nodes
  pg121_demo.ensi.NDLDS          ! applied forces
