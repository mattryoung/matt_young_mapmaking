MODULE mbound
  USE intreal_types
  use bound
  !c Min Max position values and real-integer conversion params
  !c
  !c            iwrap          flag=1 for periodic boundary conditions
  !c            ibmax          max integer value = 2**nbitwrd-1
  !c            xmin(3)        minimum bounding box x,y,z values
  !c            xmax(3)        maximum bounding box x,y,x values
  !c            dx(3)          units of integers
  !c            dx_1(3)        1/dx
  !c            rmin(3)        array of real minimum r(3,np)
  !c            rmax(3)        array of real maximum r(3,np)
  !c
END MODULE mbound
