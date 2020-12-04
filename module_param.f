
      module module_param
      implicit none

      integer,parameter ::xmax=400
      integer,parameter ::ymax=400
      integer,parameter ::nmax=1000
      integer,parameter ::sout=100
      real(8), parameter :: pi=4.0d0*datan(1.0d0)
      real(8), parameter :: bol=1.38064852d-23 ![J/K]
      real(8), parameter :: ee=1.6021767208d-19 ![As]

      real(8), parameter :: f=170d9
      real(8), parameter :: sigma=0.0d0
      real(8), parameter :: mu0=1.25663706d-6
      real(8), parameter :: eps0=8.85418782d-12
      real(8), parameter :: cc=2.99792458d8
      real(8), parameter :: amp =1.0d6

      real(8), parameter :: dx=cc/f/10d0
      real(8), parameter :: dt = 0.5d0 * 1/sqrt(2d0/dx**2)/cc
c        dt = 1d0/f/divt

      real(8), parameter :: C_ez
     &                  =(1-sigma*dt/(2*eps0))/(1+sigma*dt/(2*eps0))
      real(8), parameter :: C_ezlx
     &                  =(dt/eps0)/((1+sigma*dt/(2*eps0))*dx)
      real(8), parameter :: C_ezly
     &                  =(dt/eps0)/((1+sigma*dt/(2*eps0))*dx) !dx=dy

      real(8), parameter :: C_hxly=dt/(mu0*dx) !dy=dx
      real(8), parameter :: C_hylx=dt/(mu0*dx)
      integer :: i, j, k, n,m

      end module module_param
