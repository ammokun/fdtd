        program main
        
        integer,parameter ::xmax=100
        integer,parameter ::ymax=100
        integer,parameter ::nmax=1000
        integer,parameter ::sout=10
        real(8), parameter :: pi=4.0d0*datan(1.0d0)
        real(8), parameter :: bol=1.38064852d-23 ![J/K]
        real(8), parameter :: ee=1.6021767208d-19 ![As]

        real(8), parameter :: f=170d9
        real(8), parameter :: sigma=0.0d0
        real(8), parameter :: mu0=1.25663706d-6
        real(8), parameter :: eps0=8.85418782d-12
        real(8), parameter :: cc=2.99792458d8
        real(8), parameter :: amp =1.0d6

       ! TM
        real(8),dimension(1:xmax,1:ymax):: Ez, Hx, Hy
        real(8),dimension(1:4,1:ymax) :: Ezx
        real(8),dimension(1:4,1:xmax) :: Ezy
        real(8) :: t dt dx

        real(8) :: C_ez, C_ezlx, C_ezly, C_hxly, C_hylx

        integer :: i, j, k, n,m
        integer :: t_1, t_2 ,count_per_s,s
        real(8) :: calc_time
        character(50) :: filename

        

        t=0d0
        dx=cc/f/10d0
        dt = 0.5d0 * 1/sqrt(2d0/dx**2)/cc
c        dt = 1d0/f/divt
        
        do j=1,ymax
            do i=1,xmax
            Ez(i,j)=0.0d0
            Hx(i,j)=0.0d0
            Hy(i,j)=0.0d0
            end do
        end do

        do j=1,ymax

            Ezx(1,j)=0.0d0
            Ezx(2,j)=0.0d0
            Ezx(3,j)=0.0d0
            Ezx(4,j)=0.0d0

        end do
        do i=1,xmax

            Ezy(1,i)=0.0d0
            Ezy(2,i)=0.0d0
            Ezy(3,i)=0.0d0
            Ezy(4,i)=0.0d0

        end do

        C_ez=(1-sigma*dt/(2*eps0))/(1+sigma*dt/(2*eps0))
        C_ezlx=(dt/eps0)/((1+sigma*dt/(2*eps0))*dx)
        C_ezly=(dt/eps0)/((1+sigma*dt/(2*eps0))*dx) !dx=dy

        C_hxly=dt/(mu0*dx) !dy=dx
        C_hylx=dt/(mu0*dx)

        write(*,*) "C_ez=", C_ez
        write(*,*) "C_ezlx",  C_ezlx
        write(*,*) " C_ezly",  C_ezly
        write(*,*) "C_hxly", C_hxly
        write(*,*) "C_hylx", C_hylx
        write(*,*) "dt, dx, 1/f", dt, dx, 1d0/f


        do n=1,nmax

        call system_clock(t_1, count_per_s)
c            do m=1,divt !FDTD loop      
        ! initial
c            Ez(xmax/2,ymax/2)=amp*exp(-(4/dt)**2*(t-2*dt)**2)
            Ez(xmax/2,ymax/2)=amp*dsin(2.0d0*pi*f*t)
            Hy(xmax/2,ymax/2)=-Ez(xmax/2,ymax/2)/377d0
            Hx(xmax/2,ymax/2)=Ez(xmax/2,ymax/2)/377d0
    
c            do j=1,ymax
c                do i=1,xmax
c                 Ez(2,j)=amp*dsin(2.0d0*pi*f*t) !0.5d0*(1.0d0-dcos(pi*t/3d0/dt))
                
c                Ez(i,j)=amp*dsin(2.0d0*pi*f*t)*exp(-0.01d0/f*t**2)
c     &           *exp(-((j-ymax/2)*dx)**2)
c     &           *exp(-((i-xmax/2)*dx)**2)
c     &           *exp(-((j-ymax/2)*dx)**2 /(2.0d0*(1.0d-4*fb/f)**2))
c     &          *exp(-((i-xmax/2)*dx)**2 /(2.0d0*(1.0d-4*fb/f)**2))
c                Hy(2,j)=-Ez(2,j)/377
c                Hx(2,j)=Ez(2,j)/377
c                end do
c            end do

                ! E field


            do j=2,ymax
                do i=2,xmax
                                        
                Ez(i,j)=C_ez*Ez(i,j)
     &              +C_ezlx*(Hy(i,j)-Hy(i-1,j))
     &              -C_ezly*(Hx(i,j)-Hx(i,j-1))
                end do
            end do


            !Absorptin

            do j=2,ymax-1

                Ez(1,j)=Ezx(2,j)
     &          +(cc*dt-dx)/(cc*dt+dx)*(Ez(2,j)-Ezx(1,j))

                Ez(xmax,j)=Ezx(3,j)
     &          +(cc*dt-dx)/(cc*dt+dx)*(Ez(xmax-1,j)-Ezx(4,j))
            end do

            do i=2,xmax-1

                Ez(i,1)=Ezy(2,i)
     &          +(cc*dt-dx)/(cc*dt+dx)*(Ez(i,2)-Ezy(1,i))
                Ez(i,ymax)=Ezy(3,i)
     &          +(cc*dt-dx)/(cc*dt+dx)*(Ez(i,ymax-1)-Ezy(4,i))

            end do

            do j=1,ymax

                Ezx(1,j)=Ez(1,j)
                Ezx(2,j)=Ez(2,j)
                Ezx(3,j)=Ez(xmax-1,j)
                Ezx(4,j)=Ez(xmax,j)

            end do
            do i=1,xmax

                Ezy(1,i)=Ez(i,1)
                Ezy(2,i)=Ez(i,2)
                Ezy(3,i)=Ez(i,ymax-1)
                Ezy(4,i)=Ez(i,ymax)

            end do
            t=t+dt/2
            !H field

            do j=1,ymax-1
                do i=1,xmax
                    Hx(i,j)=Hx(i,j)-C_hxly*(Ez(i,j+1)-Ez(i,j))
                end do
            end do  
            do j=1,ymax
                do i=1,xmax-1

                    Hy(i,j)=Hy(i,j)+C_hylx*(Ez(i+1,j)-Ez(i,j))
                    
                end do
            end do  
            t=t+dt/2

c        end do !end FDTD loop
            call system_clock(t_2)
            calc_time=real(t_2-t_1)/count_per_s

            !output
c            if(mod(n,10)==0) then
c                write(*,*)"n" ,n
c                write(*,*) "time", calc_time 
c            end if

            if(mod(n,sout)==0) then
            write(*,*)"n t output" ,n,t,n/sout
            write(*,*) "time", calc_time 
            write(*,*) "|Ez=|",amp*dsin(2.0d0*pi*t*f)
            write(filename,'(a,i5.5,a)') 
     &        "./output_Ez/",n/sout,".dat"
            open(10,file=filename,FORM='FORMATTED')
            !write(10,100,advance="no") 0
            !do i=0,xmax-1
            !    write(10,100,advance="no") i
            !end do
            !write(10,100,advance="yes") xmax
            
            do i=1,xmax
            !    write(10,100,advance="no") j
                do j=1,ymax
                    write(10,100) i,j,Ez(i,j), Hx(i,j), Hy(i,j)
                end do
                write(10,*)" " 
            end do  
            close(10)
            end if

        end do    
100     format(i15.4,i15.4,1h,100(e15.7,1h,))
101     format(100(e15.7,1h,))
        
        write(*,*) n
        end program main