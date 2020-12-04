        program main
        use module_param
        use omp_lib
        implicit none


       ! TM
        real(8),dimension(1:xmax,1:ymax):: Ez, Hx, Hy
        real(8),dimension(1:xmax,1:ymax):: Ez_o, Hx_o, Hy_o
        real(8),dimension(1:4,1:ymax) :: Ezx
        real(8),dimension(1:4,1:xmax) :: Ezy
        real(8) ::t

c        real(8) :: C_ez, C_ezlx, C_ezly, C_hxly, C_hylx

        integer :: t_1, t_2 ,count_per_s,s
        double precision :: t_start=0d0,t_end=0d0
        real(8) :: calc_time
        character(50) :: filename
        
        !$ call omp_num_threads(4) 

        t=0d0

        
        do j=1,ymax
            do i=1,xmax
            Ez(i,j)=0.0d0
            Hx(i,j)=0.0d0
            Hy(i,j)=0.0d0
            Ez_o(i,j)=0.0d0
            Hx_o(i,j)=0.0d0
            Hy_o(i,j)=0.0d0
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



        write(*,*) "C_ez=", C_ez
        write(*,*) "C_ezlx",  C_ezlx
        write(*,*) " C_ezly",  C_ezly
        write(*,*) "C_hxly", C_hxly
        write(*,*) "C_hylx", C_hylx
        write(*,*) "dt, dx, 1/f", dt, dx, 1d0/f

        t_start=omp_get_wtime()
        do n=1,nmax       
            call system_clock(t_1,count_per_s)
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

            call Ez_field(Ez,Hx,Hy,Ez_o)
            
            !$OMP parallel do private(i,j)
            do j=1,ymax
                do i=1,xmax
                    Ez(i,j)=Ez_o(i,j)
                end do
            end do
            !$OMP end parallel do
                !Absorption


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

            call Hx_field(Ez,Hx,Hx_o)
            call Hy_field(Ez,Hy,Hy_o)

            !$OMP parallel do private(i,j)
            do j=1,ymax
                do i=1,xmax
                    Hx(i,j)=Hx_o(i,j)
                    Hy(i,j)=Hy_o(i,j)
                end do
            end do
            !$OMP end parallel do
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

        t_end=omp_get_wtime()
        write(*,*) "Total time=",t_end-t_start

100     format(i15.4,i15.4,1h,100(e15.7,1h,))
101     format(100(e15.7,1h,))
        
        end program main
 
        subroutine Ez_field(Ez_i, Hx_i, Hy_i,Ez_o)
        use module_param
        implicit none

        real(8), intent(in) :: Ez_i(1:xmax,1:ymax) 
        real(8), intent(in) :: Hx_i(1:xmax,1:ymax)
        real(8), intent(in) :: Hy_i(1:xmax,1:ymax)
        real(8), intent(out) :: Ez_o(1:xmax,1:ymax)  
            
        !$OMP parallel do private(i,j)
        do j=2,ymax
  
            do i=2,xmax
                            
            Ez_o(i,j)=C_ez*Ez_i(i,j)
     &             +C_ezlx*(Hy_i(i,j)-Hy_i(i-1,j))
     &              -C_ezly*(Hx_i(i,j)-Hx_i(i,j-1))
            end do
         end do
        !$OMP end parallel do
        end subroutine Ez_field

        subroutine Hx_field(Ez_i,Hx_i,Hx_o)
        use module_param
        implicit none

        real(8), intent(in) :: Ez_i(1:xmax,1:ymax) 
        real(8), intent(in) :: Hx_i(1:xmax,1:ymax)
        real(8), intent(out) :: Hx_o(1:xmax,1:ymax)  
        
        !$OMP parallel do private(i,j)
        do j=1,ymax-1
            do i=1,xmax
                Hx_o(i,j)=Hx_i(i,j)-C_hxly*(Ez_i(i,j+1)-Ez_i(i,j))
            end do
        end do  
        !$OMP end parallel do
        end subroutine Hx_field


        subroutine Hy_field(Ez_i,Hy_i,Hy_o)
        use module_param
        implicit none
    
        real(8), intent(in) :: Ez_i(1:xmax,1:ymax) 
        real(8), intent(in) :: Hy_i(1:xmax,1:ymax)
        real(8), intent(out) :: Hy_o(1:xmax,1:ymax)  

        !$OMP parallel do private(i,j)
        do j=1,ymax 
            do i=1,xmax-1

                Hy_o(i,j)=Hy_i(i,j)+C_hylx*(Ez_i(i+1,j)-Ez_i(i,j))
                
            end do
        end do
        !$OMP end parallel do
        end subroutine Hy_field