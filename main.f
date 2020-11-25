        program main
        
        use module_param
        
        integer,parameter ::xmax=200
        integer,parameter ::ymax=200

       ! TM
        real(8),dimension(0:xmax,0:ymax):: Ez, Hx, Hy
        real(8) t dt

        real(8) C_ez, C_ezlx, C_ezly, C_hxly, C_hylx
        real(8) sigma

        integer :: t_1, t_2 ,count_per_s
        real(8) :: calc_time
        
        sigma=0
        t=0
        dt = 0.99/(cc*sqrt(2/dx**2))
        do n=1,nmax

        call system_clock(t_1, count_per_s)      
        ! initial
            do j=0,ymax
                do i=0,xmax
                Ez(i,j)=amp*dsin(2.0d0*pi*f*t)
     &           *exp(-((j-ymax/2)*dx)**2 /(2.0d0*(1.0d-4*fb/f)**2))
     &          *exp(-((i-xmax/2)*dx)**2 /(2.0d0*(1.0d-4*fb/f)**2))
                end do
            end do

            
            ! E field
            do j=1,ymax-1
                do i=1,xmax-1
                    
                C_ez=(1-sigma*dt_e/2/eps0)/(1+sigma*dt_e/2/eps0)
                C_ezlx=dt_e/eps0/(1+sigma*dt_e/2/eps0)*1/dx
                C_ezly=dt_e/eps0/(1+sigma*dt_e/2/eps0)*1/dx !dx=dy
                    
                Ez(i,j)=C_ez*Ez(i,j)
     &              +C_ezlx*(Hy(i,j)-Hy(i-1,j))
     &              -C_ezly*(Hx(i,j)-Hx(i,j-1))
                end do
            end do

            !Absorptin
            do j=1,ymax-1

                Ez(1,j)=Ez(2,j)
     &          +(cc*dt_e-dx)/(cc*dt_e+dx)*(Ez(2,j)-Ez(1,j))

                Ez(xmax-1,j)=Ez(xmax-2,j)
     &          +(cc*dt_e-dx)/(cc*dt_e+dx)*(Ez(xmax-2,j)-Ez(xmax-1,j))
            end do

            do i=1,xmax-1

                Ez(i,1)=Ez(i,2)
     &          +(cc*dt_e-dx)/(cc*dt_e+dx)*(Ez(i,2)-Ez(i,1))
                Ez(i,ymax-1)=Ez(i,ymax-2)
     &          +(cc*dt_e-dx)/(cc*dt_e+dx)*(Ez(i,ymax-2)-Ez(i,ymax-1))

            end do
            t=t+dt_e/2
            !H field

            do j=1,ymax-1
                do i=1,xmax-1
                    C_hxly=dt_e/(mu0*dx) !dy=dx
                    C_hylx=dt_e/(mu0*dx)

                    Hx(i,j)=Hx(i,j)-C_hxly*(Ez(i,j+1)-Ez(i,j))
                    Hy(i,j)=Hy(i,j)+C_hylx*(Ez(i+1,j)-Ez(i,j))
                    
                end do
            end do  

            call system_clock(t_2)
            calc_time=real(t_2-t_1)/count_per_s

            !output
            if(mod(n,10)==0) then
                write(*,*)"n" ,n
                write(*,*) "time", calc_time 
            end if

            if(mod(n,sout)==0) then
            write(*,*)"n, output" ,n,n/sout
            write(*,*) "time", calc_time 
            write(filename,'(a,i5.5,a)') 
     &        "./output_Ez/",n/sout,".dat"
            open(10,file=filename,FORM='FORMATTED')
            !write(10,100,advance="no") 0
            !do i=0,xmax-1
            !    write(10,100,advance="no") i
            !end do
            !write(10,100,advance="yes") xmax
            
            do i=0,xmax
            !    write(10,100,advance="no") j
                do j=0,ymax
                    write(10,100) i,j,Ez(i,j)
                end do
                write(10,*)" " 
            end do  
            close(10)
            end if
            t=t+dt_e/2
        end do    
100     format(i15.4,i15.4,1h,100(e15.7,1h,))
101        format(100(e15.7,1h,))
        end program main