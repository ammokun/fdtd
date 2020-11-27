        program main
        
        use module_param
        
        integer,parameter ::xmax=200
        integer,parameter ::ymax=200

       ! TM
        real(8),dimension(xmax+1,ymax+1):: Ez, Hx, Hy
        real(8),dimension(1:4,ymax+1) :: Ezx
        real(8),dimension(1:4,xmax+1) :: Ezy
        real(8) t dt

        real(8) C_ez, C_ezlx, C_ezly, C_hxly, C_hylx
        real(8) sigma

        integer :: t_1, t_2 ,count_per_s
        real(8) :: calc_time
        
        sigma=0
        t=0
        dt = 0.99/(cc*sqrt(2/dx**2))

        do j=1,ymax
            do i=1,xmax
            Ez(i,j)=0.0d0
            Hx(i,j)=0.0d0
            Hy(i,j)=0.0d0
            end do
        end do
        do n=1,nmax

        call system_clock(t_1, count_per_s)      
        ! initial
            do j=1,ymax
                do i=1,xmax
                Ez(i,j)=amp*dsin(2.0d0*pi*f*t)
c     &           *exp(-((j-ymax/2)*dx)**2
c     &          *exp(-((i-xmax/2)*dx)**2)
     &           *exp(-((j-ymax/2)*dx)**2 /(2.0d0*(1.0d-4*fb/f)**2))
     &          *exp(-((i-xmax/2)*dx)**2 /(2.0d0*(1.0d-4*fb/f)**2))
                Hx(i,j)=Ez(i,j)/377
                end do
            end do
c            Ez(xmax/2, ymax/2)=amp*dsin(2.0d0*pi*f*t)
c            Hx(xmax/2,ymax/2)=Ez(xmax/2,ymax/2)/377
            ! E field
            do j=2,ymax
                do i=2,xmax
                    
                C_ez=(1-sigma*dt/2/eps0)/(1+sigma*dt/2/eps0)
                C_ezlx=dt/eps0/(1+sigma*dt/2/eps0)*1/dx
                C_ezly=dt/eps0/(1+sigma*dt/2/eps0)*1/dx !dx=dy
                    
                Ez(i,j)=C_ez*Ez(i,j)
     &              +C_ezlx*(Hy(i,j)-Hy(i-1,j))
     &              -C_ezly*(Hx(i,j)-Hx(i,j-1))
                end do
            end do

            !Absorptin
            do j=1,ymax

                Ez(1,j)=Ezx(2,j)
     &          +(cc*dt-dx)/(cc*dt+dx)*(Ez(2,j)-Ezx(1,j))

                Ez(xmax,j)=Ezx(3,j)
     &          +(cc*dt-dx)/(cc*dt+dx)*(Ez(xmax-1,j)-Ezx(4,j))
            end do

            do i=1,xmax

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

            do j=1,ymax
                do i=2,xmax
                    C_hxly=dt/(mu0*dx) !dy=dx
                    C_hylx=dt/(mu0*dx)

                    Hx(i,j)=Hx(i,j)-C_hxly*(Ez(i,j+1)-Ez(i,j))

                end do
            end do  
            do j=2,ymax
                do i=1,xmax
                    C_hxly=dt/(mu0*dx) !dy=dx
                    C_hylx=dt/(mu0*dx)

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
            t=t+dt/2
        end do    
100     format(i15.4,i15.4,1h,100(e15.7,1h,))
101        format(100(e15.7,1h,))
        end program main