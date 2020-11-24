        program main
        
        use module_param
        
        integer,parameter ::xmax=100
        integer,parameter ::ymax=100

       ! TM
        real(8),dimension(0:xmax,0:ymax):: Ez, Hx, Hy
        real(8) t

        real(8) C_ez, C_ezlx, C_ezly, C_hxly, C_hylx
        real(8) sigma
        
        sigma=0
        t=0
        do n=1,nmax
        ! initial
            do j=0,ymax
                do i=0,xmax
                Ez(i,j)=amp*dsin(2.0d0*pi*f*t)
     &           *exp(-(j-50d0)**2)*exp(-(i-50d0)**2)
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
            
            !output

            write(filename,'(a,i5.5,a)') 
     &        "./output_Ez/",n,".dat"
            open(10,file=filename,FORM='FORMATTED')
            write(10,100,advance="no") 0
            do i=0,xmax-1
                write(10,100,advance="no") i
            end do
            write(10,100,advance="yes") xmax
            
            do j=0,ymax-1
                write(10,100,advance="no") j
                do i=0,xmax-1
                    write(10,101,advance="no") Ez(i,j)
                end do
                write(10,101,advance="yes") Ez(xmax,j)
            end do  
            close(10)
            t=t+dt_e/2
        end do    
100     format(i15.4,1h,100(e15.7,1h,))
101        format(100(e15.7,1h,))
        end program main