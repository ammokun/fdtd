
      module module_param
      implicit none


!-Defined_Parameter-!

      real(8), parameter :: pi=4.0d0*datan(1.0d0)
      real(8), parameter :: bol=1.38064852d-23 ![J/K]
      real(8), parameter :: ee=1.6021767208d-19 ![As]
      real(8), parameter :: me=9.10938356d-31 ![kg]
      real(8), parameter :: eps0=8.85418782d-12   
      real(8), parameter :: mu0=1.25663706d-6      
      real(8), parameter :: cc=2.99792458d8      
      real(8),parameter :: planck=6.62607d-34 ![Js]
      real(8),parameter :: omega_n2=2358.6d2*cc ![rad/s]
      real(8),parameter :: kai_n2=14.32d2/2358.6d2
      
!-Incident_Wave_Parameter-!

      real(8), parameter :: f=170.0d9 !wave_frequency
      real(8), parameter :: fb=170.0d9 !baced wave_frequency for fixing dx
      real(8), parameter :: amp=1.6d6*dsqrt(2.0d0) !wave_amplitude
      real(8), parameter :: lam=cc/f

!-Calculation_Condition-!

      integer :: i,j,k,l,n,m,r,z,v,w,a
      integer :: num,i_dum,cont_num,cont_num2,point
      integer, parameter :: nmax=1000 !calcuration_number !about 20e-6[s]
      integer, parameter :: divx=int(100*fb/f) !cell_division_number
      integer, parameter :: divp=10 !cell_division_number
      integer, parameter :: divt=divx*1.1 !time_division_number
      integer, parameter :: imax=int(300*fb/f)+1 !cell_number

      integer, parameter :: sout=10 !output_interval
      integer, parameter :: lim_out=5001 !output_limit(n/sout)
      integer, parameter :: datamax=1000 !number_of_input_data

      real(8), parameter :: dx=lam/dble(divx) !dt_cell_length
      real(8), parameter :: dt_g=1.0d0/f !del_time_for_particles
      real(8), parameter :: dt_p=1.0d0/f !del_time_for_particles
      real(8), parameter :: dt_e=1.0d0/f/dble(divt) !del_time_for_waves

!-Numerical_Model-!

      integer, parameter :: mirror=1 !set/unset_mirror_at_right_boundary
      integer, parameter :: advection=1 !set/unset_advection_for_charged_particle
      integer, parameter :: photoion=0 !set/unset_photo-ionization_process
      integer, parameter :: photodifN2=0 !set/unset_radiation-transfer_process
      integer, parameter :: photodifO2=0 !set/unset_radiation-transfer_process
      integer, parameter :: photoexcitation=1      

      integer, parameter :: chem=1 !set/unset_rate_equation_solver
      integer, parameter :: bolsig=1 !transport_coefficient_obtaied_by_experimental-data/bolsig
      integer, parameter :: eueu=1 !set/unset_euler_equation_solver
      integer, parameter :: poisson=0 !set/unset_poisson_equation_solver
      integer, parameter :: implicit_diffusion=1 !diffusion_equation_solver_using_implicit_method
      integer, parameter :: implicit_energy=1 !energy_equation_solved_using_implicit_method
      integer, parameter :: energy_diffusion=1 !set/unset_electron_energy_diffusion

!-Neutrals_Parameter-!

      real(8), parameter :: gas0=8.314459848d0 !gas_constant[J/(K*mol)]
      real(8), parameter :: mN2=28.014d-3 !molar_mass_for_N2[kg/mol]
      real(8), parameter :: mO2=32.00d-3 !molar_mass_for_O2[kg/mol]
      real(8), parameter :: NA=6.02214086d23 !Avogadro_constant[mol-1]
      real(8), parameter :: rN2=0.80d0
      real(8), parameter :: rO2=0.20d0
      real(8), parameter :: gamm0=1.40d0
      real(8), parameter :: press0=760.0d0 !ambient_pressure[Torr]
      real(8), parameter :: temg0=300.0d0 !gas_temperature[K]
      real(8), parameter :: temv0=temg0 !vibrational_temperature[K]
      real(8), parameter :: rhog0=press0*133.322d0/temg0
     &                           /(gas0/mN2*rN2+gas0/mO2*rO2) 
      real(8), parameter :: ng0=rhog0/(mN2*rN2+mO2*rO2)*NA      
      real(8), parameter :: velg0=0.0d0
      real(8), parameter :: tauvt=1.0d-6 !vibrational_relaxation_time

!-Electrons_Parameter-!

      real(8), parameter :: teme0=temg0 !electron_temperature[K]
      real(8),parameter :: col0=5.3d9*press0 !collision_frequency[s-1]

!-Reactions_Parameter-!

      integer, parameter :: zmaxi=9 !number_of_positive_ion_types
      integer, parameter :: zmaxn=5 !number_of_negative_ion_types
      integer, parameter :: zmaxg=40 !number_of_nurtral_particle_types
      integer, parameter :: rmax=1555 !number_of_reactions_for_all
      integer, parameter :: rbmax=684 !number_of_reactions_used_by_bolsig
      integer, parameter :: vt_max=120 !number_of_VTreactions
      integer, parameter :: vv_max=380 !number_of_VVreactions
      integer, parameter :: vd_max=63 !number_of_VDreactions
      integer, parameter :: rpmax=8 !number of Photon chemical reactions

      real(8), parameter :: ene_vibN2=planck*omega_n2/bol ![K]
      real(8), parameter :: len_vib=0.16d0 !characteristic parameter [1E-10 m]      
      
!-Photo_Ionization_Parameter-!

      real(8),parameter :: del_x  = dx*1.0d2    ![cm]
      real(8),parameter :: KAIMIN = 0.035d0     ![Torr^-1cm^-1]
      real(8),parameter :: KAIMAX = 2.000d0     ![Torr^-1cm^-1]  

!-Others-!

      real(8) :: dum,Radius,func_G,PART_Sph
      character(50) :: filename,filename1,filename2,filename3,filename4


      end module module_param
