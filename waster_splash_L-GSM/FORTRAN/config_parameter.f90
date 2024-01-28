module config_parameter
implicit none
!dim : Dimension of the problem (1, 2 or 3)
public dim
integer(4),parameter :: dim=2

!maxn : Maximum number of particles
!max_interation : Maximum number of interaction pairs
public maxn,max_interaction
integer(4),parameter :: maxn=30000
integer(4),parameter :: max_interaction=30*maxn

!Parameters for the computational geometry,  
!x_maxgeom : Upper limit of allowed x-regime 
!x_mingeom : Lower limit of allowed x-regime 
!y_maxgeom : Upper limit of allowed y-regime 
!y_mingeom : Lower limit of allowed y-regime 
!z_maxgeom : Upper limit of allowed z-regime 
!z_mingeom : Lower limit of allowed z-regime
public x_maxgeom,x_mingeom,y_maxgeom,y_mingeom,z_maxgeom,z_mingeom
real(8),parameter :: x_maxgeom=3e1
real(8),parameter :: x_mingeom=-1.e1
real(8),parameter :: y_maxgeom=2e1
real(8),parameter :: y_mingeom=-1e1
real(8),parameter :: z_maxgeom=1e1
real(8),parameter :: z_mingeom=-1e1
    
!SPH algorithm for particle approximation (pa_sph)
!pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
!         2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
public pa_sph
integer(4),parameter :: pa_sph=2

!Nearest neighbor particle searching (nnps) method
!nnps = 1 : Simplest and direct searching
!       2 : Sorting grid linked list
!       3 : Tree algorithm
!       4 : GSM algorithm
public nnps
integer(4),parameter :: nnps=4

!Flow problems
!=1 with buffle
!=2 ×²»÷ÓÒ±Ú
!=3 water discharge
!=4 Splash
public flowtype
integer(4),parameter :: flowtype=4

!Smoothing length evolution (sle) algorithm
!sle = 0 : Keep unchanged
!      1 : h = fac * (m/rho)^(1/dim)
!      2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
!      3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 
public sle
integer(4),parameter :: sle=0

!Smoothing kernel function 
!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!      2, Gauss kernel   (Gingold and Monaghan 1981) 
!      3, Quintic kernel (Morris 1997)
!	   4, yang
public skf
integer(4),parameter :: skf=1

!Material Strength Model 
!msm = 0, no strength model is used
!      1, perfectly plastic yield model
!      2, Jahnson-Cook yield model 
public msm
integer(4),parameter :: msm=1

!Equation of State 
!eos = 1, Gamma law for ideal gas
!      2, Artificial EOS by Monaghan (1994) 
!      3, Artificial EOS by Morris (1997)
!      4, Mie-Gruneisen Equation of State 
!      5, Tillotson Equation of State
!      6, JWL & Mie-Gruneisen
!      7, Soil (From H.H.Bui')
public eos
integer(4),parameter :: eos=7

!     Switches for different senarios

!     summation_density = .TRUE. : Use density summation model in the code, 
!                         .FALSE.: Use continuiity equation
!     average_velocity = .TRUE. : Monaghan treatment on average velocity,
!                        .FALSE.: No average treatment.
!     config_input = .TRUE. : Load initial configuration data,
!                    .FALSE.: Generate initial configuration.
!     virtual_part1 = .TRUE. : Use type1 vritual particle,
!                    .FALSE.: No use of Monaghan vritual particle.
!     virtual_part2 = .TRUE. : Use type2 vritual particle,
!                    .FALSE.: No use of Libersky vritual particle.
!     vp_input = .TRUE. : Load virtual particle information,
!                .FALSE.: Generate virtual particle information.
!     visc = .true. : Consider viscosity,
!            .false.: No viscosity.
!     ex_force =.true. : Consider external force,
!               .false.: No external force.
!     visc_artificial = .true. : Consider artificial viscosity,
!                       .false.: No considering of artificial viscosity.
!     heat_artificial = .true. : Consider artificial heating,
!                       .false.: No considering of artificial heating.
!     self_gravity = .true. : Considering self_gravity,
!                    .false.: No considering of self_gravity
!     nor_density =  .true. : Density normalization by using CSPM,
!                    .false.: No normalization.
public kernel_correct,summation_density,average_velocity,config_input,virtual_part1,virtual_part2
public vp_input,visc,ex_force,visc_artificial,heat_artificial,self_gravity,nor_density,stress_artificial

logical,parameter :: kernel_correct=.false.
logical,parameter :: summation_density=.false.
logical,parameter :: average_velocity=.true.
logical,parameter :: config_input=.false.
logical,parameter :: virtual_part1=.true.
logical,parameter :: virtual_part2=.false.
logical,parameter :: vp_input=.false.
logical,parameter :: visc=.false.
logical,parameter :: ex_force=.true.
logical,parameter :: visc_artificial=.true.
logical,parameter :: heat_artificial=.false.
logical,parameter :: stress_artificial=.false.
logical,parameter :: self_gravity=.true.
logical,parameter :: nor_density=.false.

!Symmetry of the problem
!nsym = 0 : no symmetry,
!       1 : axis symmetry,
!       2 : center symmetry.     
public nsym
integer(4),parameter :: nsym=0

!     Control parameters for output 
!     int_stat = .true. : Print statistics about SPH particle interactions.
!                         including virtual particle information.
!     print_step: Print Timestep (On Screen)
!     save_step : Save Timestep    (To Disk File)
!	  update_step: used for restart or continue of analysis
!     moni_particle: The particle number for information monitoring.
!     
public int_stat,print_step,save_period,moni_particle,energy_check

logical,parameter :: int_stat=.true.
integer(4),parameter :: print_step = 1000
real(8),parameter :: save_period = 0.1
real(8),parameter :: moni_particle = 100
logical,parameter :: energy_check=.false.

!     save_step set to 1 for instant record of virtual particle changes      
public pi
real(8),parameter :: pi=3.14159265358979323846    

!Simulation cases
!taylor = .true. : carry out taylor bair simulation
!hvisphere = .true. : carry out HVI sphere simulation
!tungstencube = .true. : carry out Tungsten cube simulation
public taylor,hvisphere,tungstencube,TNT,landslide
logical,parameter :: taylor=.false.
logical,parameter :: hvisphere=.false.
logical,parameter :: tungstencube=.true.
logical,parameter :: TNT=.false.
logical,parameter :: landslide=.false.

!epsilon in XSPH is 0.175 (if deemed necessary)

end module

module node4linkedlist
type node
integer(4) nodeid
type(node),pointer :: nextnode
end type
end module