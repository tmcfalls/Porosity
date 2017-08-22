
#define PI 3.14159


struct data_str {
	double Koz;						//Kozeny Carman
	double Temp;					//Linear temperature
	double Ny;						// Niyama Criterion
	double pore_pres;				// Pore presure according to Young-Laplace equation
	double pore_density;			// Pore density
	double pore_frac;				// Pore Fraction
	double Coeff_diff_D1;			// Mass Diffusion Coefficient 
	double Ke;						// Equilibrium Coefficient 
	double Ap;						// Pore Liquid interface area
	double Clp;						// Equilibrium Concentration according to sievert's law
	double Cl;						// Gas concentration in the liquid
	double Coolrate;				// Cooling rate
	double G;						// Thermal Gradient
	double R;						// Solidfication front velocity
	double dmdt;					// Rate of change of the mass of the pore
	double pore_mass;				// Mass of the pore
	double liq_pres;				// Liquid pressure
	double dpdt;					// Rate of change of the liquid pressure
	double pore_rad;				// Radius of the pore
	double drdt;					// Rate of change of the pore mass
	double solid_frac;				// Calcuated solid fraction
	double init_pore_rad;			// Initial pore radius
	double Sec_dend;				// Secondary dendrite Arm Spacing
	double num_den;					// Pore number density
	double init_pore_den;			// Initial pore density
	double atm_pres;				// Atmospheric Pressure
	double gpl;						// Final Pore Fraction
};


struct properties {
	double gas_const;				// Hydrogen Gas Constant (J/kgK)
	double surface_ten;				// Surface Tension (N/m)
	double liq_temp;				// Liquidus Temp (K)
	double solid_temp;				// Solidus Temp (K)
	double h_part_co;				// Hydrogen Partition Coefficient 
	double act_co;					// Activity Coefficient 
	double liq_dyn_visc;			// Liquid Dyanamic Viscosity at Liquidus Temp (Pa.s)
	double liq_density;				// Liquid Density at Liquidus Temp (kg/m3)
	double solid_density;			// Solid Density at Solidus Temp (kg/m3)
	double solid_shrink;			// Solidification Shrinkage
	double const_equil_a1;			// Constant in Equilibrium Coefficient Equation for Hydrogen a1
	double const_equil_b1;			// Constant in Equilibrium Coefficient Equation for Hydrogen b1
	double const_dif_Do;			// Constant in Diffusion Coefficient Equation for Hydrogen Do
	double const_dif_a2;			// Constant in Diffusion Coefficient Equation for Hyd
	double const_sec_den;			// Constant in secondary dendrite arm spacing
	double Co;						// Initial Hydrogen Conentration
	
};

struct simulation {
	double time_start;				// Start Time
	double end_time;				// Time at which solidfication ends
	double time;					// Current time
	double time_flag;				// flag time
	double gpl_flag;				// flag for gpl interations 
	double dt;						// Time Step
	double Nmax;					// Number of iterations
	double gpl_min;					// Min final pore fraction guess								
	double gpl_max;					// Max final pore fraction guess
	double nuc_flag;				// Nucleation flag
	double imp_fac;					// Impingement factor
	double Tol;						// Double Tolearnce 
	double Nuke_val;				// Nucleation Standard
	double Nuke_check;				// Calculated Nucleation value
	double Nuke_flag;				// Nucleation check flag
	int i;							// imlicit iteration value
};
