#include<iostream>
#include<sstream>
#include<cmath>

#include "Header.h"

using namespace std;

// Initalize property values
void init_props(properties& props) 
{
	props.gas_const		= 4124;
	props.surface_ten	= 0.8;
	props.liq_temp		= 889;
	props.solid_temp	= 824;
	props.h_part_co		= 0.07;
	props.act_co		= 1.25;
	props.liq_dyn_visc  = 0.00158;
	props.liq_density	= 2429.4;
	props.solid_density = 2572;
	props.solid_shrink	= 0.0587;
	props.const_equil_a1 = 2691.96;
	props.const_equil_b1 = 1.32;
	props.const_dif_Do = 3.8e-6;
	props.const_dif_a2 = 2315;
	props.const_sec_den = 4.09e-5;
	props.Co = 1.68e-7;
	//props.Co = 0;

}

// Initialize data values
void init_dat(data_str& dat)
{
	dat.Ny =50;
	dat.Coolrate = 1; 
	dat.G = 1;
	dat.R = 1;
	dat.atm_pres = 101325;
	dat.liq_pres= dat.atm_pres;
	dat.pore_mass = 0;
	dat.init_pore_rad = 1e-5;
	dat.pore_rad = dat.init_pore_rad;
	dat.num_den = 1e11;
	dat.init_pore_den = 0;
	dat.pore_density = dat.init_pore_den;
	dat.gpl = .0031;
}

// Initialize simulation values
void init_sim(simulation& sim, properties& props, data_str& dat)
{
	sim.time_start = 0; 
	sim.time_flag = 1;
	sim.gpl_flag = 1;
	sim.time = 0;
	sim.dt = 1e-6;
	sim.end_time = (props.liq_temp - props.solid_temp)*dat.Coolrate;
	sim.Nmax = sim.end_time / sim.dt;
	sim.gpl_min = 1;
	sim.gpl_max = 1.5;
	sim.nuc_flag = 1;
	dat.Sec_dend = props.const_sec_den*(pow(dat.Coolrate, (-1 / 3)));
	sim.imp_fac = 3;
	sim.Tol = 1e-6;;
	sim.Nuke_val = (2 * props.surface_ten) / dat.init_pore_rad;
	sim.i = 0;
}


//Calculate solid Fraction
double gscalc(double Temp) 
{
	double gs;
	if (Temp >= 847.08) {
		gs = (-117.4336937804321) + (.2825446567378302)*Temp + (-0.0001692136074217557)*(pow(Temp, 2));
	}
	else if ((Temp < 847.08) && (Temp >= 831.45)) 
	{
		gs = (-3486364.256372994) + (9808.586485436901)*Temp + (-1.940453660798732)*(pow(Temp, 2)) + 
		(-0.01176755701724439)*(pow(Temp, 3)) + (-1.925445708386995e-06)*(pow(Temp, 4)) + 
		(1.922508061469615e-08)*(pow(Temp,5)) + (-9.936062903754908e-12)*(pow(Temp, 6));
	}
	else 
	{
		gs = (-973.4459009099721) + (2.360375373193465)*Temp + (-0.001429367924460962)*(pow(Temp, 2));
	
	}
	return gs;
		
}







// Calculating Pore Mass Change Pre Nucleation
void masscalcpre(data_str& dat, simulation& sim, properties& props)
{
	dat.dmdt = props.liq_density*dat.Coeff_diff_D1*dat.Ap*((dat.Cl - dat.Clp) / dat.pore_rad);
	dat.pore_mass = dat.pore_mass + (dat.dmdt*sim.dt);
}



// Calculating Liquid Presure Change Pre Nucleation
void Liq_prescalcpre(data_str& dat, simulation& sim, properties& props,double gs)
{
	dat.dpdt = -1 * ((props.liq_dyn_visc*dat.Coolrate) / (dat.Koz*pow(dat.Ny, 2)))*
		(((props.solid_shrink)*(1 - dat.gpl - gs)) - (dat.gpl - dat.pore_frac));
	dat.liq_pres = dat.liq_pres + dat.dpdt*sim.dt;

}


// Calculating Pore Radius Change Pre Nucleation
void radcalcpre(data_str& dat, simulation& sim, properties& props)
{
	double blake;
	blake = 1 / ((3 * dat.liq_pres*pow(dat.pore_rad, 2)) + (4 * props.surface_ten*dat.pore_rad));
	//Blake singularity fix
	//blake = 1 / ((3 * dat.liq_pres*pow(dat.pore_rad, 2)) + (4 * props.surface_ten*dat.pore_rad))+(1e-20);
	dat.drdt = blake*(((3 * props.gas_const*dat.Temp*dat.dmdt) / (4 * PI)) +
		((3 * props.gas_const*dat.pore_mass*(-1 * dat.Coolrate)) / (4 * PI)) -
		((pow(dat.pore_rad, 3)*dat.dpdt)));
	
}

//Pre nucleation calculations
void pre_nuc(data_str& dat, simulation& sim, properties& props)
{ 
		double gs;
		
		dat.Temp = props.liq_temp - (dat.Coolrate*sim.time);
		gs = gscalc(dat.Temp);
		dat.Koz = (dat.Sec_dend / 180)*((pow(1 - gs, 3)) / (pow(gs, 2)));
	
		dat.pore_frac = dat.num_den*(4 / 3)*PI*(pow(dat.pore_rad, 3));
		dat.Coeff_diff_D1 = props.const_dif_Do*exp((-props.const_dif_a2) / dat.Temp);
		dat.Ke = (pow(10, -((props.const_equil_a1 / dat.Temp) + props.const_equil_b1)))*.01;
		dat.Ap = 4 * PI*pow(dat.pore_rad, 2)*pow((1 - gs), sim.imp_fac);
		dat.Cl = (props.Co - (dat.pore_frac*(dat.pore_density / props.liq_density))) /
			(1 - dat.pore_frac + gs*((1 + props.solid_shrink)*props.h_part_co - 1));
		dat.Clp = dat.Cl;
		dat.pore_pres = dat.atm_pres*(pow(((dat.Clp*props.act_co) / dat.Ke), 2));
		dat.pore_density = dat.pore_pres / (props.gas_const*dat.Temp);
		masscalcpre(dat, sim, props);
		Liq_prescalcpre(dat, sim, props, gs);
		radcalcpre(dat, sim, props);
		sim.Nuke_check = (dat.atm_pres*pow(((dat.Clp*props.act_co) / dat.Ke), 2)) - dat.liq_pres;
			
	}









// The next set of equations is for after nucleation
// Calculating Pore Radius change post Nucleation
void radcalcpost(data_str& dat, simulation& sim, properties& props, double Pre_rad)
{
	double blake;
	blake = 1 / ((3 * dat.liq_pres*pow(dat.pore_rad, 2)) + (4 * props.surface_ten*dat.pore_rad));
	//Blake singularity fix
	//blake = 1 / ((3 * dat.liq_pres*pow(dat.pore_rad, 2)) + (4 * props.surface_ten*dat.pore_rad))+(1e-20);
	dat.drdt = blake*(((3 * props.gas_const*dat.Temp*dat.dmdt) / (4 * PI)) +
		((3 * props.gas_const*dat.pore_mass*(-1 * dat.Coolrate)) / (4 * PI)) -
		((pow(dat.pore_rad, 3)*dat.dpdt)));
	dat.pore_rad = Pre_rad + (dat.drdt*sim.dt);

}
 

// Calculating Pore Mass Change Post Nucleation
void masscalcpost(data_str& dat, simulation& sim, properties& props, double Pre_mass)
{
	dat.dmdt = props.liq_density*dat.Coeff_diff_D1*dat.Ap*((dat.Cl - dat.Clp) / dat.pore_rad);
	dat.pore_mass = dat.pore_mass + (dat.dmdt*sim.dt);
}



// Calculating Liquid Presure Change Post Nucleation
void Liq_prescalcpost(data_str& dat, simulation& sim, properties& props, double gs, double Pre_pres)
{
	dat.dpdt = -1 * ((props.liq_dyn_visc*dat.Coolrate) / (dat.Koz*pow(dat.Ny, 2)))*
		(((props.solid_shrink)*(1 - dat.gpl - gs)) - (dat.gpl - dat.pore_frac));
	dat.liq_pres = dat.liq_pres + dat.dpdt*sim.dt;

}


// Post Nucleation Conditions
void post_nuc(data_str& dat, simulation& sim, properties& props)
{
	double gs;
	double Pre_mass = dat.pore_mass;
	double Pre_pres = dat.liq_pres;
	double Pre_rad = dat.pore_rad;
	
	dat.Temp = props.liq_temp - (dat.Coolrate*sim.time);
	gs = gscalc(dat.Temp);
	dat.Koz = (dat.Sec_dend / 180)*((pow(1 - gs, 3)) / (pow(gs, 2)));

	int iflag = 1;

	while (iflag)
	{
		sim.i = sim.i + 1;
		double rad_check = dat.pore_rad;

		dat.pore_pres = dat.liq_pres + ((2 * props.surface_ten) / dat.pore_rad);
		dat.pore_density = dat.pore_pres / (props.gas_const*dat.Temp);
		dat.pore_frac = dat.num_den*(4 / 3)*PI*(pow(dat.pore_rad, 3));
		dat.Coeff_diff_D1 = props.const_dif_Do*exp((-props.const_dif_a2) / dat.Temp);
		dat.Ke = (pow(10, -((props.const_equil_a1 / dat.Temp) + props.const_equil_b1)))*.01;
		dat.Ap = 4 * PI*pow(dat.pore_rad, 2)*pow((1 - gs), sim.imp_fac);
		dat.Clp = (dat.Ke / props.act_co)*sqrt(dat.pore_pres / dat.atm_pres);
		dat.Cl = (props.Co - dat.pore_frac*(dat.pore_density / props.liq_density)) /
			(1 - dat.pore_frac + gs*((1 + props.solid_shrink)*props.h_part_co - 1));
		masscalcpost(dat, sim, props, Pre_mass);
		Liq_prescalcpost(dat, sim, props, gs, Pre_pres);
		radcalcpost(dat, sim, props, Pre_rad);
		double rad_Tol = (dat.pore_rad - rad_check) / rad_check;
		if (abs(rad_Tol) < sim.Tol)
			iflag = 0;

	}

}