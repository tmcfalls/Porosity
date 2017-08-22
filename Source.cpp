#include <cmath>
#include <iostream>
#include <fstream>

#include "Header.h"

using namespace std; 

//Function declarations
void init_props(properties&);
void init_dat(data_str&);
void init_sim(simulation&, properties&, data_str&);
void pre_nuc(data_str&, simulation&, properties&);
void post_nuc(data_str&, simulation&, properties&);




int main()
{
	ofstream file("output.csv");
	file << "T" << "," << "Liq Pres" << "," << "Pore Rad" <<
		"," << "Pore Pre" << endl;

	//Declare structs
	data_str dat;
	properties props;
	simulation sim; 

	int n = 0;
	sim.Nuke_flag = 1;
	sim.gpl_flag = 1;


	//while (sim.gpl_flag)
	//{
	//Initialize some values
	
		//Matl properties
	init_props(props);

		//Initial values for "data" struct
	init_dat(dat);

		//Initialize simulation parameters
	init_sim(sim, props, dat);

		
	
	
		while (sim.time_flag)
		{
			n = n + 1;
			sim.i = 0;
			sim.time = sim.time_start + (n*sim.dt);
			if (sim.Nuke_flag == 1)
			{
				pre_nuc(dat, sim, props);
				if ((sim.Nuke_check >= sim.Nuke_val) && (dat.drdt >= 0))
					sim.Nuke_flag = 0;
			}
			else
				post_nuc(dat, sim, props);
			
			if (n / sim.Nmax >= 0.9)
				sim.time_flag = 0;
			
			
			if (n % 100000 == 0)
			{
				cout << n/sim.Nmax << " " << dat.liq_pres
				<< " " << dat.pore_rad << " " << dat.pore_pres<<
				" " << sim.Nuke_flag << " " << sim.i << endl;
				
				file << n / sim.Nmax << "," << dat.liq_pres
				<< "," << dat.pore_rad << "," << dat.pore_pres
				<< "," << sim.Nuke_flag << endl;
			
			}
}

	
	file.close();
	cin.ignore();
	return 0;
}