#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm> 
#include <stdio.h>
#include "VCollide.H"
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
	
boost::random::mt19937 gen;	// random number generator [0,1]
boost::random::uniform_01<> uni;

#ifdef POLY
#include "Headers/Poly/Headers.h"
#else
#include "Headers/Mono/Headers.h"
#endif


int main(int argc, char** argv){

	Nc=380; //# of spherocylinder
	Ns=0; //# of spheres
	seed=42;
	rho0=0.6; //density
	pos_lambda = 0.02;
	ori_lambda = 0.02;
	vproc = 0.99;
    	Vratio = 0.005;
	finstep = 100000000;
	compstep = 100;
	savestep = 50000;
	logstep = 100;
	backupstep = 100;
	WallMove=false;

	kth=3.8;
	aspect=3;
    	Vsys=0.577525;

	for (int argi=1; argi<argc; argi++){
        std::string argvi = argv[argi];
        int next = argc-argi-1;
        if (argvi=="-Ns" && next>0) Ns      = fromString<long>(argv[++argi]);	
        else if (argvi=="-Nc" && next>0) Nc      = fromString<long>(argv[++argi]);	
		else if (argvi=="-s" && next>0) seed  = fromString<long>(argv[++argi]);
		else if (argvi=="-fs" && next>0) finstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-cs" && next>0) compstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-sa" && next>0) savestep  = fromString<long>(argv[++argi]);
		else if (argvi=="-log" && next>0) logstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-v" && next>0) vproc = fromString<double>(argv[++argi]);
		else if (argvi=="-Vr" && next>0) Vratio = fromString<double>(argv[++argi]);
		else if (argvi=="-Vp" && next>0) Vsys = fromString<double>(argv[++argi]);
		else if (argvi=="-r" && next>0) rho0 = fromString<double>(argv[++argi]);
		else if (argvi=="-pl" && next>0) pos_lambda = fromString<double>(argv[++argi]);
		else if (argvi=="-ol" && next>0) ori_lambda = fromString<double>(argv[++argi]);
		else if (argvi=="-kth" && next>0) kth = fromString<double>(argv[++argi]);
		else if (argvi=="-asp" && next>0) aspect = fromString<double>(argv[++argi]);
		else if (argvi=="-wm" && next>0) WallMove = fromString<bool>(argv[++argi]);
		else throw std::runtime_error("no such argument: "+argvi);
  	}

	a2=a1*(kth-2.0/3.0*aspect)/kth;
	a3=a1*(kth+2.0/3.0*aspect)/kth;
	b1 = aspect*a1; 
    	std::cout << "Starting Initialisation of the system!" << std::endl; 
	Init();
    	std::cout << std::endl;
    	std::cout << "Particles are placed successfully!" << std::endl;

    	Compression();

    	step++;

    	acc = 0;   
	file_pre = "Results/rho" + RHO + "/Config_Nc" + NC + "_Ns"+ NS + "_Vr" + VR + "_rho" + RHO + "_step";
	std::string newname = "Save/ConfigB" + RHO + ".dat";

	sstep = step % savestep;
	lstep = step % logstep;
	bstep = step % backupstep;

	double logstepN = 1.0/(logstep*N);
    	while(step <= finstep){

        	Move_step();
	        //if(step % 100 == 0){          
			//if(step % savestep == 0){          
			if(sstep == savestep ){          
				file_name = file_pre + toString(step) + ".dat";
				writeConfig();
				sstep=0;
			}
			if(lstep == logstep){ 
				writeLog();
				lstep =0;

				acceptance = acc*logstepN;
				std::cout <<  step << " " << acceptance << std::endl;
				if (acceptance > 0.55){
				    pos_lambda += 0.01;
				    if( pos_lambda > maxpos) pos_lambda = maxpos;
				    ori_lambda += 0.01;
				    if( ori_lambda > maxpos) ori_lambda = maxpos;
				}else{
					if (acceptance < 0.45){
					    pos_lambda -= 0.001;
					    if( pos_lambda < 0.005 ) pos_lambda = 0.005;
					    ori_lambda -= 0.001;
					    if( ori_lambda < 0.005 ) ori_lambda = 0.005;
					}
				}
				acc=0;
			}
			if(bstep == backupstep){ 
				
				rename(savefile.c_str(),newname.c_str());
				writeSave();
				bstep=0;
			}

        	//}
//		if(step % compstep == 0 && WallMove) Compression_step();
        	step++;
        	lstep++;
        	bstep++;
        	sstep++;
	}
	writeSave();

}
