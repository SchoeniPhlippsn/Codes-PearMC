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


#include "Headers/Headers.h"

int main(int argc, char** argv){

	Nc=2; //# of spherocylinder
	Ns=1298; //# of spheres
	seed=42;
	rho0=0.6; //density
	pos_lambda = 0.02;
	ori_lambda = 0.02;
	vproc = 0.99;
    	l.resize(3);
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
	distN = 0.5*a1*(aspect-1.0/aspect);
	rlistP = 2*(b1-distN)+0.01;
	
    	RHO = toString(rho0);
    	NC = toString(Nc);
    	NS = toString(Ns);
    	VR = toString(Vratio);
    	KTH = toString(kth);
    	ASP = toString(aspect);
	
    	//Vsys = (2.65072*a1*a1*b1+0.687223*a1*a2*b1+0.147262*a2*a2*b1+0.687223*a1*a3*b1+0.147262*a3*a3*b1);       

    	std::cout << "V_c="<< Vsys << std::endl;

    	rsphere = pow(3.0/4.0*Vratio*Vsys/M_PI,1.0/3.0);

    	std::cout << rsphere << std::endl;

    	Vsys = (Nc + Ns*Vratio)*Vsys;

    	N = Nc + Ns;

    	Vsys /= N;

	rlistS = rsphere + 0.01;
	if( Ns != 0 ) maxpos = 2*rlistS;
	else maxpos = rlistP;

	gen.seed(seed);

    	step = 0;
    	std::cout << "Starting Initialisation of the system!" << std::endl; 
	Init();
    	std::cout << std::endl;
    	std::cout << "Particles are placed successfully!" << std::endl;

    	Compression();

    	step++;

    	acc = 0;   
	file_pre = "Results/Config_Nc" + NC + "_Ns"+ NS + "_Vr" + VR + "_rho" + RHO + "_step";
    	while(step <= finstep){

        	Move_step();
	        if(step % 100 == 0){          
			if(step % savestep == 0){          
				std::string file_name = file_pre + toString(step) + ".dat";
				write(file_name,0);
			}
			if(step % logstep == 0) writeLog(logfile);
			if(step % backupstep == 0){ 
				
				std::string newname = "Save/ConfigB.dat";

				rename(savefile.c_str(),newname.c_str());

				write(savefile,1);
			}

           		acceptance = static_cast<double>(acc)/(100*N);
			std::cout <<  step << " " << acceptance << std::endl;
			if (acceptance > 0.55){
			    pos_lambda += 0.01;
			    if( pos_lambda > maxpos) pos_lambda = maxpos;
			    ori_lambda += 0.01;
			    if( ori_lambda > maxpos) ori_lambda = maxpos;
			}
			if (acceptance < 0.45){
			    pos_lambda -= 0.001;
			    if( pos_lambda < 0.005 ) pos_lambda = 0.005;
			    ori_lambda -= 0.001;
			    if( ori_lambda < 0.005 ) ori_lambda = 0.005;
			}
            		acc=0;
        	}
		if(step % compstep == 0 && WallMove) Compression_step();
        	step++;
	}
    	write("Results/finaldat",1);

}
