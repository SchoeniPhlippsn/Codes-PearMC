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

	std::cout << "This version contains the neighbourhood list optimized for 2 pears in a pool of spheres!" << std::endl;

	Config.Nc=2; //# of spherocylinder
	Config.Ns=2398; //# of spheres
	seed=42;
	rho0=0.5; //density
	pos_lambda = 0.02;
	ori_lambda = 0.02;
	vproc = 0.99;
    	Config.l.resize(3);
    	Vratio = 0.04;
	finstep = 1000000;
	compstep = 100;
	savestep = 100;
	logstep = 100;
	backupstep = 100;
	WallMove=false;

	kth=3.8;
	aspect=3;
    	Config.Vsys=0.577525;

	for (int argi=1; argi<argc; argi++){
        std::string argvi = argv[argi];
        int next = argc-argi-1;
        if (argvi=="-Ns" && next>0) Config.Ns      = fromString<long>(argv[++argi]);	
        else if (argvi=="-Nc" && next>0) Config.Nc      = fromString<long>(argv[++argi]);	
		else if (argvi=="-s" && next>0) seed  = fromString<long>(argv[++argi]);
		else if (argvi=="-fs" && next>0) finstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-cs" && next>0) compstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-sa" && next>0) savestep  = fromString<long>(argv[++argi]);
		else if (argvi=="-log" && next>0) logstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-v" && next>0) vproc = fromString<double>(argv[++argi]);
		else if (argvi=="-Vr" && next>0) Vratio = fromString<double>(argv[++argi]);
		else if (argvi=="-Vp" && next>0) Config.Vsys = fromString<double>(argv[++argi]);
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
	distN = 0.5*b1;
	rlist = 2*sqrt(distN*distN+a1*a1)+0.02;
	if(rlist < 2*(b1-distN)) rlist = 2*(b1-distN)+0.02;
	
    	RHO = toString(rho0);
    	NC = toString(Config.Nc);
    	NS = toString(Config.Ns);
    	VR = toString(Vratio);
    	KTH = toString(kth);
    	ASP = toString(aspect);
	
    	//Config.Vsys = (2.65072*a1*a1*b1+0.687223*a1*a2*b1+0.147262*a2*a2*b1+0.687223*a1*a3*b1+0.147262*a3*a3*b1);       

    	std::cout << "V_c="<< Config.Vsys << std::endl;

    	rsphere = pow(3.0/4.0*Vratio*Config.Vsys/M_PI,1.0/3.0);

    	std::cout << rsphere << std::endl;

    	Config.Vsys = (Config.Nc + Config.Ns*Vratio)*Config.Vsys;

    	N = Config.Nc + Config.Ns;

    	Config.Vsys /= N;

	if( Config.Ns != 0 ) rlist = 2*rsphere + 0.02;
    	rlist_2 = rlist*rlist;  
    	maxpos = 0.5*rlist;


	gen.seed(seed);

    	Config.step = 0;
    	std::cout << "Starting Initialisation of the system!" << std::endl; 
	Init();
    	std::cout << std::endl;
    	std::cout << "Particles are placed successfully!" << std::endl;

    	Compression();

    	Config.step++;

    	acc = 0;   
	file_pre = "Results/Config_Nc" + NC + "_Ns"+ NS + "_Vr" + VR + "_rho" + RHO + "_step";
    	while(Config.step <= finstep){

        	Move_step();
	        if(Config.step % 100 == 0){          
			if(Config.step % savestep == 0){          
				std::string file_name = file_pre + toString(Config.step) + ".dat";
				Config.write(file_name,0);
			}
			if(Config.step % logstep == 0) Config.writeLog(logfile);
			if(Config.step % backupstep == 0){ 
				
				char oldname[] = "Save/Config.dat";
				char newname[] = "Save/ConfigB.dat";
				rename(oldname,newname);

				std::ofstream File ("Save/status.txt");

				File << 2 << std::endl;

				File.close();
				Config.write("Save/Config.dat",1);

				File.open("Save/status.txt");

				File << 1 << std::endl;

				File.close();
			}

           		acceptance = static_cast<double>(acc)/(100*N);
			std::cout << Config.step << " " << acceptance << std::endl;
            		if (acceptance > 0.55){
                		if( pos_lambda < maxpos) pos_lambda += 0.05;
                		else pos_lambda = maxpos;
                		if( ori_lambda < maxpos) ori_lambda += 0.05;
                		else ori_lambda = maxpos;
            		}
            		if (acceptance < 0.45){
                		if( pos_lambda < 0.01 ) pos_lambda = 0.01;
                		else pos_lambda -= 0.01;
                		if( ori_lambda < 0.01 ) ori_lambda = 0.01;
                		else ori_lambda -= 0.01;
            		}
            		acc=0;
        	}
		if(Config.step % compstep == 0 && WallMove) Compression_step();
        	Config.step++;
	}
    	Config.write("Results/finalConfig.dat",1);

    	std::ofstream File ("Save/status.txt");

    	File << 0 << std::endl;

    	File.close();
}
