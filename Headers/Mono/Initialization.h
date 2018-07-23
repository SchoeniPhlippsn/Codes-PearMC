void Init(){

    	RHO = toString(rho0);
    	NC = toString(Nc);
    	NS = toString(Ns);
    	VR = toString(Vratio);
    	KTH = toString(kth);
    	ASP = toString(aspect);
	
    	//Vsys = (2.65072*a1*a1*b1+0.687223*a1*a2*b1+0.147262*a2*a2*b1+0.687223*a1*a3*b1+0.147262*a3*a3*b1);       

    	std::cout << "V_c="<< Vsys << std::endl;

    	Vsys = Nc*Vsys;

    	N = Nc + Ns;

    	Vsys /= N;

	gen.seed(seed);

    	step = 0;

	for (int i=0; i<2; i++){
		file_pre = "Meshes/Pear" + ASP + "-" + KTH + ".inp";
		std::ifstream iFile(file_pre.c_str());

		if(!iFile){
			std::cerr << "The file " << file_pre << " is missing!" << std::endl;
			exit(-1);
		}
		iFile >> num_tri;

		std::cout<<"Reading object "<<i<<"\n";
		pear_mesh.NewObject(&(id[i]));
		std::cout<<"Adding triangles\n";
		      
		for (int j=1; j<=num_tri; j++){
			double v1[3], v2[3], v3[3];
			iFile >> v1[0] >> v1[1] >> v1[2];
			iFile >> v2[0] >> v2[1] >> v2[2];
			iFile >> v3[0] >> v3[1] >> v3[2];
			
			v1[0] *= a1*2;	
			v1[1] *= a1*2;	
			v1[2] *= a1*2;	
			v2[0] *= a1*2;	
			v2[1] *= a1*2;	
			v2[2] *= a1*2;	
			v3[0] *= a1*2;	
			v3[1] *= a1*2;	
			v3[2] *= a1*2;	

			pear_mesh.AddTri(v1, v2, v3, j);
		}
		      
		iFile.close();
	      
		std::cout<<"Calling finish_object with " << num_tri << " triangles\n";
		pear_mesh.EndObject();      
	     
		std::cout<<"Inserted object "<<i<<"\n";
	}

	file_pre = "Meshes/Bezier" + ASP + "-" + KTH + ".txt";
	std::ifstream iFile(file_pre.c_str());
	if(!iFile){
		std::cerr << "The file " << file_pre << " is missing!" << std::endl;
		exit(-1);
	}
	std::cout<<"Reading sphere-particle interface\n";

	max_x = 0;
	max_z = 0;
	while( iFile.eof() == false ){
		double Besx, Besz;
		iFile >> Besx >> Besz;
		if( iFile.eof()== false ){
			if(max_x < fabs(2*a1*Besx)) max_x = fabs(2*a1*Besx);
			if(max_z < fabs(2*a1*Besz)) max_z = fabs(2*a1*Besz);
		}
	}
	iFile.close();

    	rcut_P = 4*max_z*max_z;
    	rcut_P_T = (max_z+max_x)*(max_z+max_x);
    	rcut_P_II = 4*max_x*max_x;

	rlistP = 0.5*(sqrt(3*b1*b1-2*b1*max_x+3*max_x*max_x)-b1+max_x);
	distN = b1-rlistP;
	maxpos = 2*rlistP;

	std::cout << "rlist " << rlistP << std::endl;

	
	
	ln.resize(3);
    
	savefile = "Save/Config.dat";

	iFile.open(savefile.c_str());
	if(!iFile){
		 std::cerr << "There is no " << savefile << " file!" << std::endl;
        	exit(-1);
	}
	iFile.close();


	logfile = "Results/Log_Nc" + NC + "_Ns"+ NS + "_Vr" + VR + "_rho" + RHO + ".dat";

    	std::ifstream DatFile(logfile.c_str() );
	if(!DatFile){ 
    		std::ofstream DataFile(logfile.c_str() );
		DataFile << "step\tmsdP" << std::endl;
		DataFile.close();
	}
	DatFile.close();

    	l.resize(3);
    	l_2.resize(3);


        initCompress=false;

    	WP.resize(3,0); 
    	wP.resize(3,0); 
        read(savefile);

    	RenewList();
}
