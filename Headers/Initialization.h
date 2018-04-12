void Init(){

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
			bezier_x.push_back(2*a1*Besx);
			bezier_z.push_back(2*a1*Besz);


			if(max_x < fabs(2*a1*Besx)) max_x = fabs(2*a1*Besx);
			if(max_z < fabs(2*a1*Besz)) max_z = fabs(2*a1*Besz);
		}
	}
	iFile.close();

    	rcut_SPH = 4*rsphere*rsphere;
    	rcut_PSPH = (max_z+rsphere)*(max_z+rsphere);
    	rcut_PSPH_T = (rsphere + max_x)*(rsphere + max_x);
    	rcut_P = 4*max_z*max_z;
    	rcut_PII = (max_z+0.5*rlist)*(max_z+0.5*rlist);
    	rcut_PT = (0.5*rlist + max_x)*(0.5*rlist + max_x);
	
	ln.resize(3);
    
	std::cout << rcut_PII << " " << rcut_PT << std::endl;

	std::cout<<"Interface inserted with " << bezier_x.size() << " vertices\n";

    	std::ifstream status("Save/status.txt");

    	if(status){ 
		status >> restart;
		iFile.open("Save/Config.dat");
		if(!iFile){
			 std::cerr << "There is no Save/Config.dat file!" << std::endl;
        		exit(-1);
		}
    	}else{
        	std::cerr << "There is no Save/status.txt file! Please indicate if I have to start from the beginning" << std::endl;
        	exit(-1);
    	}

	logfile = "Results/Log_Nc" + NC + "_Ns"+ NS + "_Vr" + VR + "_rho" + RHO + ".dat";

    	std::ifstream DatFile(logfile.c_str() );
	if(!DatFile){ 
    		std::ofstream DataFile(logfile.c_str() );
		DataFile << "step\tu\tRR\trr\tRr\tmsd" << std::endl;
		DataFile.close();
	}
	DatFile.close();

    	Config.l.resize(3);
    	Config.l_2.resize(3);
    	if(restart!=0){ 
        	if(restart==1) Config.read("Save/Config.dat");
        	else Config.read("Save/ConfigB.dat");
    	}else{

        	Config.rhoV = 0.01;
        	Config.rhoN = Config.rhoV/Config.Vsys;
        	Config.Vbox = N/Config.rhoN;
            
        	Config.l[0] =pow(Config.Vbox,1.0/3.0);
        	Config.l[1] = Config.l[0];
        	Config.l[2] = Config.l[0];

		Config.l_2[0] = Config.l[0]*0.5;
		Config.l_2[1] = Config.l[1]*0.5;
		Config.l_2[2] = Config.l[2]*0.5;

        
        	std::cout << "System size\tlx=" << Config.l[0] << "\tly=" << Config.l[1] << "\tlz=" << Config.l[2] << std::endl;



        	int grid = pow(N,1.0/3.0);
        	grid++; 
        
        	int i_min = 1;
        	int i = 1;
        	int j = 1;
        	int k = 1;
        	int bk = 1;
        	double stepi = Config.l[0]/grid;
        	double stepj = Config.l[1]/grid;
        	double stepk = Config.l[2]/grid;
		
        	for (int ii=0; ii < N; ii++){			
            		bool inside =false;
            		double rmin = 1000;
            		int h=0;
            		class teilchen newpart;
			newpart.dist = distN;
			newpart.already = false;
            		newpart.pos.push_back( i*stepi );
            		newpart.pos.push_back( j*stepj );
            		newpart.pos.push_back( k*stepk );

            		newpart.pos_msd.push_back( 0 );
            		newpart.pos_msd.push_back( 0 );
            		newpart.pos_msd.push_back( 0 );

            		if( k % 2 == 0 ){ 
                		newpart.pos[0] += 0.5*stepi; 
                		newpart.pos[1] += 0.5*stepj; 
            		}

			if(newpart.pos[0] > Config.l[0]) newpart.pos[0] -= Config.l[0];
			else if(newpart.pos[0] < 0 ) newpart.pos[0] += Config.l[0];

			if(newpart.pos[1] > Config.l[1]) newpart.pos[1] -= Config.l[1];
			else if(newpart.pos[1] < 0 ) newpart.pos[1] += Config.l[1];

			if(newpart.pos[2] > Config.l[2]) newpart.pos[2] -= Config.l[2];
			else if(newpart.pos[2] < 0 ) newpart.pos[2] += Config.l[2];

            		if(ii < Config.Nc ){
                		newpart.ori.push_back(0); 
                		newpart.ori.push_back(0); 
                		newpart.ori.push_back(bk); 
            
                		for( int v=0; v< Config.part.size(); v++){
					if(v < Config.Nc ){
						if(overlapP(newpart,Config.part[v],Config.l)){
							inside = true;
							break;   
						}
					}else{
						if(overlapPSPH(newpart,Config.part[v],Config.l)){
							inside = true;
							break;   
						}
					}
				}	

				newpart.trans[0][0]=1;
				newpart.trans[0][1]=0;
				newpart.trans[0][2]=0;
				newpart.trans[0][3]=0;

				newpart.trans[1][0]=0;
				newpart.trans[1][1]=newpart.ori[2];
				newpart.trans[1][2]=0;
				newpart.trans[1][3]=0;

				newpart.trans[2][0]=0;
				newpart.trans[2][1]=0;
				newpart.trans[2][2]=newpart.ori[2];
				newpart.trans[2][3]=0;

				newpart.trans[3][0]=0;
				newpart.trans[3][1]=0;
				newpart.trans[3][2]=0;
				newpart.trans[3][3]=1;
            		}else{
                		for( int v=0; v< Config.part.size(); v++){
                    			if(v < Config.Nc ){
                        			if(overlapPSPH(Config.part[v],newpart,Config.l)){
                            				inside = true;
                            				break;   
                        			}
                    			}else{
                        			if(overlapSPH(newpart,Config.part[v],Config.l) ){
                            				inside = true;
                            				break;   
                        			}
                    			}
                		}	
            		}
            		if(inside){
               	 		std::cerr << "Failed to place particle! Intersection! " << grid << " " << i << "  " << j << " " << k << std::endl;
                		exit(1);
            		}


            		Config.part.push_back(newpart);

            		i+=2;
            		if (i > grid && ii != N-1){
                		j++;
                		i=i_min;
                		if (j > grid){
                    			k++;
                    			bk*=-1;
                    			j=1;
                    			if(k>grid){
                        			if( i_min == 2 ){
                            				std::cerr << "There is a problem with the placement of the particles!" << std::endl;
                            				exit(1);
                        			}
                        			std::cout << ii << std::endl;
                        			j = 1;
                        			k = 1;
                        			i_min =2;
                        			i = 2;
                    			}
                		}
            		}
        	}
    	}
        initCompress=false;

    	Config.W.resize(3,0); 
    	Config.w.resize(3,0); 

    	RenewList();
}
