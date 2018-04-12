void system::write (std::string file, bool save){

    std::ofstream oFile (file.c_str() );
    if(!save){
        oFile << Nc + Ns << std::endl;
        oFile << -l_2[0] << " " << l_2[0] << std::endl;
        oFile << -l_2[1] << " " << l_2[1] << std::endl;
        oFile << -l_2[2] << " " << l_2[2] << std::endl;
    }else{
        oFile << Nc << " " << Ns << std::endl;
        oFile << l[0] << " " << l[1] << " " << l[2] << std::endl;
        oFile << step << std::endl;

        std::ofstream ooFile("Save/status.txt");
        ooFile << 1 << std::endl;
        ooFile.close();
    }

    for( v = 0 ; v < Nc; v++){

        if(!save){ 

            oFile << "PEAR(" << aspect << "," << kth << ") ";
            oFile <<  part[v].pos[0] - int(part[v].pos[0]/l_2[0])*l[0] << " " << part[v].pos[1] - int(part[v].pos[1]/l_2[1])*l[1] << " " << part[v].pos[2] - int(part[v].pos[2]/l_2[2])*l[2]<< " ";
            oFile << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " " << 1 << std::endl; 
        }else{
            oFile <<  part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " "; 
            oFile <<  part[v].pos_msd[0] << " " << part[v].pos_msd[1] << " " << part[v].pos_msd[2] << " "; 
            oFile <<  part[v].ori[0] << " " << part[v].ori[1] << " " << part[v].ori[2] << std::endl; 
        }
    }
    for( v = Nc ; v < part.size(); v++){

        if(!save){ 
            oFile << "SPHERE ";
            oFile <<  part[v].pos[0] - int(part[v].pos[0]/l_2[0])*l[0] << " " << part[v].pos[1] - int(part[v].pos[1]/l_2[1])*l[1] << " " << part[v].pos[2] - int(part[v].pos[2]/l_2[2])*l[2]<< " ";
            oFile << 2*rsphere << " 0 0 0 " << 2*rsphere << " 0 0 0 " << 2*rsphere << " " << v << std::endl;
        }else{
            oFile << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " "; 
            oFile <<  part[v].pos_msd[0] << " " << part[v].pos_msd[1] << " " << part[v].pos_msd[2] << " "; 
            oFile << rsphere <<  std::endl; 
        }
    }
    oFile.close();
}

void system::read (std::string file){

    std::ifstream iFile (file.c_str() );

    if(!iFile){
        std::cerr << "Can not read " << file << "!" << std::endl;
        exit(-1);
    }
    iFile >> Nc >> Ns;
    part.resize(Nc+Ns);
    l.resize(3);
    
    iFile >> l[0] >> l[1] >> l[2];

    l_2[0] = l[0]*0.5;
    l_2[1] = l[1]*0.5;
    l_2[2] = l[2]*0.5;

    iFile >> step;

    Vbox=l[0]*l[1]*l[2];

    rhoN = part.size()/Vbox;
    rhoV = Vsys*rhoN;

    for( v = 0 ; v < Nc; v++){
        part[v].pos.resize(3);
        part[v].pos_msd.resize(3);
        part[v].ori.resize(3);
	part[v].dist = distN;
	part[v].already = false;

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2] >> part[v].pos_msd[0] >> part[v].pos_msd[1] >> part[v].pos_msd[2] >> part[v].ori[0] >> part[v].ori[1] >>  part[v].ori[2]; 

	if(part[v].pos[0] < 0 ) part[v].pos[0] += l[0];
	if(part[v].pos[1] < 0 ) part[v].pos[1] += l[1];
	if(part[v].pos[2] < 0 ) part[v].pos[2] += l[2];

	if(step < 100) part[v].pos_msd[0] =0;
	if(step < 100) part[v].pos_msd[1] =0;
	if(step < 100) part[v].pos_msd[2] =0;

	sqrtz = 1-part[v].ori[2]*part[v].ori[2];
	if(sqrtz > 1e-4){
		sqrtz = 1/sqrt(sqrtz);
		part[v].trans[0][0] = part[v].ori[2]*part[v].ori[0]*sqrtz;
		part[v].trans[0][1] = -part[v].ori[1]*sqrtz;
		part[v].trans[0][2] = part[v].ori[0];

		part[v].trans[1][0] = part[v].ori[2]*part[v].ori[1]*sqrtz;
		part[v].trans[1][1] = part[v].ori[0]*sqrtz;
		part[v].trans[1][2] = part[v].ori[1];

		part[v].trans[2][0] = -1/sqrtz;
		part[v].trans[2][1] = 0;
		part[v].trans[2][2] = part[v].ori[2];
	}else{
		part[v].trans[0][0] = 1;
		part[v].trans[0][1] = 0;
		part[v].trans[0][2] = 0;
		
		part[v].trans[1][0] = 0;
		part[v].trans[1][2] = 0;


		part[v].trans[2][0] = 0;
		part[v].trans[2][1] = 0;

		if(part[v].ori[2]>0){
			part[v].trans[1][1] = 1;
			part[v].trans[2][2] = 1;
		}else{
			part[v].trans[1][1] = -1;
			part[v].trans[2][2] = -1;
		}
	}

	part[v].trans[0][3] = 0;
	part[v].trans[1][3] = 0;
	part[v].trans[2][3] = 0;

	part[v].trans[3][0] = 0;
	part[v].trans[3][1] = 0;
	part[v].trans[3][2] = 0;
	part[v].trans[3][3] = 1;

	for(vv = 0; vv < v; vv++ ){ 
		if (overlapP ( part[v], part[vv], l) ){
			std::cerr << "There is still a problem (" << vv << "," << v << ")!" << std::endl;

        		std::cerr << "PEAR " << part[vv].pos[0] << " " << part[vv].pos[1] << " " << part[vv].pos[2] << " " << part[vv].trans[0][0] << " " << part[vv].trans[1][0] << " " << part[vv].trans[2][0]  << " " << part[vv].trans[0][1] << " " << part[vv].trans[1][1] << " " << part[vv].trans[2][1] << " " <<  part[vv].trans[0][2] << " " << part[vv].trans[1][2] << " " << part[vv].trans[2][2] << " 1" << std::endl;

        		std::cerr << "PEAR " << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " " << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " 1" << std::endl;
			//exit(1);
		}
	}
    }

    for( v = Nc ; v < part.size(); v++){
        part[v].pos.resize(3);
        part[v].ori.resize(3);

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2]  >> part[v].pos_msd[0] >> part[v].pos_msd[1] >> part[v].pos_msd[2]>> rsphere; 

	if(part[v].pos[0] < 0 ) part[v].pos[0] += l[0];
	if(part[v].pos[1] < 0 ) part[v].pos[1] += l[1];
	if(part[v].pos[2] < 0 ) part[v].pos[2] += l[2];

	for(vv = 0; vv < Nc; vv++ ){ 
		if (overlapPSPH ( part[vv], part[v], l) ){
			std::cerr << "There is still a problem (" << vv << "," << v << ")!" << std::endl;
			exit(1);
		}
	}
	for(vv = Nc; vv < v; vv++ ){ 
		if (overlapSPH( part[v], part[vv], l) ){
			std::cerr << "There is still a problem (" << vv << "," << v << ")!" << std::endl;
    			std::vector<double> RTest = DeltaR(part[vv].pos,part[v].pos,l);
			double testing = scal_p(RTest,RTest);
			std::cerr << testing << " " << rcut_SPH << " " << rsphere << std::endl;
//			exit(1);
		}
	}
    }

    iFile.close();
}

void system::writeLog (std::string file){

    std::fstream fFile (file.c_str(), std::fstream::app );

    sqrtz = 0;
    a11 = 0;
    a22 = 0;
    a12 = 0;
    msd = 0;

    RR.resize(3);

    for( v = 0 ; v < Nc; v++){
    	for( vv = v+1 ; vv < Nc; vv++){
		sqrtz += part[v].ori[0]*part[vv].ori[0] + part[v].ori[1]*part[vv].ori[1] + part[v].ori[2]*part[vv].ori[2];

		RR[0] = part[v].pos[0]-part[vv].pos[0]-(part[v].ori[0]-part[vv].ori[0])*b1;
		if( RR[0] > l[0]/2 ) RR[0] -= l[0];
		else if( RR[0] < -l[0]/2 ) RR[0] += l[0];

		RR[1] = part[v].pos[1]-part[vv].pos[1]-(part[v].ori[1]-part[vv].ori[1])*b1;
		if( RR[1] > l[1]/2 ) RR[1] -= l[1];
		else if( RR[1] < -l[1]/2 ) RR[1] += l[1];

		RR[2] = part[v].pos[2]-part[vv].pos[2]-(part[v].ori[2]-part[vv].ori[2])*b1;
		if( RR[2] > l[2]/2 ) RR[2] -= l[2];
		else if( RR[2] < -l[2]/2 ) RR[2] += l[2];
		a11 += sqrt(RR[0]*RR[0] + RR[1]*RR[1] + RR[2]*RR[2]);


		RR[0] = part[v].pos[0]-part[vv].pos[0]+(part[v].ori[0]-part[vv].ori[0])*b1;
		if( RR[0] > l[0]/2 ) RR[0] -= l[0];
		else if( RR[0] < -l[0]/2 ) RR[0] += l[0];

		RR[1] = part[v].pos[1]-part[vv].pos[1]+(part[v].ori[1]-part[vv].ori[1])*b1;
		if( RR[1] > l[1]/2 ) RR[1] -= l[1];
		else if( RR[1] < -l[1]/2 ) RR[1] += l[1];

		RR[2] = part[v].pos[2]-part[vv].pos[2]+(part[v].ori[2]-part[vv].ori[2])*b1;
		if( RR[2] > l[2]/2 ) RR[2] -= l[2];
		else if( RR[2] < -l[2]/2 ) RR[2] += l[2];
		a22 += sqrt(RR[0]*RR[0] + RR[1]*RR[1] + RR[2]*RR[2]);


		RR[0] = part[v].pos[0]-part[vv].pos[0]+(part[v].ori[0]+part[vv].ori[0])*b1;
		if( RR[0] > l[0]/2 ) RR[0] -= l[0];
		else if( RR[0] < -l[0]/2 ) RR[0] += l[0];

		RR[1] = part[v].pos[1]-part[vv].pos[1]+(part[v].ori[1]+part[vv].ori[1])*b1;
		if( RR[1] > l[1]/2 ) RR[1] -= l[1];
		else if( RR[1] < -l[1]/2 ) RR[1] += l[1];

		RR[2] = part[v].pos[2]-part[vv].pos[2]+(part[v].ori[2]+part[vv].ori[2])*b1;
		if( RR[2] > l[2]/2 ) RR[2] -= l[2];
		else if( RR[2] < -l[2]/2 ) RR[2] += l[2];
		a12 += sqrt(RR[0]*RR[0] + RR[1]*RR[1] + RR[2]*RR[2]);


		RR[0] = part[v].pos[0]-part[vv].pos[0]-(part[v].ori[0]+part[vv].ori[0])*b1;
		if( RR[0] > l[0]/2 ) RR[0] -= l[0];
		else if( RR[0] < -l[0]/2 ) RR[0] += l[0];

		RR[1] = part[v].pos[1]-part[vv].pos[1]-(part[v].ori[1]+part[vv].ori[1])*b1;
		if( RR[1] > l[1]/2 ) RR[1] -= l[1];
		else if( RR[1] < -l[1]/2 ) RR[1] += l[1];

		RR[2] = part[v].pos[2]-part[vv].pos[2]-(part[v].ori[2]+part[vv].ori[2])*b1;
		if( RR[2] > l[2]/2 ) RR[2] -= l[2];
		else if( RR[2] < -l[2]/2 ) RR[2] += l[2];
		a12 += sqrt(RR[0]*RR[0] + RR[1]*RR[1] + RR[2]*RR[2]);

	}
	msd += part[v].pos_msd[0]*part[v].pos_msd[0]+part[v].pos_msd[1]*part[v].pos_msd[1]+part[v].pos_msd[2]*part[v].pos_msd[2];
	
    }
    msd = msd/Nc;
    
    sqrtz = 2*sqrtz/(Nc*Nc-Nc);
   
    a11 = 2*a11/(Nc*Nc-Nc);
    a22 = 2*a22/(Nc*Nc-Nc);
    a12 = a12/(Nc*Nc-Nc);

    fFile << step << "\t" << sqrtz << "\t" << a11 << "\t" << a22 << "\t" << a12 <<  "\t" << msd << std::endl;
	
    fFile.close();
}
