void write (std::string file, bool save){

    std::ofstream oFile (file.c_str() );
    if(!save){
        oFile << Nc + Ns << std::endl;
        oFile << -l_2[0] << " " << l_2[0] << std::endl;
        oFile << -l_2[1] << " " << l_2[1] << std::endl;
        oFile << -l_2[2] << " " << l_2[2] << std::endl;
    }else{
        oFile << Nc << " " << Ns << std::endl;
        oFile << step << " " << pos_lambda << " " << ori_lambda << std::endl;
        oFile << l[0] << " " << l[1] << " " << l[2] << std::endl;
    }

    for( int v = 0 ; v < Nc; v++){
        if(!save){ 

            oFile << "PEAR(" << aspect << "," << kth << ") ";
            oFile <<  part[v].pos[0] - int(part[v].pos[0]/l_2[0])*l[0] << " " << part[v].pos[1] - int(part[v].pos[1]/l_2[1])*l[1] << " " << part[v].pos[2] - int(part[v].pos[2]/l_2[2])*l[2] << " ";
            oFile << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " " << 1 << std::endl; 
        }else{
            oFile <<  part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " "; 
            oFile <<  part[v].pos_msd[0] << " " << part[v].pos_msd[1] << " " << part[v].pos_msd[2] << " "; 
            oFile <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << std::endl; 
        }
    }
    oFile.close();
}

void read (std::string file){

    std::ifstream iFile (file.c_str() );

    if(!iFile){
        std::cerr << "Can not read " << file << "!" << std::endl;
        exit(-1);
    }
    iFile >> Nc >> Ns;

    if(Ns != 0){
	std::cerr << "Save/Config.dat does not contain monodisperse data! Please compile file again :)" << std::endl;
	exit(0);
    }

    part.resize(Nc);
    l.resize(3);
    
    iFile >> step >> pos_lambda >> ori_lambda;


    if( step < 100 && step != 0){
	std::cerr << "Old Save file version! No information for the step sizes in second line [step move_step rot_step]" << std::endl;
	std::cerr << "I try to convert it :)" << std::endl;
	
	iFile >> l[2] >> step;
	l[0]= l[2];
	l[1]= l[2];
    	
	pos_lambda = -1;
	ori_lambda = 0.02;

	std:: cout << step << " " << pos_lambda << " " << ori_lambda << " " << l[2] << std::endl;
    }else iFile >> l[0] >> l[1] >> l[2];

    l_2[0] = l[0]*0.5;
    l_2[1] = l[1]*0.5;
    l_2[2] = l[2]*0.5;

    Vbox=l[0]*l[1]*l[2];

    rhoN = part.size()/Vbox;
    rhoV = Vsys*rhoN;

    WP[0] = (int)(l[0]/rlistP); 
    WP[1] = (int)(l[1]/rlistP); 
    WP[2] = (int)(l[2]/rlistP); 

    wP[0] = WP[0]/l[0];
    wP[1] = WP[1]/l[1];
    wP[2] = WP[2]/l[2];
	
    headP.resize(WP[0]*WP[1]*WP[2]); 
    usedCell.resize(WP[0]*WP[1]*WP[2],false); 
    for( int i=0; i<headP.size(); i++) headP[i] = -1;

    linkP.resize(Nc+part.size()); 
    for( int i=0; i<linkP.size(); i++) linkP[i] = -1;

	s_n.resize(WP[0]);

	s_nn[0] = 0;
	s_nn[1] = 1;
	s_nn[2] = WP[0]-1;
	s_nn[3] = 2;
	s_nn[4] = WP[0]-2;
	s_n[0]=s_nn;

	s_nn[0] = 1;
	s_nn[1] = 2;
	s_nn[2] = 0;
	s_nn[3] = 3;
	s_nn[4] = WP[0]-1;
	s_n[1]=s_nn;
	for( int i=2; i<s_n.size()-2; i++){ 
		s_nn[0] = i;
		s_nn[1] = i+1;
		s_nn[3] = i+2;
		s_nn[2] = i-1;
		s_nn[4] = i-2;

		s_n[i]=s_nn;
	}

	s_nn[0] = WP[0]-2;
	s_nn[1] = WP[0]-1;
	s_nn[2] = WP[0]-3;
	s_nn[3] = 0;
	s_nn[4] = WP[0]-4;

	s_n[WP[0]-2]=s_nn;

	s_nn[0] = WP[0]-1;
	s_nn[1] = 0;
	s_nn[2] = WP[0]-2;
	s_nn[3] = 1;
	s_nn[4] = WP[0]-3;

	s_n[WP[0]-1]=s_nn;

    for( int v = 0 ; v < Nc; v++){
        part[v].pos.resize(3);
        part[v].pos_msd.resize(3);
        part[v].dist_ori.resize(3);
	part[v].already = false;

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2] >> part[v].pos_msd[0] >> part[v].pos_msd[1] >> part[v].pos_msd[2] >> part[v].trans[0][2]  >> part[v].trans[1][2] >> part[v].trans[2][2]; 

//	std::cout << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " " << part[v].pos_msd[0] << " " << part[v].pos_msd[1] << " " << part[v].pos_msd[2] << " " << part[v].ori[0]  << " " << part[v].ori[1] << " " << part[v].ori[2] << std::endl;
//	exit(0);

	if(part[v].pos[0] < 0 ) part[v].pos[0] += l[0];
	if(part[v].pos[1] < 0 ) part[v].pos[1] += l[1];
	if(part[v].pos[2] < 0 ) part[v].pos[2] += l[2];

	if(step < 100) part[v].pos_msd[0] =0;
	if(step < 100) part[v].pos_msd[1] =0;
	if(step < 100) part[v].pos_msd[2] =0;

	sqrtz = 1-part[v].trans[2][2]*part[v].trans[2][2];
	if(sqrtz > 1e-5){
		sqrtz = 1/sqrt(sqrtz);
		part[v].trans[0][0] = part[v].trans[2][2]*part[v].trans[0][2]*sqrtz;
		part[v].trans[0][1] = -part[v].trans[1][2]*sqrtz;

		part[v].trans[1][0] = part[v].trans[2][2]*part[v].trans[1][2]*sqrtz;
		part[v].trans[1][1] = part[v].trans[0][2]*sqrtz;

		part[v].trans[2][0] = -1/sqrtz;
		part[v].trans[2][1] = 0;
	}else{
		part[v].trans[0][0] = 1;
		part[v].trans[0][1] = 0;
		part[v].trans[0][2] = 0;
		
		part[v].trans[1][0] = 0;
		part[v].trans[1][2] = 0;


		part[v].trans[2][0] = 0;
		part[v].trans[2][1] = 0;

		if(part[v].trans[2][2]>0){
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

	part[v].dist_ori[0]= distN*part[v].trans[0][2];
	part[v].dist_ori[1]= distN*part[v].trans[1][2];
	part[v].dist_ori[2]= distN*part[v].trans[2][2];

	NCell.resize(0);
	NPart.resize(0);

	dsx = part[v].pos[0]-part[v].dist_ori[0];
	if(dsx < 0) part[v].s[0] = WP[0]-1;
	else{
		if(dsx > l[0]) part[v].s[0] = 0;
		else part[v].s[0] = dsx*wP[0];
	}

	dsy = part[v].pos[1]-part[v].dist_ori[1];
	if(dsy < 0) part[v].s[1] = WP[1]-1;
	else{
		if(dsy > l[1]) part[v].s[1] = 0;
		else part[v].s[1] = dsy*wP[1];
	}

	dsz = part[v].pos[2]-part[v].dist_ori[2];
	if(dsz < 0) part[v].s[2] = WP[2]-1;
	else{
		if(dsz > l[2]) part[v].s[2] = 0;
		else part[v].s[2] = dsz*wP[2];
	}

	for (int iz=0; iz < 5 && !inside; iz++){
		kz = WP[1]*s_n[part[v].s[2]][iz];
		for (int iy=0; iy < 5 && !inside; iy++){
			ky = WP[0]*(s_n[part[v].s[1]][iy]+kz);
			for (int ix=0; ix < 5 && !inside; ix++){
				k = s_n[part[v].s[0]][ix] + ky;
				if( usedCell[k] ) continue;

				vv = headP[k];
			    
				while(vv != -1){
					newvv = vv;
					if(newvv >= part.size() ) newvv -= part.size(); 
					if(!part[newvv].already){ 
						NPart.push_back(newvv);
						part[newvv].already=true;	
						if (overlapP ( part[v], part[newvv], l) ){
							std::cerr << "There is still a problem (" << newvv << "," << v << ")!" << std::endl;

							std::cerr << "PEAR " << part[newvv].pos[0] << " " << part[newvv].pos[1] << " " << part[newvv].pos[2] << " " << part[newvv].trans[0][0] << " " << part[newvv].trans[1][0] << " " << part[newvv].trans[2][0]  << " " << part[newvv].trans[0][1] << " " << part[newvv].trans[1][1] << " " << part[newvv].trans[2][1] << " " <<  part[newvv].trans[0][2] << " " << part[newvv].trans[1][2] << " " << part[newvv].trans[2][2] << " 1" << std::endl;
							std::cerr << "PEAR " << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " " << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " 1" << std::endl;
						}
					}
					vv = linkP[vv];
				}
				NCell.push_back(k);
				usedCell[k]=true;
			}
		}
	}

	dsxN = part[v].pos[0]+part[v].dist_ori[0];
	if(dsxN < 0) part[v].sN[0] = WP[0]-1;
	else{
		if(dsxN > l[0]) part[v].sN[0] = 0;
		else part[v].sN[0] = dsxN*wP[0];
	}

	dsyN = part[v].pos[1]+part[v].dist_ori[1];
	if(dsyN < 0) part[v].sN[1] = WP[1]-1;
	else{
		if(dsyN > l[1]) part[v].sN[1] = 0;
		else part[v].sN[1] = dsyN*wP[1];
	}

	dszN = part[v].pos[2]+part[v].dist_ori[2];
	if(dszN < 0) part[v].sN[2] = WP[2]-1;
	else{
		if(dszN > l[2]) part[v].sN[2] = 0;
		else part[v].sN[2] = dszN*wP[2];
	}

	if(part[v].s[0] != part[v].sN[0] || part[v].s[1] != part[v].sN[1] || part[v].s[2] != part[v].sN[2] ){
		for (int iz=0; iz < 5 && !inside; iz++){
			kz = WP[1]*s_n[MovedParticle.sN[2]][iz];
			for (int iy=0; iy < 5 && !inside; iy++){
				ky = WP[0]*(s_n[MovedParticle.sN[1]][iy]+kz);
				for (int ix=0; ix < 5 && !inside; ix++){

					k = s_n[MovedParticle.sN[0]][ix] + ky;
					if( usedCell[k] ) continue;

					vv = headP[k];
			    
					while(vv != -1){
						newvv = vv;
						if(newvv >= part.size() ) newvv -= part.size(); 
						if(!part[newvv].already){ 
							NPart.push_back(newvv);
							part[newvv].already=true;	
							if (overlapP ( part[v], part[newvv], l) ){
								std::cerr << "There is still a problem (" << newvv << "," << v << ")!" << std::endl;

								std::cerr << "PEAR " << part[newvv].pos[0] << " " << part[newvv].pos[1] << " " << part[newvv].pos[2] << " " << part[newvv].trans[0][0] << " " << part[newvv].trans[1][0] << " " << part[newvv].trans[2][0]  << " " << part[newvv].trans[0][1] << " " << part[newvv].trans[1][1] << " " << part[newvv].trans[2][1] << " " <<  part[newvv].trans[0][2] << " " << part[newvv].trans[1][2] << " " << part[newvv].trans[2][2] << " 1" << std::endl;

								std::cerr << "PEAR " << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " " << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " 1" << std::endl;
							}
						}
						vv = linkP[vv];
					}
				}
			}
		}
	}

	for( vv = 0; vv < NPart.size(); vv++) part[NPart[vv]].already=false;
	for( vv = 0; vv < NCell.size(); vv++) usedCell[NCell[vv]]=false;

	int k = part[v].s[0] + WP[0]*(part[v].s[1]+WP[1]*part[v].s[2]); 

	part[v].cell = k;
	if(headP[k]==-1) headP[k] = v;
	else{
	   linkP[v] = headP[k];
	   headP[k] = v; 
	}

	k = part[v].sN[0] + WP[0]*(part[v].sN[1]+WP[1]*part[v].sN[2]); 

	part[v].cellN = k;
	v += part.size();
	if(headP[k]==-1) headP[k] = v;
	else{
	   linkP[v] = headP[k];
	   headP[k] = v; 
	}
	v -= part.size();
    }
    iFile.close();
}

void writeLog (std::string file){

    std::fstream fFile (file.c_str(), std::fstream::app );

    msd = 0;

    for( int v = 0 ; v < Nc; v++) msd += part[v].pos_msd[0]*part[v].pos_msd[0]+part[v].pos_msd[1]*part[v].pos_msd[1] + part[v].pos_msd[2]*part[v].pos_msd[2];
    msd = msd/Nc;
    
    fFile << step << "\t" << msd << std::endl;

    fFile.close();
}
