void writeSave (){

    oFile.open(savefile.c_str() );
    oFile << Nc << " " << Ns << std::endl;
    oFile << step << " " << pos_lambda << " " << ori_lambda << std::endl;
    oFile << lx << " " << ly << " " << lz << std::endl;

    for( int v = 0 ; v < Nc; v++){
            oFile <<  posx[v] << " " << posy[v] << " " << posz[v] << " "; 
            oFile <<  pos_msdx[v] << " " << pos_msdy[v] << " " << pos_msdz[v] << " "; 
            oFile <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << std::endl; 
    }
    oFile.close();
}


void writeConfig (){

    oFile.open(file_name.c_str() );
    oFile << Nc + Ns << std::endl;
    oFile << -lx_2 << " " << lx_2 << std::endl;
    oFile << -ly_2 << " " << ly_2 << std::endl;
    oFile << -lz_2 << " " << lz_2 << std::endl;

    for( int v = 0 ; v < Nc; v++){
        oFile << "PEAR(" << aspect << "," << kth << ") ";
        oFile <<  posx[v] - int(posx[v]/lx_2)*lx << " " << posy[v] - int(posy[v]/ly_2)*ly << " " << posz[v] - int(posz[v]/lz_2)*lz << " ";
        oFile << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " " << 1 << std::endl; 
    }
    oFile.close();
}

void read (){

    iFile.open(savefile.c_str() );

    if(!iFile){
        std::cerr << "Can not read " << savefile << "!" << std::endl;
        exit(-1);
    }
    iFile >> Nc >> Ns;

    if(Ns != 0){
	std::cerr << "Save/Config.dat does not contain monodisperse data! Please compile version for polydisperse data [make] :)" << std::endl;
	exit(0);
    }

    part.resize(Nc);
    
    iFile >> step >> pos_lambda >> ori_lambda;


    if( step < 100 && step != 0){
	std::cerr << "Old Save file version! No information for the step sizes in second line [step move_step rot_step]" << std::endl;
	std::cerr << "I try to convert it :)" << std::endl;
	
	iFile >> lz >> step;
	lx = lz;
	ly = lz;
    	
	pos_lambda = -1;
	ori_lambda = 0.02;

	std:: cout << step << " " << pos_lambda << " " << ori_lambda << " " << lz << std::endl;
    }else iFile >> lx >> ly >> lz;

    lx_2 = lx*0.5;
    ly_2 = ly*0.5;
    lz_2 = lz*0.5;

    Vbox=lx*ly*lz;

    rhoN = Nc/Vbox;
    rhoV = Vsys*rhoN;

    WPx = (int)(lx/rlistP); 
    WPy = (int)(ly/rlistP); 
    WPz = (int)(lz/rlistP); 

    wPx = WPx/lx;
    wPy = WPy/ly;
    wPz = WPz/lz;
	
    headP.resize(WPx*WPy*WPz); 
    usedCell.resize(WPx*WPy*WPz,false); 
    for( int i=0; i<headP.size(); i++) headP[i] = -1;

    linkP.resize(Nc*2); 
    for( int i=0; i<linkP.size(); i++) linkP[i] = -1;

	s_n.resize(WPx);

	s_nn[0] = 0;
	s_nn[1] = 1;
	s_nn[2] = WPx-1;
	s_nn[3] = 2;
	s_nn[4] = WPx-2;
	s_n[0]=s_nn;

	s_nn[0] = 1;
	s_nn[1] = 2;
	s_nn[2] = 0;
	s_nn[3] = 3;
	s_nn[4] = WPx-1;
	s_n[1]=s_nn;
	for( int i=2; i<s_n.size()-2; i++){ 
		s_nn[0] = i;
		s_nn[1] = i+1;
		s_nn[3] = i+2;
		s_nn[2] = i-1;
		s_nn[4] = i-2;

		s_n[i]=s_nn;
	}

	s_nn[0] = WPx-2;
	s_nn[1] = WPx-1;
	s_nn[2] = WPx-3;
	s_nn[3] = 0;
	s_nn[4] = WPx-4;

	s_n[WPx-2]=s_nn;

	s_nn[0] = WPx-1;
	s_nn[1] = 0;
	s_nn[2] = WPx-2;
	s_nn[3] = 1;
	s_nn[4] = WPx-3;

	s_n[WPx-1]=s_nn;

    for( int v = 0 ; v < Nc; v++){
	already[v] = false;

        iFile >> posx[v] >> posy[v] >> posz[v] >> pos_msdx[v] >> pos_msdy[v] >> pos_msdz[v] >> part[v].trans[0][2]  >> part[v].trans[1][2] >> part[v].trans[2][2]; 
        //std::cout << " " << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " " << part[v].pos_msd[0] << " " << part[v].pos_msd[1] << " " << part[v].pos_msd[2] << " " << part[v].trans[0][2]  << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << std::endl;

	if(posx[v] < 0 ) posx[v] += lx;
	if(posy[v] < 0 ) posy[v] += ly;
	if(posz[v] < 0 ) posz[v] += lz;

	if(step < 100) pos_msdx[v] =0;
	if(step < 100) pos_msdy[v] =0;
	if(step < 100) pos_msdz[v] =0;

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

	dist_orix[v] = distN*part[v].trans[0][2];
	dist_oriy[v] = distN*part[v].trans[1][2];
	dist_oriz[v] = distN*part[v].trans[2][2];

	NCell.resize(0);
	NPart.resize(0);

	dsx = posx[v]-dist_orix[v];
	if(dsx < 0) sx[v] = WPx-1;
	else{
		if(dsx > lx) sx[v] = 0;
		else sx[v] = dsx*wPx;
	}

	dsy = posy[v]-dist_oriy[v];
	if(dsy < 0) sy[v] = WPy-1;
	else{
		if(dsy > ly) sy[v] = 0;
		else sy[v] = dsy*wPy;
	}

	dsz = posz[v]-dist_oriz[v];
	if(dsz < 0) sz[v] = WPz-1;
	else{
		if(dsz > lz) sz[v] = 0;
		else sz[v] = dsz*wPz;
	}

	for (int iz=0; iz < 5 && !inside; iz++){
		kz = WPy*s_n[sz[v]][iz];
		for (int iy=0; iy < 5 && !inside; iy++){
			ky = WPx*(s_n[sy[v]][iy]+kz);
			for (int ix=0; ix < 5 && !inside; ix++){
				k = s_n[sx[v]][ix] + ky;
				if( usedCell[k] ) continue;

				vv = headP[k];
			    
				while(vv != -1){
					newvv = vv;
					if(newvv >= part.size() ) newvv -= part.size(); 
					if(!already[newvv]){ 
						NPart.push_back(newvv);
						already[newvv]=true;	
						Rx = posx[newvv]-posx[v]; 
						if (Rx > lx_2) Rx =- lx;
						else{
							if(Rx< -lx_2) Rx += lx;
						}
						Ry = posy[newvv]-posy[v]; 
						if (Ry > ly_2) Ry =- ly;
						else{
							if(Ry< -ly_2) Ry += ly;
						}
						Rz = posz[newvv]-posz[v]; 
						if (Rz > lz_2) Rz =- lz;
						else{
							if(Rz< -lz_2) Rz += lz;
						}
						if (overlapP ( Rx, Ry, Rz, part[v], part[newvv]) ){
							std::cerr << "There is still a problem (" << newvv << "," << v << ")!" << std::endl;

							std::cerr << "PEAR " << posx[newvv] << " " << posy[newvv] << " " << posz[newvv] << " " << part[newvv].trans[0][0] << " " << part[newvv].trans[1][0] << " " << part[newvv].trans[2][0]  << " " << part[newvv].trans[0][1] << " " << part[newvv].trans[1][1] << " " << part[newvv].trans[2][1] << " " <<  part[newvv].trans[0][2] << " " << part[newvv].trans[1][2] << " " << part[newvv].trans[2][2] << " 1" << std::endl;
							std::cerr << "PEAR " << posx[v] << " " << posy[v] << " " << posz[v] << " " << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " 1" << std::endl;
						}
					}
					vv = linkP[vv];
				}
				NCell.push_back(k);
				usedCell[k]=true;
			}
		}
	}

	dsxN = posx[v]+dist_orix[v];
	if(dsxN < 0) sNx[v] = WPx-1;
	else{
		if(dsxN > lz) sNx[v] = 0;
		else sNx[v] = dsxN*wPx;
	}

	dsyN = posy[v]+dist_oriy[v];
	if(dsyN < 0) sNy[v] = WPy-1;
	else{
		if(dsyN > lz) sNy[v] = 0;
		else sNy[v] = dsyN*wPy;
	}

	dszN = posz[v]+dist_oriz[v];
	if(dszN < 0) sNz[v] = WPz-1;
	else{
		if(dszN > lz) sNz[v] = 0;
		else sNz[v] = dszN*wPz;
	}

	if(sx[v] != sNx[v] || sy[v] != sNy[v] || sz[v] != sNz[v] ){
		for (int iz=0; iz < 5 && !inside; iz++){
			kz = WPy*s_n[sNz[v]][iz];
			for (int iy=0; iy < 5 && !inside; iy++){
				ky = WPx*(s_n[sNy[v]][iy]+kz);
				for (int ix=0; ix < 5 && !inside; ix++){

					k = s_n[sNx[v]][ix] + ky;
					if( usedCell[k] ) continue;

					vv = headP[k];
			    
					while(vv != -1){
						newvv = vv;
						if(newvv >= part.size() ) newvv -= part.size(); 
						if(!already[newvv]){ 
							NPart.push_back(newvv);
							already[newvv]=true;	
							Rx = posx[newvv]-posx[v]; 
							if (Rx > lx_2) Rx =- lx;
							else{
								if(Rx< -lx_2) Rx += lx;
							}
							Ry = posy[newvv]-posy[v]; 
							if (Ry > ly_2) Ry =- ly;
							else{
								if(Ry< -ly_2) Ry += ly;
							}
							Rz = posz[newvv]-posz[v]; 
							if (Rz > lz_2) Rz =- lz;
							else{
								if(Rz< -lz_2) Rz += lz;
							}
							if (overlapP ( Rx, Ry, Rz, part[v], part[newvv]) ){
								std::cerr << "There is still a problem (" << newvv << "," << v << ")!" << std::endl;

							std::cerr << "PEAR " << posx[newvv] << " " << posy[newvv] << " " << posz[newvv] << " " << part[newvv].trans[0][0] << " " << part[newvv].trans[1][0] << " " << part[newvv].trans[2][0]  << " " << part[newvv].trans[0][1] << " " << part[newvv].trans[1][1] << " " << part[newvv].trans[2][1] << " " <<  part[newvv].trans[0][2] << " " << part[newvv].trans[1][2] << " " << part[newvv].trans[2][2] << " 1" << std::endl;
							std::cerr << "PEAR " << posx[v] << " " << posy[v] << " " << posz[v] << " " << part[v].trans[0][0] << " " << part[v].trans[1][0] << " " << part[v].trans[2][0]  << " " << part[v].trans[0][1] << " " << part[v].trans[1][1] << " " << part[v].trans[2][1] << " " <<  part[v].trans[0][2] << " " << part[v].trans[1][2] << " " << part[v].trans[2][2] << " 1" << std::endl;
							}
						}
						vv = linkP[vv];
					}
				}
			}
		}
	}

	for( vv = 0; vv < NPart.size(); vv++) already[NPart[vv]]=false;
	for( vv = 0; vv < NCell.size(); vv++) usedCell[NCell[vv]]=false;

	int k = sx[v] + WPx*(sy[v]+WPy*sz[v]); 

	cell[v] = k;
	if(headP[k]==-1) headP[k] = v;
	else{
	   linkP[v] = headP[k];
	   headP[k] = v; 
	}

	k = sNx[v] + WPx*(sNy[v]+WPy*sNz[v]); 

	cellN[v] = k;
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

void writeLog (){

    fFile.open(logfile.c_str(), std::fstream::app );

    msd = 0;

    for( int v = 0 ; v < Nc; v++) msd += pos_msdx[v]*pos_msdx[v]+pos_msdy[v]*pos_msdy[v] + pos_msdz[v]*pos_msdz[v];
    msd = msd/Nc;
    
    fFile << step << "\t" << msd << std::endl;

    fFile.close();
}
