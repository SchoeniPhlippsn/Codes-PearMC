void system::write (std::string file, bool save){

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

    for( v = 0 ; v < Nc; v++){

        if(!save){ 

            oFile << "PEAR(" << aspect << "," << kth << ") ";
            oFile <<  part[v].pos[0] - int(part[v].pos[0]/l_2[0])*l[0] << " " << part[v].pos[1] - int(part[v].pos[1]/l_2[1])*l[1] << " " << part[v].pos[2] - int(part[v].pos[2]/l_2[2])*l[2] << " ";
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
            oFile <<  part[v].pos[0] - int(part[v].pos[0]/l_2[0])*l[0] << " " << part[v].pos[1] - int(part[v].pos[1]/l_2[1])*l[1] << " " << part[v].pos[2] - int(part[v].pos[2]/l_2[2])*l[2] << " ";
            oFile << 2*rsphere << " 0 0 0 " << 2*rsphere << " 0 0 0 " << 2*rsphere << " " << v << std::endl;
        }else{
            oFile <<  part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " "; 
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
    
    iFile >> step >> pos_lambda >> ori_lambda;

    if( pos_lambda > 2 ){
	std::cerr << "Old Save file version! No information for the step sizes in second line [step move_step rot_step]" << std::endl;
	exit(0);
    }

    iFile >> l[0] >> l[1] >> l[2];

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

    for( v = 0 ; v < Nc; v++){
        part[v].pos.resize(3);
        part[v].pos_msd.resize(3);
        part[v].ori.resize(3);
	part[v].already = false;

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2] >> part[v].pos_msd[0] >> part[v].pos_msd[1] >> part[v].pos_msd[2] >> part[v].ori[0]  >> part[v].ori[1] >> part[v].ori[2]; 

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

	NCell.resize(0);
	NPart.resize(0);

	dsx = part[v].pos[0]-distN*part[v].ori[0];
	if(dsx < 0){
		sx[0] = WP[0]-1;
		sx[1] = sx[0] -1;
		sx[2] = 0;
	}else{
		if(dsx > l[0]){
			sx[0] = 0;
			sx[1] = WP[0]-1;
			sx[2] = 1;
		}else{
			sx[0] = dsx*wP[0];

			sx[1] = sx[0]-1;
			
			if(sx[1] < 0 ) sx[1] = WP[0]-1;

			sx[2] = sx[0]+1;
			if(sx[2]> WP[0]-1) sx[2] = 0;
		}
	}

	dsy = part[v].pos[1]-distN*part[v].ori[1];
	if(dsy < 0){
		sy[0] = WP[1]-1;
		sy[1] = sy[0] -1;
		sy[2] = 0;
	}else{
		if(dsy > l[1]){
			sy[0] = 0;
			sy[1] = WP[1]-1;
			sy[2] = 1;
		}else{
			sy[0] = dsy*wP[1];

			sy[1] = sy[0]-1;
			if(sy[1] < 0 ) sy[1] = WP[1]-1;

			sy[2] = sy[0]+1;
			if(sy[2]> WP[1]-1) sy[2] = 0;
		}
	}

	dsz = part[v].pos[2]-distN*part[v].ori[2];
	if(dsz < 0){
		sz[0] = WP[2]-1;
		sz[1] = sz[0] -1;
		sz[2] = 0;
	}else{
		if(dsz > l[2]){
			sz[0] = 0;
			sz[1] = WP[2]-1;
			sz[2] = 1;
		}else{
			sz[0] = dsz*wP[2];

			sz[1] = sz[0]-1;
			if(sz[1] < 0 ) sz[1] = WP[2]-1;

			sz[2] = sz[0]+1;
			if(sz[2]> WP[2]-1) sz[2] = 0;
		}
	}

	for (int ix=0; ix < 3; ix++){

		for (int iy=0; iy < 3; iy++){

			for (int iz=0; iz < 3; iz++){
				int k = sx[ix] + WP[0]*(sy[iy]+WP[1]*sz[iz]);
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

	dsxN = part[v].pos[0]+distN*part[v].ori[0];
	if(dsxN < 0){
		sxN[0] = WP[0]-1;
		sxN[1] = sxN[0] -1;
		sxN[2] = 0;
	}else{
		if(dsxN > l[0]){
			sxN[0] = 0;
			sxN[1] = WP[0]-1;
			sxN[2] = 1;
		}else{
			sxN[0] = dsxN*wP[0];

			sxN[1] = sxN[0]-1;
			if(sxN[1] < 0 ) sxN[1] = WP[0]-1;

			sxN[2] = sxN[0]+1;
			if(sxN[2]> WP[0]-1) sxN[2] = 0;
		}
	}

	dsyN = part[v].pos[1]+distN*part[v].ori[1];
	if(dsyN < 0){
		syN[0] = WP[1]-1;
		syN[1] = syN[0] -1;
		syN[2] = 0;
	}else{
		if(dsyN > l[1]){
			syN[0] = 0;
			syN[1] = WP[1]-1;
			syN[2] = 1;
		}else{
			syN[0] = dsyN*wP[1];

			syN[1] = syN[0]-1;
			if(syN[1] < 0 ) syN[1] = WP[1]-1;

			syN[2] = syN[0]+1;
			if(syN[2]> WP[1]-1) syN[2] = 0;
		}
	}

	dszN = part[v].pos[1]+distN*part[v].ori[1];
	if(dszN < 0){
		szN[0] = WP[2]-1;
		szN[1] = szN[0] -1;
		szN[2] = 0;
	}else{
		if(dszN > l[2]){
			szN[0] = 0;
			szN[1] = WP[2]-1;
			szN[2] = 1;
		}else{
			szN[0] = dszN*wP[2];

			szN[1] = szN[0]-1;
			if(szN[1] < 0 ) szN[1] = WP[2]-1;

			szN[2] = szN[0]+1;
			if(szN[2]> WP[2]-1) szN[2] = 0;
		}
	}

	if(sx[1] != sxN[1] || sy[1] != syN[1] || sz[1] != szN[1] ){
		for (int ix=0; ix < 3; ix++){

			for (int iy=0; iy < 3; iy++){

				for (int iz=0; iz < 3; iz++){
					int k = sxN[ix] + WP[0]*(syN[iy]+WP[1]*szN[iz]);
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

	int k = sx[0] + WP[0]*(sy[0]+WP[1]*sz[0]); 

	part[v].cell = k;
	if(headP[k]==-1) headP[k] = v;
	else{
	   linkP[v] = headP[k];
	   headP[k] = v; 
	}

	k = sxN[0] + WP[0]*(syN[0]+WP[1]*szN[0]); 

	part[v].cellN = k;
	v += part.size();
	if(headP[k]==-1) headP[k] = v;
	else{
	   linkP[v] = headP[k];
	   headP[k] = v; 
	}
	v -= part.size();
    }


    WS[0] = (int)(l[0]/(rlistS)); 
    WS[1] = (int)(l[1]/(rlistS)); 
    WS[2] = (int)(l[2]/(rlistS)); 

    wS[0] = WS[0]/l[0];
    wS[1] = WS[1]/l[1];
    wS[2] = WS[2]/l[2];

    headS.resize(WS[0]*WS[1]*WS[2]); 
    for( int i=0; i<headS.size(); i++) headS[i] = -1;

    headPS.resize(WP[0]*WP[1]*WP[2]); 
    for( int i=0; i<headPS.size(); i++) headPS[i] = -1;


    for( v = Nc ; v < part.size(); v++){
        part[v].pos.resize(3);
        part[v].pos_msd.resize(3);
        part[v].ori.resize(3);

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2]  >> part[v].pos_msd[0] >> part[v].pos_msd[1] >> part[v].pos_msd[2]>> rsphere; 

	if(part[v].pos[0] < 0 ) part[v].pos[0] += l[0];
	if(part[v].pos[1] < 0 ) part[v].pos[1] += l[1];
	if(part[v].pos[2] < 0 ) part[v].pos[2] += l[2];

	if(step < 100) part[v].pos_msd[0] =0;
	if(step < 100) part[v].pos_msd[1] =0;
	if(step < 100) part[v].pos_msd[2] =0;

	NPart.resize(0);

	sx[0] = part[v].pos[0]*wS[0];
	if(sx[0] == WS[0]) sx[0] = WS[0]-1;
	sy[0] = part[v].pos[1]*wS[1];
	if(sy[0] == WS[1]) sy[0] = WS[1]-1;
	sz[0] = part[v].pos[2]*wS[2];
	if(sz[0] == WS[2]) sz[0] = WS[2]-1;

	int kN = sx[0] + WS[0]*(sy[0]+WS[1]*sz[0]);
		
	if(headS[kN] != -1){ 
		std::cerr << "There is still a problem:" << vv << " and " << v << " in same cell!!!" << std::endl;
		exit(0);
	}

	
	sxN[0] = part[v].pos[0]*wP[0];
	if(sxN[0] == WP[0]) sxN[0] = WP[0]-1;
	syN[0] = part[v].pos[1]*wP[1];
	if(syN[0] == WP[1]) syN[0] = WP[1]-1;
	szN[0] = part[v].pos[2]*wP[2];
	if(szN[0] == WP[2]) syN[0] = WP[2]-1;

	sxN[1] = sxN[0]-1;
	if(sxN[1] < 0 ) sxN[1] = WP[0]-1;

	sxN[2] = sxN[0]+1;
	if(sxN[2]> WP[0]-1) sxN[2] = 0;

	syN[1] = syN[0]-1;
	if(syN[1] < 0 ) syN[1] = WP[1]-1;

	syN[2] = syN[0]+1;
	if(syN[2]> WP[1]-1) syN[2] = 0;

	szN[1] = szN[0]-1;
	if(szN[1] < 0 ) szN[1] = WP[2]-1;

	szN[2] = szN[0]+1;
	if(szN[2]> WP[2]-1) szN[2] = 0;

	for (int ix=0; ix < 3; ix++){

		for (int iy=0; iy < 3; iy++){
			for (int iz=0; iz < 3; iz++){

				int k = sxN[ix] + WP[0]*(syN[iy]+WP[1]*szN[iz]);

				vv = headP[k];
			    
				while(vv != -1){
					newvv = vv;
					if(newvv >= part.size() ) newvv -= part.size(); 
					if(!part[newvv].already){ 
						NPart.push_back(newvv);
						part[newvv].already=true;	
						if (overlapPSPH ( part[newvv], part[v], l) ){ 
							std::cerr << "There is still a problem (" << newvv << "," << v << ")!" << std::endl;

							std::cerr << "PEAR " << part[newvv].pos[0] << " " << part[newvv].pos[1] << " " << part[newvv].pos[1] << " " << part[newvv].trans[0][0] << " " << part[newvv].trans[1][0] << " " << part[newvv].trans[2][0]  << " " << part[newvv].trans[0][1] << " " << part[newvv].trans[1][1] << " " << part[newvv].trans[2][1] << " " <<  part[newvv].trans[0][2] << " " << part[newvv].trans[1][2] << " " << part[newvv].trans[2][2] << " 1" << std::endl;

							std::cerr << "SPHERE " << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " " << 2*rsphere << " 0 0 0 " << 2*rsphere << " 0 0 0 " << 2*rsphere << " 1" << std::endl;
						}
					}
					vv = linkP[vv];
				}
			}
		}
	}
	for( vv = 0; vv < NPart.size(); vv++) part[NPart[vv]].already=false;


	bx[0] =0;
	bx[1] =0;
	bx[2] =0;
	bx[3] =0;
	bx[4] =0;

	sx[1] = sx[0]-1;
	if(sx[1] < 0 ){ 
		sx[1] = WS[0]-1;
		bx[1]=-1;
		bx[3]=-1;
	}

	sx[2] = sx[0]+1;
	if(sx[2]> WS[0]-1){ 
		sx[2] = 0;
		bx[2]=1;
		bx[4]=1;
	}

	sx[3] = sx[1]-1;
	if(sx[3] < 0 ){ 
		sx[3] = WS[0]-1;
		bx[3]=-1;
	}

	sx[4] = sx[2]+1;
	if(sx[4]> WS[0]-1){ 
		sx[4] = 0;
		bx[4]=1;
	}

	by[0] =0;
	by[1] =0;
	by[2] =0;
	by[3] =0;
	by[4] =0;

	sy[1] = sy[0]-1;
	if(sy[1] < 0 ){ 
		sy[1] = WS[1]-1;
		by[1]=-1;
		by[3]=-1;
	}

	sy[2] = sy[0]+1;
	if(sy[2]> WS[1]-1){ 
		sy[2] = 0;
		by[2]=1;
		by[4]=1;
	}

	sy[3] = sy[1]-1;
	if(sy[3] < 0 ){ 
		sy[3] = WS[1]-1;
		by[3]=-1;
	}

	sy[4] = sy[2]+1;
	if(sy[4]> WS[1]-1){ 
		sy[4] = 0;
		by[4]=1;
	}

	bz[0] =0;
	bz[1] =0;
	bz[2] =0;
	bz[3] =0;
	bz[4] =0;

	sz[1] = sz[0]-1;
	if(sz[1] < 0 ){ 
		sz[1] = WS[2]-1;
		bz[1]=-1;
		bz[3]=-1;
	}

	sz[2] = sz[0]+1;
	if(sz[2]> WS[2]-1){ 
		sz[2] = 0;
		bz[2]=1;
		bz[4]=1;
	}

	sz[3] = sz[1]-1;
	if(sz[3] < 0 ){ 
		sz[3] = WS[2]-1;
		bz[3]=-1;
	}

	sz[4] = sz[2]+1;
	if(sz[4]> WS[2]-1){ 
		sz[4] = 0;
		bz[4]=1;
	}

	for (int ix=0; ix < 5; ix++){

		for (int iy=0; iy < 5; iy++){
			for (int iz=0; iz < 5; iz++){

				int k = sx[ix] + WS[0]*(sy[iy]+WS[1]*sz[iz]);

				vv = headS[k];
		    
				if(vv != -1){
					if (overlapSPH( part[v], part[vv], l) ){
						std::cerr << "There is still a problem (" << vv << "," << v << ")!" << std::endl;
						std::vector<double> RTest = DeltaR(part[vv].pos,part[v].pos,l);
						double testing = scal_p(RTest,RTest);
						std::cerr << testing << " " << rcut_SPH << " " << rsphere << std::endl;
					}
				}
			}
		}
	}

	int k = sx[0] + WS[0]*(sy[0]+WS[1]*sz[0]);

	part[v].cell = k;
	headS[k] = v;

	 k = sxN[0] + WP[0]*(syN[0]+WP[1]*szN[0]);

	part[v].cellN = k;
	if(headPS[k]==-1) headPS[k] = v;
	else{
		linkP[v] = headPS[k];
		headPS[k] = v; 
	}
    }
    iFile.close();
}

void system::writeLog (std::string file){

    std::fstream fFile (file.c_str(), std::fstream::app );

    msd = 0;

    for( v = 0 ; v < Nc; v++) msd += part[v].pos_msd[0]*part[v].pos_msd[0]+part[v].pos_msd[1]*part[v].pos_msd[1] + part[v].pos_msd[2]*part[v].pos_msd[2];
    msd = msd/Nc;
    
    fFile << step << msd;

    msd = 0;
    for( v = Nc ; v < part.size() ; v++) msd += part[v].pos_msd[0]*part[v].pos_msd[0]+part[v].pos_msd[1]*part[v].pos_msd[1]+ part[v].pos_msd[2]*part[v].pos_msd[2];
    msd = msd/Ns;

    fFile <<  "\t" << msd << std::endl;

    fFile.close();
}
