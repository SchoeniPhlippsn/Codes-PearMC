void Move_step(){
	for(int v = 0;  v < N; v++){
        	randP = uni(gen)*N;
		
        	if (randP==N){
            		v--;
            		continue;
		}
        	MovedParticle = part[randP];         

        	inside = false;

		phi = 2*M_PI*uni(gen);
		theta = acos(2*uni(gen) - 1);				
		cosphi = uni(gen);

		costheta = cos(theta)*sin(phi)*cosphi;
		sintheta = sin(theta)*sin(phi)*cosphi;
		cosphi = cos(phi)*cosphi;

	    if( 2*uni(gen) < 1 ){
			MovedParticle.pos[0] = MovedParticle.pos[0] + costheta*pos_lambda;
			MovedParticle.pos[1] = MovedParticle.pos[1] + sintheta*pos_lambda;
			MovedParticle.pos[2] = MovedParticle.pos[2] + cosphi*pos_lambda;

			MovedParticle.pos_msd[0] += (MovedParticle.pos[0]-part[randP].pos[0]);
			MovedParticle.pos_msd[1] += (MovedParticle.pos[1]-part[randP].pos[1]);
			MovedParticle.pos_msd[2] += (MovedParticle.pos[2]-part[randP].pos[2]);

	    }else{
			MovedParticle.trans[0][2] = MovedParticle.trans[0][2] + costheta*ori_lambda;
			MovedParticle.trans[1][2] = MovedParticle.trans[1][2] + sintheta*ori_lambda;
			MovedParticle.trans[2][2] = MovedParticle.trans[2][2] + cosphi*ori_lambda;

			norm = MovedParticle.trans[0][2]*MovedParticle.trans[0][2] + MovedParticle.trans[1][2]*MovedParticle.trans[1][2] + MovedParticle.trans[2][2]*MovedParticle.trans[2][2];
			norm = 1.0/sqrt(norm);

			MovedParticle.trans[0][2] *= norm;
			MovedParticle.trans[1][2] *= norm;
			MovedParticle.trans[2][2] *= norm;

			sqrtz = 1-MovedParticle.trans[2][2]*MovedParticle.trans[2][2];
			if(sqrtz > 1e-4){
				sqrtz = sqrt(sqrtz);

				MovedParticle.trans[2][0] = -sqrtz;

				sqrtz = 1/sqrtz;
				MovedParticle.trans[0][1] = -MovedParticle.trans[1][2]*sqrtz;
				MovedParticle.trans[1][1] = MovedParticle.trans[0][2]*sqrtz;

				MovedParticle.trans[0][0] = MovedParticle.trans[2][2]*MovedParticle.trans[1][1];
				MovedParticle.trans[1][0] = -MovedParticle.trans[2][2]*MovedParticle.trans[0][1];

			}else{
				MovedParticle.trans[0][0] = -1;
				MovedParticle.trans[0][1] = 0;
				MovedParticle.trans[0][2] = 0;
				
				MovedParticle.trans[1][0] = 0;
				MovedParticle.trans[1][2] = 0;


				MovedParticle.trans[2][0] = 0;

				if(MovedParticle.trans[2][2]>0){
					MovedParticle.trans[1][1] = -1;
					MovedParticle.trans[2][2] = 1;
				}else{
					MovedParticle.trans[1][1] = 1;
					MovedParticle.trans[2][2] = -1;
				}
			}

			MovedParticle.dist_ori[0] = distN*MovedParticle.trans[0][2];
			MovedParticle.dist_ori[1] = distN*MovedParticle.trans[1][2];
			MovedParticle.dist_ori[2] = distN*MovedParticle.trans[2][2];

		}
		pear_mesh.UpdateTrans(id[0], MovedParticle.trans);


		NCell.resize(0);
		NPart.resize(0);

		vv = headP[MovedParticle.cell];
	    
		while(vv != -1 && !inside){
			newvv = vv;
			if(newvv >= N ) newvv -= N; 
			if(randP != newvv && !part[newvv].already){ 
				NPart.push_back(newvv);
				part[newvv].already=true;	

				R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
				R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
				R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
				if (R[0] > l_2[0])  R[0] -= l[0];
				else if (R[0] < -l_2[0])  R[0] += l[0];
				if (R[1] > l_2[1])  R[1] -= l[1];
				else if (R[1] < -l_2[1])  R[1] += l[1];
				if (R[2] > l_2[2])  R[2] -= l[2];
				else if (R[2] < -l_2[2])  R[2] += l[2];
				overlapPear(); 
			}

			vv = linkP[vv];
		}
		NCell.push_back(MovedParticle.cell);
		usedCell[MovedParticle.cell]=true;

		if( MovedParticle.cellN != MovedParticle.cell && !inside){

			vv = headP[MovedParticle.cellN];
		    
			while(vv != -1 && !inside){
				newvv = vv;
				if(newvv >= N ) newvv -= N; 
				if(randP != newvv && !part[newvv].already){ 
					NPart.push_back(newvv);
					part[newvv].already=true;
					R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
					R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
					R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
					if (R[0] > l_2[0])  R[0] -= l[0];
					else if (R[0] < -l_2[0])  R[0] += l[0];
					if (R[1] > l_2[1])  R[1] -= l[1];
					else if (R[1] < -l_2[1])  R[1] += l[1];
					if (R[2] > l_2[2])  R[2] -= l[2];
					else if (R[2] < -l_2[2])  R[2] += l[2];
					overlapPear(); 
				}

				vv = linkP[vv];
			}
			NCell.push_back(MovedParticle.cellN);
			usedCell[MovedParticle.cellN]=true;
		}


		if(!inside){
			if( MovedParticle.pos[0] > l[0] ) MovedParticle.pos[0] -= l[0];
			else if( MovedParticle.pos[0] < 0 ) MovedParticle.pos[0] += l[0];
		
			if( MovedParticle.pos[1] > l[1] ) MovedParticle.pos[1] -= l[1];
			else if( MovedParticle.pos[1] < 0 ) MovedParticle.pos[1] += l[1];

			if( MovedParticle.pos[2] > l[2] ) MovedParticle.pos[2] -= l[2];
			else if( MovedParticle.pos[2] < 0 ) MovedParticle.pos[2] += l[2];

			dsx = MovedParticle.pos[0]-MovedParticle.dist_ori[0];
			dsy = MovedParticle.pos[1]-MovedParticle.dist_ori[1];
			dsz = MovedParticle.pos[2]-MovedParticle.dist_ori[2];

			in_test = false;
			if(dsx < 0){
				MovedParticle.s[0] = WP[0]-1;
				in_test = true;
			}else{
				if(dsx > l[0]){
					MovedParticle.s[0] = 0;
					in_test = true;
				}else{
					MovedParticle.s[0] = dsx*wP[0];
					if( MovedParticle.s[0] >= WP[0]-3 || MovedParticle.s[0] <= 2) in_test = true;
				}
			}
			

			if(dsy < 0){
				MovedParticle.s[1] = WP[1]-1;
				in_test = true;
			}else{
				if(dsy > l[1]){
					MovedParticle.s[1] = 0;
					in_test = true;
				}else{
					MovedParticle.s[1] = dsy*wP[1];
					if( MovedParticle.s[1] >= WP[1]-3 || MovedParticle.s[1] <= 2) in_test = true;
				}
			}

			if(dsz < 0){
				MovedParticle.s[2] = WP[2]-1;
				in_test = true;
			}else{
				if(dsz > l[2]){
					MovedParticle.s[2] = 0;
					in_test = true;
				}else{
					MovedParticle.s[2] = dsz*wP[2];
					if( MovedParticle.s[2] >= WP[2]-3 || MovedParticle.s[2] <= 2) in_test = true;
				}
			}

			if(in_test){
				for (int iz=0; iz < 5 && !inside; iz++){
					kz = WP[1]*s_n[MovedParticle.s[2]][iz];
					for (int iy=0; iy < 5 && !inside; iy++){
						ky = WP[0]*(s_n[MovedParticle.s[1]][iy]+kz);
						for (int ix=0; ix < 5 && !inside; ix++){
							k = s_n[MovedParticle.s[0]][ix] + ky;
							if( usedCell[k] ) continue;
				
							vv = headP[k];

							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= N ) newvv -= N; 
								if(!part[newvv].already){ 
									NPart.push_back(newvv);
									part[newvv].already=true;	
									
									R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
									R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
									R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
									if (R[0] > l_2[0])  R[0] -= l[0];
									else if (R[0] < -l_2[0])  R[0] += l[0];
									if (R[1] > l_2[1])  R[1] -= l[1];
									else if (R[1] < -l_2[1])  R[1] += l[1];
									if (R[2] > l_2[2])  R[2] -= l[2];
									else if (R[2] < -l_2[2])  R[2] += l[2];
									overlapPear(); 
								}
								vv = linkP[vv];
							}
							NCell.push_back(k);
							usedCell[k]=true;
						}
					}
				}
			}else{
				for (int iz=0; iz < 5 && !inside; iz++){
					kz = WP[1]*s_n[MovedParticle.s[2]][iz];
					for (int iy=0; iy < 5 && !inside; iy++){
						ky = WP[0]*(s_n[MovedParticle.s[1]][iy]+kz);
						for (int ix=0; ix < 5 && !inside; ix++){
							k = s_n[MovedParticle.s[0]][ix] + ky;
							if( usedCell[k] ) continue;
				
							vv = headP[k];

							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= N ) newvv -= N; 
								if(!part[newvv].already){ 
									NPart.push_back(newvv);
									part[newvv].already=true;	
									
									R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
									R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
									R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
									overlapPear(); 
								}
								vv = linkP[vv];
							}
							NCell.push_back(k);
							usedCell[k]=true;
						}
					}
				}
			}
		}
	
		if(!inside){

			dsxN = MovedParticle.pos[0]+MovedParticle.dist_ori[0];

			dsyN = MovedParticle.pos[1]+MovedParticle.dist_ori[1];

			dszN = MovedParticle.pos[2]+MovedParticle.dist_ori[2];

			in_test = false;
			if(dsxN < 0){
				MovedParticle.sN[0] = WP[0]-1;
				in_test = true;
			}else{
				if(dsxN > l[0]){
					MovedParticle.sN[0] = 0;
					in_test = true;
				}else{
					MovedParticle.sN[0] = dsxN*wP[0];
					if( MovedParticle.sN[0] >= WP[0]-3 || MovedParticle.sN[0] <= 2) in_test = true;
				}
			}

			if(dsyN < 0){
				MovedParticle.sN[1] = WP[1]-1;
				in_test = true;
			}else{
				if(dsyN > l[1]){
					MovedParticle.sN[1] = 0;
					in_test = true;
				}else{
					MovedParticle.sN[1] = dsyN*wP[1];
					if( MovedParticle.sN[1] >= WP[1]-3 || MovedParticle.sN[1] <= 2) in_test = true;
				}
			}

			if(dszN < 0){
				MovedParticle.sN[2] = WP[2]-1;
				in_test = true;
			}else{
				if(dszN > l[2]){
					MovedParticle.sN[2] = 0;
					in_test = true;
				}else{
					MovedParticle.sN[2] = dszN*wP[2];
					if( MovedParticle.sN[2] >= WP[2]-3 || MovedParticle.sN[2] <= 2) in_test = true;
				}
			}
			MovedParticle.sN[2]=MovedParticle.sN[2];

			if(MovedParticle.s[0] != MovedParticle.sN[0] || MovedParticle.s[1] != MovedParticle.sN[1] || MovedParticle.s[2] != MovedParticle.sN[2] ){

				if(in_test){
					for (int iz=0; iz < 5 && !inside; iz++){
						kz = WP[1]*s_n[MovedParticle.sN[2]][iz];
						for (int iy=0; iy < 5 && !inside; iy++){
							ky = WP[0]*(s_n[MovedParticle.sN[1]][iy]+kz);
							for (int ix=0; ix < 5 && !inside; ix++){

								k = s_n[MovedParticle.sN[0]][ix] + ky;
								if( usedCell[k] ) continue;
				
								vv = headP[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= N ) newvv -= N; 
									if(!part[newvv].already){ 
										NPart.push_back(newvv);
										part[newvv].already=true;	

										R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
										if (R[0] > l_2[0])  R[0] -= l[0];
										else if (R[0] < -l_2[0])  R[0] += l[0];
										if (R[1] > l_2[1])  R[1] -= l[1];
										else if (R[1] < -l_2[1])  R[1] += l[1];
										if (R[2] > l_2[2])  R[2] -= l[2];
										else if (R[2] < -l_2[2])  R[2] += l[2];
										overlapPear(); 
									}
									vv = linkP[vv];
								}
							}
						}
					}
				}else{
					for (int iz=0; iz < 5 && !inside; iz++){
						kz = WP[1]*s_n[MovedParticle.sN[2]][iz];
						for (int iy=0; iy < 5 && !inside; iy++){
							ky = WP[0]*(s_n[MovedParticle.sN[1]][iy]+kz);
							for (int ix=0; ix < 5 && !inside; ix++){

								k = s_n[MovedParticle.sN[0]][ix] + ky;
								if( usedCell[k] ) continue;
				
								vv = headP[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= N ) newvv -= N; 
									if(!part[newvv].already){ 
										NPart.push_back(newvv);
										part[newvv].already=true;	

										R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
										overlapPear(); 
									}
									vv = linkP[vv];
								}
							}
						}
					}
				}
			}
		}

		for( vv = 0; vv < NPart.size(); vv++) part[NPart[vv]].already=false;
		for( vv = 0; vv < NCell.size(); vv++) usedCell[NCell[vv]]=false;



		if(!inside){
	            	acc++;
			k = MovedParticle.s[0] + WP[0]*(MovedParticle.s[1]+WP[1]*MovedParticle.s[2]); 

			if(MovedParticle.cell != k){
				if(headP[MovedParticle.cell] != randP ){
					int vvv = headP[MovedParticle.cell];
					vv = linkP[vvv];
					while(vv!=randP){ 
						vvv = vv;
						vv = linkP[vvv];
					}
			
					linkP[vvv] = linkP[randP];
				}else{
					headP[MovedParticle.cell] = linkP[randP];
				}
	    
				linkP[randP] = headP[k];
				headP[k] = randP;
				MovedParticle.cell = k;
			}

			k = MovedParticle.sN[0] + WP[0]*(MovedParticle.sN[1]+WP[1]*MovedParticle.sN[2]); 

			if(MovedParticle.cellN != k){
				if(headP[MovedParticle.cellN] != randP+N ){
					int vvv = headP[MovedParticle.cellN];
					vv = linkP[vvv];
					while(vv!=randP+N){ 
						vvv = vv;
						vv = linkP[vvv];
					}
			
					linkP[vvv] = linkP[randP+N];
				}else{
					headP[MovedParticle.cellN] = linkP[randP+N];
				}
	    
				linkP[randP+N] = headP[k];
				headP[k] = randP+N;
				MovedParticle.cellN = k;
			}
			part[randP] = MovedParticle;
        	}
	}
}



void Compression_step(){	

    if(Compressing){
	Vn = log(Vbox)*vproc;
	Vn = exp(Vn);
	powVn = pow(Vn,1.0/3.0);

    	ln[0] = powVn/l[0];
    	ln[1] = powVn/l[1];
    	ln[2] = powVn/l[2];

    }else{
	Vn = (1+uni(gen)*vproc);
	powVn = sqrt(Vn);

	Case = 2*uni(gen);
	if(Case==0) powVn = 1/powVn;
	else Vn = 1/Vn;

	Case = 3*uni(gen);
	switch(Case){
		case 0:
		ln[0]=Vn;
		ln[1]=powVn;
		ln[2]=powVn;
		break;
		case 1:
		ln[0]=powVn;
		ln[1]=Vn;
		ln[2]=powVn;
		break;
		case 2:
		ln[0]=powVn;
		ln[1]=powVn;
		ln[2]=Vn;
		break;
	}
   }

    l[0] *= ln[0];
    l[1] *= ln[1];
    l[2] *= ln[2];

    l_2[0] = l[0]*0.5;
    l_2[1] = l[1]*0.5;
    l_2[2] = l[2]*0.5;


    for( int i=0; i < N; i++){
            part[i].pos[0] *= ln[0];
            part[i].pos[1] *= ln[1];
            part[i].pos[2] *= ln[2];
    }
    RenewList();
    inside = false;

    for( int v=0; v < N-1 && !inside; v++){
	MovedParticle=part[v];

	pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

	NCell.resize(0);
	NPart.resize(0);

	k = MovedParticle.cell;

	vv = headP[k];

	while(vv != -1 && !inside){
		newvv = vv;
		if(newvv >= N ) newvv -= N; 
		if(v < newvv && !part[newvv].already){ 
			NPart.push_back(newvv);
			part[newvv].already=true;

			R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
			R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
			R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
			if (R[0] > l_2[0])  R[0] -= l[0];
			else if (R[0] < -l_2[0])  R[0] += l[0];
			if (R[1] > l_2[1])  R[1] -= l[1];
			else if (R[1] < -l_2[1])  R[1] += l[1];
			if (R[2] > l_2[2])  R[2] -= l[2];
			else if (R[2] < -l_2[2])  R[2] += l[2];
			overlapPear(); 
		}

		vv = linkP[vv];
	}
	NCell.push_back(k);
	usedCell[k]=true;

	if( MovedParticle.cellN != MovedParticle.cell && !inside){
		k = MovedParticle.cellN;

		vv = headP[k];
	    
		while(vv != -1 && !inside){
			newvv = vv;
			if(newvv >= N ) newvv -= N; 
			if(v < newvv && !part[newvv].already){
				NPart.push_back(newvv);
				part[newvv].already=true;
				R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
				R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
				R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
				if (R[0] > l_2[0])  R[0] -= l[0];
				else if (R[0] < -l_2[0])  R[0] += l[0];
				if (R[1] > l_2[1])  R[1] -= l[1];
				else if (R[1] < -l_2[1])  R[1] += l[1];
				if (R[2] > l_2[2])  R[2] -= l[2];
				else if (R[2] < -l_2[2])  R[2] += l[2];
				overlapPear(); 
			}
			vv = linkP[vv];
		}
		NCell.push_back(k);
		usedCell[k]=true;
	}

	if(!inside){

		in_test = false;
		if( MovedParticle.s[0] >= WP[0]-3 || MovedParticle.s[0] <= 2) in_test = true;
		if( MovedParticle.s[1] >= WP[1]-3 || MovedParticle.s[1] <= 2) in_test = true;
		if( MovedParticle.s[2] >= WP[2]-3 || MovedParticle.s[2] <= 2) in_test = true;

		if(in_test){
			for (int iz=0; iz < 5 && !inside; iz++){
				kz = WP[1]*s_n[MovedParticle.s[2]][iz];
				for (int iy=0; iy < 5 && !inside; iy++){
					ky = WP[0]*(s_n[MovedParticle.s[1]][iy]+kz);
					for (int ix=0; ix < 5 && !inside; ix++){
						k = s_n[MovedParticle.s[0]][ix] + ky;
						if( usedCell[k] ) continue;
		
						vv = headP[k];

						while(vv != -1 && !inside){
							newvv = vv;
							if(newvv >= N ) newvv -= N; 
							if(v < newvv && !part[newvv].already){
								NPart.push_back(newvv);
								part[newvv].already=true;

								R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
								R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
								R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
								if (R[0] > l_2[0])  R[0] -= l[0];
								else if (R[0] < -l_2[0])  R[0] += l[0];
								if (R[1] > l_2[1])  R[1] -= l[1];
								else if (R[1] < -l_2[1])  R[1] += l[1];
								if (R[2] > l_2[2])  R[2] -= l[2];
								else if (R[2] < -l_2[2])  R[2] += l[2];
								overlapPear(); 
							}
							vv = linkP[vv];
						}
						NCell.push_back(k);
						usedCell[k]=true;
					}
				}
			}
		}else{
			for (int iz=0; iz < 5 && !inside; iz++){
				kz = WP[1]*s_n[MovedParticle.s[2]][iz];
				for (int iy=0; iy < 5 && !inside; iy++){
					ky = WP[0]*(s_n[MovedParticle.s[1]][iy]+kz);
					for (int ix=0; ix < 5 && !inside; ix++){
						k = s_n[MovedParticle.s[0]][ix] + ky;
						if( usedCell[k] ) continue;
		
						vv = headP[k];

						while(vv != -1 && !inside){
							newvv = vv;
							if(newvv >= N ) newvv -= N; 
							if(v < newvv && !part[newvv].already){
								NPart.push_back(newvv);
								part[newvv].already=true;

								R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
								R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
								R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
								overlapPear(); 
							}
							vv = linkP[vv];
						}
						NCell.push_back(k);
						usedCell[k]=true;
					}
				}
			}
		}
	}

	if(!inside){
		if(MovedParticle.s[0] != MovedParticle.sN[0] || MovedParticle.s[1] != MovedParticle.sN[1] || MovedParticle.s[2] != MovedParticle.sN[2]){
			in_test = false;
					
			if( MovedParticle.sN[0] >= WP[0]-3 || MovedParticle.sN[0] <= 2) in_test = true;
			if( MovedParticle.sN[1] >= WP[1]-3 || MovedParticle.sN[1] <= 2) in_test = true;
			if( MovedParticle.sN[2] >= WP[2]-3 || MovedParticle.sN[2] <= 2) in_test = true;

			if(in_test){
				for (int iz=0; iz < 5 && !inside; iz++){
					kz = WP[1]*s_n[MovedParticle.sN[2]][iz];
					for (int iy=0; iy < 5 && !inside; iy++){
						ky = WP[0]*(s_n[MovedParticle.sN[1]][iy]+kz);
						for (int ix=0; ix < 5 && !inside; ix++){

							k = s_n[MovedParticle.sN[0]][ix] + ky;
							if( usedCell[k] ) continue;
			
							vv = headP[k];
					    
							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= N ) newvv -= N; 
								if(v < newvv && !part[newvv].already){
									NPart.push_back(newvv);
									part[newvv].already=true;

									R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
									R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
									R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
									if (R[0] > l_2[0])  R[0] -= l[0];
									else if (R[0] < -l_2[0])  R[0] += l[0];
									if (R[1] > l_2[1])  R[1] -= l[1];
									else if (R[1] < -l_2[1])  R[1] += l[1];
									if (R[2] > l_2[2])  R[2] -= l[2];
									else if (R[2] < -l_2[2])  R[2] += l[2];
									overlapPear(); 
								}
								vv = linkP[vv];
							}
						}
					}
				}
			}else{
				for (int iz=0; iz < 5 && !inside; iz++){
					kz = WP[1]*s_n[MovedParticle.sN[2]][iz];
					for (int iy=0; iy < 5 && !inside; iy++){
						ky = WP[0]*(s_n[MovedParticle.sN[1]][iy]+kz);
						for (int ix=0; ix < 5 && !inside; ix++){

							k = s_n[MovedParticle.sN[0]][ix] + ky;
							if( usedCell[k] ) continue;
			
							vv = headP[k];
					    
							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= N ) newvv -= N; 
								if(v < newvv && !part[newvv].already){
									NPart.push_back(newvv);
									part[newvv].already=true;

									R[0] = part[newvv].pos[0] - MovedParticle.pos[0];
									R[1] = part[newvv].pos[1] - MovedParticle.pos[1];
									R[2] = part[newvv].pos[2] - MovedParticle.pos[2];
									overlapPear(); 
								}
								vv = linkP[vv];
							}
						}
					}
				}

			}
		}
	}
	for( vv = 0; vv < NPart.size(); vv++) part[NPart[vv]].already=false;
	for( vv = 0; vv < NCell.size(); vv++) usedCell[NCell[vv]]=false;
    }

    if(!inside){ 
	if(Compressing){
		Vbox = Vn;
		rhoN = N/Vbox; 
		rhoV = rhoN*Vsys; 
		vproc = 0.5*(3*vproc-1);

		if(vproc < 0.9 ) vproc=0.9;
		std::cout << "Compession successful! New rho = " << rhoV << " and new vproc=" << vproc << std::endl;
		std::cout << "(lx,ly) = (" << l[0] << "," << l[1] << ")"<< std::endl;

		write(savefile,1);

		if( rhoV - rhoInit < 0.01 ){
			write("Results/Config" + toString(rhoV) + ".dat",0);
			rhoInit = rhoV;
		}
	}else{ 
		std::cout << "Wall-Move successful!" << std::endl;
		zahl1++;
		zahl=0;
		if(zahl1==5){
			vproc = 2*vproc;
			if(vproc > 0.01) vproc = 0.01;
			else std::cout << "New vproc " << vproc << std::endl;
			zahl1 = 0;
		}
	}
	
    }else{
	ln[0] = 1/ln[0];
	ln[1] = 1/ln[1];
	ln[2] = 1/ln[2];
        l[0] *= ln[0];
        l[1] *= ln[1];
        l[2] *= ln[2];

	l_2[0] = l[0]*0.5;
	l_2[1] = l[1]*0.5;
	l_2[2] = l[2]*0.5;

        for( int i=0; i < N; i++){
                part[i].pos[0] *= ln[0];
                part[i].pos[1] *= ln[1];
        	part[i].pos[2] *= ln[2];
        }
	RenewList();

	if(Compressing){
		vproc = 0.5*(1+vproc);
		if(vproc > 0.9999){ 
		    vproc = 0.9999;
		    initCompress = true;
		}
		std::cout << "Compession unsuccessful! New vproc=" << vproc << std::endl;
		std::cout << "(lx,ly,lz) = (" << l[0] << "," << l[1] << ")"<< std::endl;
		if (acceptance > 0.7){
		    if( pos_lambda < maxpos) pos_lambda += 0.05;
		    else pos_lambda = maxpos;
		    if( ori_lambda < maxpos) ori_lambda += 0.05;
		    else ori_lambda = maxpos;
		    std::cout << "and changed pos_lambda=" << pos_lambda << " and ori_lambda=" <<  ori_lambda << std::endl;
		}
		if (acceptance < 0.5){
		    if( pos_lambda < 0.01 ) pos_lambda = 0.01;
		    else pos_lambda -= 0.01;
		    if( ori_lambda < 0.01 ) ori_lambda = 0.01;
		    else ori_lambda -= 0.01;
		    std::cout << "and changed pos_lambda=" << pos_lambda << " and ori_lambda=" <<  ori_lambda << std::endl;
		}
	}else{ 
		zahl++;
		zahl1=0;
		if(zahl==5){
			vproc = vproc*0.5;
			if(vproc < 1e-6) vproc = 1e-6;
			else std::cout << "New vproc " << vproc << std::endl;
			zahl = 0;
		}
	}
    }
}
