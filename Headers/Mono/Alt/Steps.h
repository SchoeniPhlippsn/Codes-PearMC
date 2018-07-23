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
				sx[0] = WP[0]-1;
				sx[1] = sx[0] -1;
				sx[2] = 0;
				in_test = true;
			}else{
				if(dsx > l[0]){
					sx[0] = 0;
					sx[1] = WP[0]-1;
					sx[2] = 1;
					in_test = true;
				}else{
					sx[0] = dsx*wP[0];

					if( sx[0] >= WP[0]-2 || sx[0] <= 1){ 
						in_test = true;

						sx[1] = sx[0]-1;
						if(sx[1] < 0 ) sx[1] = WP[0]-1;

						sx[2] = sx[0]+1;
						if(sx[2]> WP[0]-1) sx[2] = 0;
					}else{
						sx[1] = sx[0]-1;
						sx[2] = sx[0]+1;
					}
				}
			}
		
			

			if(dsy < 0){
				sy[0] = WP[1]-1;
				sy[1] = sy[0] -1;
				sy[2] = 0;
				in_test = true;
			}else{
				if(dsy > l[1]){
					sy[0] = 0;
					sy[1] = WP[1]-1;
					sy[2] = 1;
					in_test = true;
				}else{
					sy[0] = dsy*wP[1];
					if( sy[0] >= WP[1]-2 || sy[0] <= 1){ 
						in_test = true;

						sy[1] = sy[0]-1;
						if(sy[1] < 0 ) sy[1] = WP[1]-1;

						sy[2] = sy[0]+1;
						if(sy[2]> WP[1]-1) sy[2] = 0;
					}else{
						sy[1] = sy[0]-1;
						sy[2] = sy[0]+1;
					}
				}
			}

			if(dsz < 0){
				sz[0] = WP[2]-1;
				sz[1] = sz[0] -1;
				sz[2] = 0;
				in_test = true;
			}else{
				if(dsz > l[2]){
					sz[0] = 0;
					sz[1] = WP[2]-1;
					sz[2] = 1;
					in_test = true;
				}else{
					sz[0] = dsz*wP[2];
					if( sz[0] >= WP[2]-2 || sz[0] <= 1){ 
						in_test = true;

						sz[1] = sz[0]-1;
						if(sz[1] < 0 ) sz[1] = WP[2]-1;

						sz[2] = sz[0]+1;
						if(sz[2]> WP[2]-1) sz[2] = 0;
					}else{
						sz[1] = sz[0]-1;
						sz[2] = sz[0]+1;
					}
				}
			}

			if(in_test){
				for (int iz=0; iz < 3 && !inside; iz++){
					kz = WP[1]*sz[iz];
					for (int iy=0; iy < 3 && !inside; iy++){
						ky = WP[0]*(sy[iy]+kz);
						for (int ix=0; ix < 3 && !inside; ix++){
							k = sx[ix] + ky;
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
				for (int iz=0; iz < 3 && !inside; iz++){
					kz = WP[1]*sz[iz];
					for (int iy=0; iy < 3 && !inside; iy++){
						ky = WP[0]*(sy[iy]+kz);
						for (int ix=0; ix < 3 && !inside; ix++){
							k = sx[ix] + ky;
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
				sxN[0] = WP[0]-1;
				in_test = true;
			}else{
				if(dsxN > l[0]){
					sxN[0] = 0;
					in_test = true;
				}else{
					sxN[0] = dsxN*wP[0];
					if( sxN[0] >= WP[0]-2 || sxN[0] <= 1) in_test = true;
				}
			}

			if(dsyN < 0){
				syN[0] = WP[1]-1;
				in_test = true;
			}else{
				if(dsyN > l[1]){
					syN[0] = 0;
					in_test = true;
				}else{
					syN[0] = dsyN*wP[1];
					if( syN[0] >= WP[1]-2 || syN[0] <= 1) in_test = true;
				}
			}

			if(dszN < 0){
				szN[0] = WP[2]-1;
				in_test = true;
			}else{
				if(dszN > l[2]){
					szN[0] = 0;
					in_test = true;
				}else{
					szN[0] = dszN*wP[2];
					if( szN[0] >= WP[2]-2 || szN[0] <= 1) in_test = true;
				}
			}

			if(sx[0] != sxN[0] || sy[0] != syN[0] || sz[0] != szN[0] ){

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
				
				if(in_test){
					for (int iz=0; iz < 3 && !inside; iz++){
						kz = WP[1]*szN[iz];
						for (int iy=0; iy < 3 && !inside; iy++){
							ky = WP[0]*(syN[iy]+kz);
							for (int ix=0; ix < 3 && !inside; ix++){

								k = sxN[ix] + ky;
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
					for (int iz=0; iz < 3 && !inside; iz++){
						kz = WP[1]*szN[iz];
						for (int iy=0; iy < 3 && !inside; iy++){
							ky = WP[0]*(syN[iy]+kz);
							for (int ix=0; ix < 3 && !inside; ix++){

								k = sxN[ix] + ky;
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
			k = sx[0] + WP[0]*(sy[0]+WP[1]*sz[0]); 

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

			k = sxN[0] + WP[0]*(syN[0]+WP[1]*szN[0]); 

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

		dsx = MovedParticle.pos[0]-MovedParticle.dist_ori[0];
		dsy = MovedParticle.pos[1]-MovedParticle.dist_ori[1];
		dsz = MovedParticle.pos[2]-MovedParticle.dist_ori[2];

		in_test = true;
		if(dsx < 0){
			sx[0] = WP[0]-1;
			sx[1] = sx[0] -1;
			sx[2] = 0;
			in_test = true;
		}else{
			if(dsx > l[0]){
				sx[0] = 0;
				sx[1] = WP[0]-1;
				sx[2] = 1;
				in_test = true;
			}else{
				sx[0] = dsx*wP[0];

				if( sx[0] >= WP[0]-2 || sx[0] <= 1){ 
					in_test = true;

					sx[1] = sx[0]-1;
					if(sx[1] < 0 ) sx[1] = WP[0]-1;

					sx[2] = sx[0]+1;
					if(sx[2]> WP[0]-1) sx[2] = 0;
				}else{
					sx[1] = sx[0]-1;
					sx[2] = sx[0]+1;
				}
			}
		}
	
		

		if(dsy < 0){
			sy[0] = WP[1]-1;
			sy[1] = sy[0] -1;
			sy[2] = 0;
			in_test = true;
		}else{
			if(dsy > l[1]){
				sy[0] = 0;
				sy[1] = WP[1]-1;
				sy[2] = 1;
				in_test = true;
			}else{
				sy[0] = dsy*wP[1];
				if( sy[0] >= WP[1]-2 || sy[0] <= 1){ 
					in_test = true;

					sy[1] = sy[0]-1;
					if(sy[1] < 0 ) sy[1] = WP[1]-1;

					sy[2] = sy[0]+1;
					if(sy[2]> WP[1]-1) sy[2] = 0;
				}else{
					sy[1] = sy[0]-1;
					sy[2] = sy[0]+1;
				}
			}
		}

		if(dsz < 0){
			sz[0] = WP[2]-1;
			sz[1] = sz[0] -1;
			sz[2] = 0;
			in_test = true;
		}else{
			if(dsz > l[2]){
				sz[0] = 0;
				sz[1] = WP[2]-1;
				sz[2] = 1;
				in_test = true;
			}else{
				sz[0] = dsz*wP[2];
				if( sz[0] >= WP[2]-2 || sz[0] <= 1){ 
					in_test = true;

					sz[1] = sz[0]-1;
					if(sz[1] < 0 ) sz[1] = WP[2]-1;

					sz[2] = sz[0]+1;
					if(sz[2]> WP[2]-1) sz[2] = 0;
				}else{
					sz[1] = sz[0]-1;
					sz[2] = sz[0]+1;
				}
			}
		}

		if(in_test){
			for (int ix=0; ix < 3 && !inside; ix++){
				for (int iy=0; iy < 3 && !inside; iy++){
					for (int iz=0; iz < 3 && !inside; iz++){
						k = sx[ix] + WP[0]*(sy[iy]+WP[1]*sz[iz]);
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
			for (int ix=0; ix < 3 && !inside; ix++){
				for (int iy=0; iy < 3 && !inside; iy++){
					for (int iz=0; iz < 3 && !inside; iz++){
						k = sx[ix] + WP[0]*(sy[iy]+WP[1]*sz[iz]);
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
		dsxN = MovedParticle.pos[0]+MovedParticle.dist_ori[0];

		dsyN = MovedParticle.pos[1]+MovedParticle.dist_ori[1];

		dszN = MovedParticle.pos[2]+MovedParticle.dist_ori[2];

		in_test = true;
		if(dsxN < 0){
			sxN[0] = WP[0]-1;
			in_test = true;
		}else{
			if(dsxN > l[0]){
				sxN[0] = 0;
				in_test = true;
			}else{
				sxN[0] = dsxN*wP[0];
				if( sxN[0] >= WP[0]-2 || sxN[0] <= 1) in_test = true;
			}
		}

		if(dsyN < 0){
			syN[0] = WP[1]-1;
			in_test = true;
		}else{
			if(dsyN > l[1]){
				syN[0] = 0;
				in_test = true;
			}else{
				syN[0] = dsyN*wP[1];
				if( syN[0] >= WP[1]-2 || syN[0] <= 1) in_test = true;
			}
		}

		if(dszN < 0){
			szN[0] = WP[2]-1;
			in_test = true;
		}else{
			if(dszN > l[2]){
				szN[0] = 0;
				in_test = true;
			}else{
				szN[0] = dszN*wP[2];
				if( szN[0] >= WP[2]-2 || szN[0] <= 1) in_test = true;
			}
		}

		if(sx[0] != sxN[0] || sy[0] != syN[0] || sz[0] != szN[0]){

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

			if(in_test){
				for (int ix=0; ix < 3 && !inside; ix++){

					for (int iy=0; iy < 3 && !inside; iy++){

						for (int iz=0; iz < 3 && !inside; iz++){

							k = sxN[ix] + WP[0]*(syN[iy]+WP[1]*szN[iz]);
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
				for (int ix=0; ix < 3 && !inside; ix++){

					for (int iy=0; iy < 3 && !inside; iy++){

						for (int iz=0; iz < 3 && !inside; iz++){

							k = sxN[ix] + WP[0]*(syN[iy]+WP[1]*szN[iz]);
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
