void Translation(){

	theta = two_pi*uni(gen);
	phi = 2*uni(gen) - 1;				
	rand_lambda = uni(gen)*pos_lambda;

	cosphi=phi*rand_lambda;
        phi = sqrt(1-phi*phi)*rand_lambda;
        costheta = cos(theta)*phi;
        sintheta = sin(theta)*phi;

	newposx += costheta;
	newposy += sintheta;
	newposz += cosphi;

	if( newposx > lx ) newposx -= lx;
	else if( newposx < 0 ) newposx += lx;

	if( newposy > ly ) newposy -= ly;
	else if( newposy < 0 ) newposy += ly;

	if( newposz > lz ) newposz -= lz;
	else if( newposz < 0 ) newposz += lz;
}

void Rotation(){
	theta = two_pi*uni(gen);
        phi = 2*uni(gen) - 1;
        rand_lambda = ori_lambda;

        cosphi=phi*rand_lambda;
        phi = sqrt(1-phi*phi)*rand_lambda;
        costheta = cos(theta)*phi;
        sintheta = sin(theta)*phi;

	MovedParticle.trans[0][2] +=  costheta;
	MovedParticle.trans[1][2] +=  sintheta;
	MovedParticle.trans[2][2] +=  cosphi;

	sqrtz = MovedParticle.trans[2][2]*MovedParticle.trans[2][2];
	norm = MovedParticle.trans[0][2]*MovedParticle.trans[0][2] + MovedParticle.trans[1][2]*MovedParticle.trans[1][2] + sqrtz;
	norm = 1.0/norm;


	sqrtz = 1-sqrtz*norm;
	if(sqrtz > 1e-4){
		norm = sqrt(norm);
		MovedParticle.trans[0][2] *= norm;
		MovedParticle.trans[1][2] *= norm;
		MovedParticle.trans[2][2] *= norm;
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

	newdist_orix = distN*MovedParticle.trans[0][2];
	newdist_oriy = distN*MovedParticle.trans[1][2];
	newdist_oriz = distN*MovedParticle.trans[2][2];
}


void Check_overlap(){

	NCell.resize(0);
	NPart.resize(0);
	NPart.push_back(randP);
	already[randP]=true;
/*	vv = headP[MovedParticle.cell];

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

//	if( MovedParticle.cellN != MovedParticle.cell && !inside){
	if( !inside){
		vv = headP[MovedParticle.cellN];
	    
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
		NCell.push_back(MovedParticle.cellN);
		usedCell[MovedParticle.cellN]=true;
		NCell.push_back(MovedParticle.cell);
		usedCell[MovedParticle.cell]=true;
	}
*/

	if(!inside){
//		if( MovedParticle.pos[0] > l[0] ) MovedParticle.pos[0] -= l[0];
//		else if( MovedParticle.pos[0] < 0 ) MovedParticle.pos[0] += l[0];
//
//		if( MovedParticle.pos[1] > l[1] ) MovedParticle.pos[1] -= l[1];
//		else if( MovedParticle.pos[1] < 0 ) MovedParticle.pos[1] += l[1];
//
//		if( MovedParticle.pos[2] > l[2] ) MovedParticle.pos[2] -= l[2];
//		else if( MovedParticle.pos[2] < 0 ) MovedParticle.pos[2] += l[2];


		in_test = false;
		if( newsx >= WPx-3 || newsx <= 2) in_test = true;
		else{
 			if( newsy >= WPy-3 || newsy <= 2) in_test = true;
			else if( newsz >= WPz-3 || newsz <= 2) in_test = true;
		}


		if(in_test){
			for (int iz=0; iz < 5 && !inside; iz++){
				kz = WPy*s_n[newsz][iz];
				for (int iy=0; iy < 5 && !inside; iy++){
					ky = WPx*(s_n[newsy][iy]+kz);
					for (int ix=0; ix < 5 && !inside; ix++){
						k = s_n[newsx][ix] + ky;
						//if( usedCell[k] ) continue;
			
						vv = headP[k];

						while(vv != -1 && !inside){
							newvv = vv;
							if(newvv >= N ) newvv -= N; 
							if(!already[newvv]){ 
								NPart.push_back(newvv);
								already[newvv]=true;	
								
								Rx = posx[newvv] - newposx;
								Ry = posy[newvv] - newposy;
								Rz = posz[newvv] - newposz;
								if (Rx > lx_2)  Rx -= lx;
								else if (Rx < -lx_2)  Rx += lx;
								if (Ry > ly_2)  Ry -= ly;
								else if (Ry < -ly_2)  Ry += ly;
								if (Rz > lz_2)  Rz -= lz;
								else if (Rz < -lz_2)  Rz += lz;
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
				kz = WPy*s_n[newsz][iz];
				for (int iy=0; iy < 5 && !inside; iy++){
					ky = WPx*(s_n[newsy][iy]+kz);
					for (int ix=0; ix < 5 && !inside; ix++){
						k = s_n[newsx][ix] + ky;
						if( usedCell[k] ) continue;
			
						vv = headP[k];

						while(vv != -1 && !inside){
							newvv = vv;
							if(newvv >= N ) newvv -= N; 
							if(!already[newvv]){ 
								NPart.push_back(newvv);
								already[newvv]=true;	
								
								Rx = posx[newvv] - newposx;
								Ry = posy[newvv] - newposy;
								Rz = posz[newvv] - newposz;
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

		if(!inside){
			//if(news[0] != newsN[0] || newsy != newsNy || newsz != newsNz ){
				in_test=false;
				if( newsNx >= WPx-3 || newsNx <= 2) in_test = true;
				else{
					if( newsNy >= WPy-3 || newsNy <= 2) in_test = true;
					else if( newsNz >= WPz-3 || newsNz <= 2) in_test = true;
				}

				if(in_test){
					for (int iz=0; iz < 5 && !inside; iz++){
						kz = WPy*s_n[newsNz][iz];
						for (int iy=0; iy < 5 && !inside; iy++){
							ky = WPx*(s_n[newsNy][iy]+kz);
							for (int ix=0; ix < 5 && !inside; ix++){

								k = s_n[newsNx][ix] + ky;
								if( usedCell[k] ) continue;
				
								vv = headP[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= N ) newvv -= N; 
									if(!already[newvv]){ 
										NPart.push_back(newvv);
										already[newvv]=true;	

										Rx = posx[newvv] - newposx;
										Ry = posy[newvv] - newposy;
										Rz = posz[newvv] - newposz;
										if (Rx > lx_2)  Rx -= lx;
										else if (Rx < -lx_2)  Rx += lx;
										if (Ry > ly_2)  Ry -= ly;
										else if (Ry < -ly_2)  Ry += ly;
										if (Rz > lz_2)  Rz -= lz;
										else if (Rz < -lz_2)  Rz += lz;
										overlapPear(); 
									}
									vv = linkP[vv];
								}
							}
						}
					}
				}else{
					for (int iz=0; iz < 5 && !inside; iz++){
						kz = WPy*s_n[newsNz][iz];
						for (int iy=0; iy < 5 && !inside; iy++){
							ky = WPx*(s_n[newsNy][iy]+kz);
							for (int ix=0; ix < 5 && !inside; ix++){

								k = s_n[newsNx][ix] + ky;
								if( usedCell[k] ) continue;
				
								vv = headP[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= N ) newvv -= N; 
									if(!already[newvv]){ 
										NPart.push_back(newvv);
										already[newvv]=true;	

										Rx = posx[newvv] - newposx;
										Ry = posy[newvv] - newposy;
										Rz = posz[newvv] - newposz;
										overlapPear(); 
									}
									vv = linkP[vv];
								}
							}
						}
					}
				}
		//	}
		}
	}

	for( vv = 0; vv < NPart.size(); vv++) already[NPart[vv]]=false;
	for( vv = 0; vv < NCell.size(); vv++) usedCell[NCell[vv]]=false;
}


void Replace_part(){
	acc++;

	dsx = newposx-newdist_orix;
	dsy = newposy-newdist_oriy;
	dsz = newposz-newdist_oriz;

	if(dsx < 0) newsx = WPx-1;
	else{
		if(dsx > lx) newsx = 0;
		else newsx = dsx*wPx;
	}
	

	if(dsy < 0) newsy = WPy-1;
	else{
		if(dsy > ly) newsy = 0;
		else newsy = dsy*wPy;
	}

	if(dsz < 0) newsz = WPz-1;
	else{
		if(dsz > lz) newsz = 0;
		else newsz = dsz*wPz;
	}


	k = newsx + WPx*(newsy+WPy*newsz); 

	if(newcell != k){
		if(headP[newcell] != randP ){
			int vvv = headP[newcell];
			vv = linkP[vvv];
			while(vv!=randP){ 
				vvv = vv;
				vv = linkP[vvv];
			}
	
			linkP[vvv] = linkP[randP];
		}else{
			headP[newcell] = linkP[randP];
		}

		linkP[randP] = headP[k];
		headP[k] = randP;
		newcell = k;
	}

	dsxN = newposx+newdist_orix;

	dsyN = newposy+newdist_oriy;

	dszN = newposz+newdist_oriz;

	if(dsxN < 0) newsNx = WPx-1;
	else{
		if(dsxN > lx) newsNx = 0;
		else newsNx = dsxN*wPx;
	}

	if(dsyN < 0) newsNy = WPy-1;
	else{
		if(dsyN > ly) newsNy = 0;
		else newsNy = dsyN*wPy;
	}

	if(dszN < 0) newsNz = WPz-1;
	else{
		if(dszN > lz) newsNz = 0;
		else newsNz = dszN*wPz;
	}

	k = newsNx + WPx*(newsNy+WPy*newsNz); 

	if(newcellN != k){
		if(headP[newcellN] != randP+N ){
			int vvv = headP[newcellN];
			vv = linkP[vvv];
			while(vv!=randP+N){ 
				vvv = vv;
				vv = linkP[vvv];
			}
	
			linkP[vvv] = linkP[randP+N];
		}else{
			headP[newcellN] = linkP[randP+N];
		}

		linkP[randP+N] = headP[k];
		headP[k] = randP+N;
		newcellN = k;
	}
	posx[randP] = newposx;
	posy[randP] = newposy;
	posz[randP] = newposz;

	pos_msdx[randP] = newpos_msdx;
	pos_msdy[randP] = newpos_msdy;
	pos_msdz[randP] = newpos_msdz;

	dist_orix[randP] = newdist_orix;
	dist_oriy[randP] = newdist_oriy;
	dist_oriz[randP] = newdist_oriz;

	sx[randP] = newsx;
	sy[randP] = newsy;
	sz[randP] = newsz;

	sNx[randP] = newsNx;
	sNy[randP] = newsNy;
	sNz[randP] = newsNz;

	cell[randP]=newcell;
	cellN[randP]=newcellN;

	part[randP]=MovedParticle;
}



void Move_step(){
	for(int v = 0;  v < N; v++){
        	randP = uni(gen)*N;
		
        	if (randP==N){
            		v--;
            		continue;
		}
        	MovedParticle = part[randP];         
		
		newposx = posx[randP];
		newposy = posy[randP];
		newposz = posz[randP];

		newpos_msdx = pos_msdx[randP];
		newpos_msdy = pos_msdy[randP];
		newpos_msdz = pos_msdz[randP];

		newdist_orix = dist_orix[randP];
		newdist_oriy = dist_oriy[randP];
		newdist_oriz = dist_oriz[randP];

		newcell = cell[randP];
		newcellN = cellN[randP];

		newsx = sx[randP];
		newsy = sy[randP];
		newsz = sz[randP];

		newsNx = sNx[randP];
		newsNy = sNy[randP];
		newsNz = sNz[randP];
		
        	inside = false;


		if( 2*uni(gen) < 1 ){
			Translation();

			newpos_msdx += costheta;
			newpos_msdy += sintheta;
			newpos_msdz += cosphi;

		}else Rotation();

		pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

		Check_overlap();

		if(!inside) Replace_part();
	}
}

void Trans_step(){
	for(int v = 0;  v < N; v++){
        	randP = uni(gen)*N;
		
        	if (randP==N){
            		v--;
            		continue;
		}
        	MovedParticle = part[randP];         
		
		newposx = posx[randP];
		newposy = posy[randP];
		newposz = posz[randP];

		newpos_msdx = pos_msdx[randP];
		newpos_msdy = pos_msdy[randP];
		newpos_msdz = pos_msdz[randP];

		newdist_orix = dist_orix[randP];
		newdist_oriy = dist_oriy[randP];
		newdist_oriz = dist_oriz[randP];

		newcell = cell[randP];
		newcellN = cellN[randP];

		newsx = sx[randP];
		newsy = sy[randP];
		newsz = sz[randP];

		newsNx = sNx[randP];
		newsNy = sNy[randP];
		newsNz = sNz[randP];

        	inside = false;

		Translation();

		pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

		Check_overlap();

		if(!inside) Replace_part();
	}
}

void Rot_step(){
	for(int v = 0;  v < N; v++){
        	randP = uni(gen)*N;
		
        	if (randP==N){
            		v--;
            		continue;
		}
        	MovedParticle = part[randP];         
		
		newposx = posx[randP];
		newposy = posy[randP];
		newposz = posz[randP];

		newpos_msdx = pos_msdx[randP];
		newpos_msdy = pos_msdy[randP];
		newpos_msdz = pos_msdz[randP];

		newdist_orix = dist_orix[randP];
		newdist_oriy = dist_oriy[randP];
		newdist_oriz = dist_oriz[randP];

		newcell = cell[randP];
		newcellN = cellN[randP];

		newsx = sx[randP];
		newsy = sy[randP];
		newsz = sz[randP];

		newsNx = sNx[randP];
		newsNy = sNy[randP];
		newsNz = sNz[randP];

        	inside = false;

		Rotation();

		pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

		Check_overlap();

		if(!inside) Replace_part();
	}
}



void Compression_step(){	

    if(Compressing){
	Vn = log(Vbox)*vproc;
	Vn = exp(Vn);
	powVn = pow(Vn,1.0/3.0);

    	lxn = powVn/lx;
    	lyn = powVn/ly;
    	lzn = powVn/lz;

    }else{
	Vn = (1+uni(gen)*vproc);
	powVn = sqrt(Vn);

	Case = 2*uni(gen);
	if(Case==0) powVn = 1/powVn;
	else Vn = 1/Vn;

	Case = 3*uni(gen);
	switch(Case){
		case 0:
		lxn=Vn;
		lyn=powVn;
		lzn=powVn;
		break;
		case 1:
		lxn=powVn;
		lyn=Vn;
		lzn=powVn;
		break;
		case 2:
		lxn=powVn;
		lyn=powVn;
		lzn=Vn;
		break;
	}
   }

    lx *= lxn;
    ly *= lyn;
    lz *= lzn;

    lx_2 = lx*0.5;
    ly_2 = ly*0.5;
    lz_2 = lz*0.5;


    for( int i=0; i < N; i++){
            posx[i] *= lxn;
            posy[i] *= lzn;
            posz[i] *= lyn;
    }
    RenewList();
    inside = false;

    for( int v=0; v < N-1 && !inside; v++){
	MovedParticle = part[v];         
	
	newposx = posx[v];
	newposy = posy[v];
	newposz = posz[v];

	newpos_msdx = pos_msdx[v];
	newpos_msdy = pos_msdy[v];
	newpos_msdz = pos_msdz[v];

	newdist_orix = dist_orix[v];
	newdist_oriy = dist_oriy[v];
	newdist_oriz = dist_oriz[v];

	newcell = cell[v];
	newcellN = cellN[v];

	newsx = sx[v];
	newsy = sy[v];
	newsz = sz[v];

	newsNx = sNx[v];
	newsNy = sNy[v];
	newsNz = sNz[v];

	pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

	NCell.resize(0);
	NPart.resize(0);

	k = newcell;

	vv = headP[k];

	while(vv != -1 && !inside){
		newvv = vv;
		if(newvv >= N ) newvv -= N; 
		if(v < newvv && !already[newvv]){ 
			NPart.push_back(newvv);
			already[newvv]=true;

			Rx = posx[newvv] - newposx;
			Ry = posy[newvv] - newposy;
			Rz = posz[newvv] - newposz;
			if (Rx > lx_2)  Rx -= lx;
			else if (Rx < -lx_2)  Rx += lx;
			if (Ry > ly_2)  Ry -= ly;
			else if (Ry < -ly_2)  Ry += ly;
			if (Rz > lz_2)  Rz -= lz;
			else if (Rz < -lz_2)  Rz += lz;
			overlapPear(); 
		}

		vv = linkP[vv];
	}
	NCell.push_back(k);
	usedCell[k]=true;

	if( newcellN != newcell && !inside){
		k = newcellN;

		vv = headP[k];
	    
		while(vv != -1 && !inside){
			newvv = vv;
			if(newvv >= N ) newvv -= N; 
			if(v < newvv && !already[newvv]){ 
				NPart.push_back(newvv);
				already[newvv]=true;

				Rx = posx[newvv] - newposx;
				Ry = posy[newvv] - newposy;
				Rz = posz[newvv] - newposz;
				if (Rx > lx_2)  Rx -= lx;
				else if (Rx < -lx_2)  Rx += lx;
				if (Ry > ly_2)  Ry -= ly;
				else if (Ry < -ly_2)  Ry += ly;
				if (Rz > lz_2)  Rz -= lz;
				else if (Rz < -lz_2)  Rz += lz;
				overlapPear(); 
			}
			vv = linkP[vv];
		}
		NCell.push_back(k);
		usedCell[k]=true;
	}

	if(!inside){

		in_test = false;
		if( newsx >= WPx-3 || newsx <= 2) in_test = true;
		if( newsy >= WPy-3 || newsy <= 2) in_test = true;
		if( newsz >= WPz-3 || newsz <= 2) in_test = true;

		if(in_test){
			for (int iz=0; iz < 5 && !inside; iz++){
				kz = WPy*s_n[newsz][iz];
				for (int iy=0; iy < 5 && !inside; iy++){
					ky = WPx*(s_n[newsy][iy]+kz);
					for (int ix=0; ix < 5 && !inside; ix++){
						k = s_n[newsx][ix] + ky;
						if( usedCell[k] ) continue;
		
						vv = headP[k];

						while(vv != -1 && !inside){
							newvv = vv;
							if(newvv >= N ) newvv -= N; 
							if(v < newvv && !already[newvv]){
								NPart.push_back(newvv);
								already[newvv]=true;

								Rx = posx[newvv] - newposx;
								Ry = posy[newvv] - newposy;
								Rz = posz[newvv] - newposz;
								if (Rx > lx_2)  Rx -= lx;
								else if (Rx < -lx_2)  Rx += lx;
								if (Ry > ly_2)  Ry -= ly;
								else if (Ry < -ly_2)  Ry += ly;
								if (Rz > lz_2)  Rz -= lz;
								else if (Rz < -lz_2)  Rz += lz;
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
				kz = WPy*s_n[newsz][iz];
				for (int iy=0; iy < 5 && !inside; iy++){
					ky = WPx*(s_n[newsy][iy]+kz);
					for (int ix=0; ix < 5 && !inside; ix++){
						k = s_n[newsx][ix] + ky;
						if( usedCell[k] ) continue;
		
						vv = headP[k];

						while(vv != -1 && !inside){
							newvv = vv;
							if(newvv >= N ) newvv -= N; 
							if(v < newvv && !already[newvv]){
								NPart.push_back(newvv);
								already[newvv]=true;

								Rx = posx[newvv] - newposx;
								Ry = posy[newvv] - newposy;
								Rz = posz[newvv] - newposz;
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
		if(newsx != newsNx || newsy != newsNy || newsz != newsNz){
			in_test = false;
					
			if( newsNx >= WPx-3 || newsNx <= 2) in_test = true;
			if( newsNy >= WPy-3 || newsNy <= 2) in_test = true;
			if( newsNz >= WPz-3 || newsNz <= 2) in_test = true;

			if(in_test){
				for (int iz=0; iz < 5 && !inside; iz++){
					kz = WPy*s_n[newsNz][iz];
					for (int iy=0; iy < 5 && !inside; iy++){
						ky = WPx*(s_n[newsNy][iy]+kz);
						for (int ix=0; ix < 5 && !inside; ix++){

							k = s_n[newsNx][ix] + ky;
							if( usedCell[k] ) continue;
			
							vv = headP[k];
					    
							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= N ) newvv -= N; 
								if(v < newvv && !already[newvv]){
									NPart.push_back(newvv);
									already[newvv]=true;

									Rx = posx[newvv] - newposx;
									Ry = posy[newvv] - newposy;
									Rz = posz[newvv] - newposz;
									if (Rx > lx_2)  Rx -= lx;
									else if (Rx < -lx_2)  Rx += lx;
									if (Ry > ly_2)  Ry -= ly;
									else if (Ry < -ly_2)  Ry += ly;
									if (Rz > lz_2)  Rz -= lz;
									else if (Rz < -lz_2)  Rz += lz;
									overlapPear(); 
								}
								vv = linkP[vv];
							}
						}
					}
				}
			}else{
				for (int iz=0; iz < 5 && !inside; iz++){
					kz = WPy*s_n[newsNz][iz];
					for (int iy=0; iy < 5 && !inside; iy++){
						ky = WPx*(s_n[newsNy][iy]+kz);
						for (int ix=0; ix < 5 && !inside; ix++){

							k = s_n[newsNx][ix] + ky;
							if( usedCell[k] ) continue;
			
							vv = headP[k];
					    
							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= N ) newvv -= N; 
								if(v < newvv && !already[newvv]){
									NPart.push_back(newvv);
									already[newvv]=true;

									Rx = posx[newvv] - newposx;
									Ry = posy[newvv] - newposy;
									Rz = posz[newvv] - newposz;
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
	for( vv = 0; vv < NPart.size(); vv++) already[NPart[vv]]=false;
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
		std::cout << "(lx,ly,lz) = (" << lx << "," << ly << "," << lz << ")"<< std::endl;

		writeSave();

		if( rhoV - rhoInit < 0.01 ){
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
	lxn = 1/lxn;
	lyn = 1/lyn;
	lzn = 1/lzn;
        lx *= lxn;
        ly *= lyn;
        lz *= lzn;

	lx_2 = lx*0.5;
	ly_2 = ly*0.5;
	lz_2 = lz*0.5;

        for( int i=0; i < N; i++){
                posx[i] *= lxn;
                posy[i] *= lyn;
        	posz[i] *= lzn;
        }
	RenewList();

	if(Compressing){
		vproc = 0.5*(1+vproc);
		if(vproc > 0.9999){ 
		    vproc = 0.9999;
		    initCompress = true;
		}
		std::cout << "Compession unsuccessful! New vproc=" << vproc << std::endl;
		std::cout << "(lx,ly,lz) = (" << lx << "," << ly << "," << lz <<  ")"<< std::endl;
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
