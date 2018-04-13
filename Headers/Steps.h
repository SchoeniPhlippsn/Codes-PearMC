void Move_step(){
	for(int v = 0;  v < N; v++){

      
        	randP = uni(gen)*N;
        	if (randP==N){
            		v--;
            		continue;
		}
       
        	MovedParticle = Config.part[randP];         

        	inside = false;

        	if( randP < Config.Nc ){
            		if( 2*uni(gen) < 1 ){
                		MovedParticle.pos[0] = MovedParticle.pos[0] + (uni(gen)-0.5)*pos_lambda;
                		MovedParticle.pos[1] = MovedParticle.pos[1] + (uni(gen)-0.5)*pos_lambda;
                		MovedParticle.pos[2] = MovedParticle.pos[2] + (uni(gen)-0.5)*pos_lambda;
            		}else{
                		phi = (uni(gen)-0.5)*ori_lambda;
                		cosphi = cos(phi);
                		sinphi = sin(phi);

                		 Case = 3*uni(gen);
                		switch(Case){
                			case 0:
                    			MovedParticle.ori[0] = cosphi*MovedParticle.ori[0] + sinphi*MovedParticle.ori[1];
                    			MovedParticle.ori[1] = -sinphi*MovedParticle.ori[0] + cosphi*MovedParticle.ori[1];
                    			break;

                			case 1:
                    			MovedParticle.ori[2] = cosphi*MovedParticle.ori[2] + sinphi*MovedParticle.ori[0];
                    			MovedParticle.ori[0] = -sinphi*MovedParticle.ori[2] + cosphi*MovedParticle.ori[0];
                    			break;

                			case 2:
                    			MovedParticle.ori[1] = cosphi*MovedParticle.ori[1] + sinphi*MovedParticle.ori[2];
                    			MovedParticle.ori[2] = -sinphi*MovedParticle.ori[1] + cosphi*MovedParticle.ori[2];
                    			break;
                		}

                		double norm = 1/sqrt(MovedParticle.ori[0]*MovedParticle.ori[0] + MovedParticle.ori[1]*MovedParticle.ori[1] + MovedParticle.ori[2]*MovedParticle.ori[2]);
                		MovedParticle.ori[0] *= norm;
                		MovedParticle.ori[1] *= norm;
                		MovedParticle.ori[2] *= norm;

				sqrtz = 1-MovedParticle.ori[2]*MovedParticle.ori[2];
				if(sqrtz > 1e-4){
					sqrtz = sqrt(sqrtz);

					MovedParticle.trans[2][0] = -sqrtz;

					MovedParticle.trans[0][2] = MovedParticle.ori[0];
					MovedParticle.trans[1][2] = MovedParticle.ori[1];
					MovedParticle.trans[2][2] = MovedParticle.ori[2];

					sqrtz = 1/sqrtz;
					MovedParticle.trans[0][1] = -MovedParticle.ori[1]*sqrtz;
					MovedParticle.trans[1][1] = MovedParticle.ori[0]*sqrtz;

					MovedParticle.trans[0][0] = MovedParticle.ori[2]*MovedParticle.trans[1][1];
					MovedParticle.trans[1][0] = -MovedParticle.ori[2]*MovedParticle.trans[0][1];

				}else{
					MovedParticle.trans[0][0] = 1;
					MovedParticle.trans[0][1] = 0;
					MovedParticle.trans[0][2] = 0;
					
					MovedParticle.trans[1][0] = 0;
					MovedParticle.trans[1][2] = 0;


					MovedParticle.trans[2][0] = 0;

					if(MovedParticle.ori[2]>0){
						MovedParticle.trans[1][1] = 1;
						MovedParticle.trans[2][2] = 1;
					}else{
						MovedParticle.trans[1][1] = -1;
						MovedParticle.trans[2][2] = -1;
					}
				}
            		}
			pear_mesh.UpdateTrans(id[0], MovedParticle.trans);
			if(Config.Ns == 0 ){

				NCell.resize(0);
				NPart.resize(0);
				first=true;
				vorzeichen = 1;

				int k = MovedParticle.cell;

		    		vv = Config.head[k];
		    
		    		while(vv != -1 && !inside){
					newvv = vv;
					if(newvv >= Config.Nc ) newvv -= Config.Nc; 
					if(randP != newvv && !Config.part[newvv].already) overlapPear(); 

					vv = Config.link[vv];
		    		}
				NCell.push_back(k);
				Config.usedCell[k]=true;

				int kN = MovedParticle.cell;
				if( MovedParticle.cellN != MovedParticle.cell && !inside){
					k = MovedParticle.cellN;

					vv = Config.head[k];
			    
					while(vv != -1 && !inside){
						newvv = vv;
						if(newvv >= Config.Nc ) newvv -= Config.Nc; 
						if(randP != newvv && !Config.part[newvv].already) overlapPear(); 

						vv = Config.link[vv];
					}
					NCell.push_back(k);
					Config.usedCell[k]=true;
				}


				if(!inside){

				dsx = MovedParticle.pos[0]-MovedParticle.dist*MovedParticle.ori[0];
				dsxN = MovedParticle.pos[0]+MovedParticle.dist*MovedParticle.ori[0];

				dsy = MovedParticle.pos[1]-MovedParticle.dist*MovedParticle.ori[1];
				dsyN = MovedParticle.pos[1]+MovedParticle.dist*MovedParticle.ori[1];

				dsz = MovedParticle.pos[2]-MovedParticle.dist*MovedParticle.ori[2];
				dszN = MovedParticle.pos[2]+MovedParticle.dist*MovedParticle.ori[2];

					if(dsx < 0){
						sx[0] = Config.W[0]-1;
						sx[1] = sx[0] -1;
						sx[2] = 0;
					}else{
						sx[0] = dsx*Config.w[0];
						if(sx[0]>=Config.W[0]) sx[0]-=Config.W[0];

						sx[1] = sx[0]-1;
						if(sx[1] < 0 ) sx[1]+=Config.W[0];

						sx[2] = sx[0]+1;
						if(sx[2]>=Config.W[0]) sx[2]-=Config.W[0];
					}

					if(dsy < 0){
						sy[0] = Config.W[1]-1;
						sy[1] = sy[0] -1;
						sy[2] = 0;
					}else{
						sy[0] = dsy*Config.w[1];
						if(sy[0]>=Config.W[1]) sy[0]-=Config.W[1];

						sy[1] = sy[0]-1;
						if(sy[1] < 0 ) sy[1]+=Config.W[1];

						sy[2] = sy[0]+1;
						if(sy[2]>=Config.W[1]) sy[2]-=Config.W[1];
					}

					if(dsz < 0){
						sz[0] = Config.W[2]-1;
						sz[1] = sz[0] -1;
						sz[2] = 0;
					}else{
						sz[0] = dsz*Config.w[2];
						if(sz[0]>=Config.W[2]) sz[0]-=Config.W[2];

						sz[1] = sz[0]-1;
						if(sz[1] < 0 ) sz[1]+=Config.W[2];

						sz[2] = sz[0]+1;
						if(sz[2]>=Config.W[2]) sz[2]-=Config.W[2];
					}

					for (int ix=0; ix < 3 && !inside; ix++){

						for (int iy=0; iy < 3 && !inside; iy++){

							for (int iz=0; iz <= 3 && !inside; iz++){
								int k = sx[ix] + Config.W[0]*(sy[iy] + Config.W[1]*sz[iz]);
								if( Config.usedCell[k] ) continue;
				
								vv = Config.head[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.Nc ) newvv -= Config.Nc; 
									if(randP != newvv && !Config.part[newvv].already) overlapPear(); 
									vv = Config.link[vv];
								}
								NCell.push_back(k);
								Config.usedCell[k]=true;
							}
						}
					}
				}

			
				if(!inside){
					first=false;

					dsxN = MovedParticle.pos[0]+MovedParticle.dist*MovedParticle.ori[0];

					dsyN = MovedParticle.pos[1]+MovedParticle.dist*MovedParticle.ori[1];

					dszN = MovedParticle.pos[2]+MovedParticle.dist*MovedParticle.ori[2];

					if(dsxN < 0){
						sxN[0] = Config.W[0]-1;
						sxN[1] = sxN[0] -1;
						sxN[2] = 0;
					}else{
						sxN[0] = dsxN*Config.w[0];
						if(sxN[0]>=Config.W[0]) sxN[0]-=Config.W[0];

						sxN[1] = sxN[0]-1;
						if(sxN[1] < 0 ) sxN[1]+=Config.W[0];

						sxN[2] = sxN[0]+1;
						if(sxN[2]>=Config.W[0]) sxN[2]-=Config.W[0];
					}

					if(dsyN < 0){
						syN[0] = Config.W[1]-1;
						syN[1] = syN[0] -1;
						syN[2] = 0;
					}else{
						syN[0] = dsyN*Config.w[1];
						if(syN[0]>=Config.W[1]) syN[0]-=Config.W[1];

						syN[1] = syN[0]-1;
						if(syN[1] < 0 ) syN[1]+=Config.W[1];

						syN[2] = syN[0]+1;
						if(syN[2]>=Config.W[1]) syN[2]-=Config.W[1];
					}

					if(dszN < 0){
						szN[0] = Config.W[2]-1;
						szN[1] = szN[0] -1;
						szN[2] = 0;
					}else{
						szN[0] = dszN*Config.w[2];
						if(szN[0]>=Config.W[2]) szN[0]-=Config.W[2];

						szN[1] = szN[0]-1;
						if(szN[1] < 0 ) szN[1]+=Config.W[2];

						szN[2] = szN[0]+1;
						if(szN[2]>=Config.W[2]) szN[2]-=Config.W[2];
					}

					if(sx[1] != sxN[1] || sy[1] != syN[1] || sz[1] != szN[1] ){
						
						vorzeichen=-1;
						for (int ix=0; ix < 3 && !inside; ix++){

							for (int iy=0; iy < 3 && !inside; iy++){

								for (int iz=0; iz <= 3 && !inside; iz++){
									int k = sxN[ix] + Config.W[0]*(syN[iy] + Config.W[1]*szN[iz]);
									if( Config.usedCell[k] ) continue;
					
									vv = Config.head[k];
							    
									while(vv != -1 && !inside){
										newvv = vv;
										if(newvv >= Config.Nc ) newvv -= Config.Nc; 
										if(randP != newvv && !Config.part[newvv].already) overlapPear(); 
										vv = Config.link[vv];
									}
								}
							}
						}
					}
				}

				for( vv = 0; vv < NPart.size(); vv++) Config.part[NPart[vv]].already=false;
				for( vv = 0; vv < NCell.size(); vv++) Config.usedCell[NCell[vv]]=false;
	        	}else{
				for( vv = 0; vv < Config.Nc && !inside; vv++ ){
					if(randP==vv) continue;
					if(overlapP(MovedParticle,Config.part[vv], Config.l)) inside=true; 
				}

				for( vv = Config.Nc; vv < Config.part.size() && !inside; vv++ ){
					if(overlapPSPH(MovedParticle,Config.part[vv], Config.l)) inside=true;
				}
	        	}

        	}else{
			MovedParticle.pos[0] = MovedParticle.pos[0] + (uni(gen)-0.5)*pos_lambda;
			MovedParticle.pos[1] = MovedParticle.pos[1] + (uni(gen)-0.5)*pos_lambda;
			MovedParticle.pos[2] = MovedParticle.pos[2] + (uni(gen)-0.5)*pos_lambda;

            		for( vv = 0; vv < Config.Nc && !inside; vv++ ){
                		if(overlapPSPH(Config.part[vv],MovedParticle, Config.l)) inside=true; 
            		}

	    		if(!inside){

				int k = MovedParticle.cell;

		    		vv = Config.head[k];
		    
		    		while(vv != -1 && !inside){
					if(randP != vv ) overlapSphere();

					vv = Config.link[vv];
		    		}
			

		    		sx[0] = MovedParticle.pos[0]*Config.w[0];
		    		if(sx[0]>=Config.W[0]-1){ 
					sx[0] = Config.W[0]-1;
					sx[1] = Config.W[0]-2;
					sx[2] = 0;
		    		}else{
					sx[2]=sx[0]+1;
					if(sx[0]==0) sx[1] = Config.W[0]-1;
					else sx[1] = sx[0] -1;
	    	   		}

		    		sy[0] = MovedParticle.pos[1]*Config.w[1];
		    		if(sy[0]>=Config.W[1]-1){ 
					sy[0] = Config.W[1]-1;
					sy[1] = Config.W[1]-2;
					sy[2] = 0;
		    		}else{
					sy[2]=sy[0]+1;
					if(sy[0]==0) sy[1] = Config.W[1]-1;
					else sy[1] = sy[0] -1;
		    		}

		    		sz[0] = MovedParticle.pos[2]*Config.w[2];
		    		if(sz[0]>=Config.W[2]-1){ 
					sz[0] = Config.W[2]-1;
					sz[1] = Config.W[2]-2;
					sz[2] = 0;
		    		}else{
					sz[2]=sz[0]+1;
					if(sz[0]==0) sz[1] = Config.W[2]-1;
					else sz[1] = sz[0] -1;
		    		}

		    		for (int ix=0; ix < 3 && !inside; ix++){

					for (int iy=0; iy < 3 && !inside; iy++){

						for (int iz=0; iz <= 3 && !inside; iz++){
					    		k = sx[ix] + Config.W[0]*(sy[iy] + Config.W[1]*sz[iz]);
					    		if( k == MovedParticle.cell ) continue;
		    	
					    		vv = Config.head[k];
					    
					    		while(vv != -1 && !inside){
								overlapSphere(); 
								vv = Config.link[vv];
					    		}
						}
					}
		    		}
			}
		}


		if(!inside){
	            	acc++;

                	MovedParticle.pos_msd[0] += (MovedParticle.pos[0]-Config.part[randP].pos[0]);
                	MovedParticle.pos_msd[1] += (MovedParticle.pos[1]-Config.part[randP].pos[1]);
                	MovedParticle.pos_msd[2] += (MovedParticle.pos[2]-Config.part[randP].pos[2]);

	            	if( MovedParticle.pos[0] > Config.l[0] ) MovedParticle.pos[0] -= Config.l[0];
        	    	else if( MovedParticle.pos[0] < 0 ) MovedParticle.pos[0] += Config.l[0];
		
	            	if( MovedParticle.pos[1] > Config.l[1] ) MovedParticle.pos[1] -= Config.l[1];
        	    	else if( MovedParticle.pos[1] < 0 ) MovedParticle.pos[1] += Config.l[1];

	            	if( MovedParticle.pos[2] > Config.l[2] ) MovedParticle.pos[2] -= Config.l[2];
        	    	else if( MovedParticle.pos[2] < 0 ) MovedParticle.pos[2] += Config.l[2];
	
           		if (Config.Ns == 0){
				int k = sx[0] + Config.W[0]*(sy[0] + Config.W[1]*sz[0]); 

				if(MovedParticle.cell != k){
					if(Config.head[MovedParticle.cell] != randP ){
						int vvv = Config.head[MovedParticle.cell];
						vv = Config.link[vvv];
						while(vv!=randP){ 
							vvv = vv;
							vv = Config.link[vvv];
						}
				
						Config.link[vvv] = Config.link[randP];
					}else{
						Config.head[MovedParticle.cell] = Config.link[randP];
					}
		    
					Config.link[randP] = Config.head[k];
					Config.head[k] = randP;
					MovedParticle.cell = k;
				}

				k = sxN[0] + Config.W[0]*(syN[0] + Config.W[1]*szN[0]); 

				if(MovedParticle.cellN != k){
					if(Config.head[MovedParticle.cellN] != randP+Config.Nc ){
						int vvv = Config.head[MovedParticle.cellN];
						vv = Config.link[vvv];
						while(vv!=randP+Config.Nc){ 
							vvv = vv;
							vv = Config.link[vvv];
						}
				
						Config.link[vvv] = Config.link[randP+Config.Nc];
					}else{
						Config.head[MovedParticle.cellN] = Config.link[randP+Config.Nc];
					}
		    
					Config.link[randP+Config.Nc] = Config.head[k];
					Config.head[k] = randP+Config.Nc;
					MovedParticle.cellN = k;
				}


			}else{ 
				if( randP >= Config.Nc){
					int k = sx[0] + Config.W[0]*(sy[0] + Config.W[1]*sz[0]); 

					if(MovedParticle.cell != k){
						if(Config.head[MovedParticle.cell] != randP ){
							int vvv = Config.head[MovedParticle.cell];
							vv = Config.link[vvv];
							while(vv!=randP){ 
								vvv = vv;
								vv = Config.link[vvv];
							}
					
							Config.link[vvv] = Config.link[randP];
						}else{
							Config.head[MovedParticle.cell] = Config.link[randP];
						}
			    
						Config.link[randP] = Config.head[k];
						Config.head[k] = randP;
						MovedParticle.cell = k;
					}
				}
			}
			Config.part[randP] = MovedParticle;
        	}
	}
}



void Compression_step(){	

    if(Compressing){
	Vn = log(Config.Vbox)*vproc;
	Vn = exp(Vn);
	powVn = pow(Vn,1.0/3.0);

    	ln[0] = powVn/Config.l[0];
    	ln[1] = powVn/Config.l[1];
    	ln[2] = powVn/Config.l[2];

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

    Config.l[0] *= ln[0];
    Config.l[1] *= ln[1];
    Config.l[2] *= ln[2];

    Config.l_2[0] = Config.l[0]*0.5;
    Config.l_2[1] = Config.l[1]*0.5;
    Config.l_2[2] = Config.l[2]*0.5;


    for( int i=0; i < Config.part.size(); i++){
            Config.part[i].pos[0] *= ln[0];
            Config.part[i].pos[1] *= ln[1];
            Config.part[i].pos[2] *= ln[2];
    }
    RenewList();
    inside = false;

    for( int v=0; v < Config.part.size()-1 && !inside; v++){
        if(v < Config.Nc){
		if(Config.Ns == 0 ){
			MovedParticle=Config.part[v];
			pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

			NCell.resize(0);
			NPart.resize(0);
			first=true;


			dsx = MovedParticle.pos[0]-MovedParticle.dist*MovedParticle.ori[0];
			dsxN = MovedParticle.pos[0]+MovedParticle.dist*MovedParticle.ori[0];

			dsy = MovedParticle.pos[1]-MovedParticle.dist*MovedParticle.ori[1];
			dsyN = MovedParticle.pos[1]+MovedParticle.dist*MovedParticle.ori[1];

			dsz = MovedParticle.pos[2]-MovedParticle.dist*MovedParticle.ori[2];
			dszN = MovedParticle.pos[2]+MovedParticle.dist*MovedParticle.ori[2];

			int k = MovedParticle.cell;

			vv = Config.head[k];
	    
			while(vv != -1 && !inside){
				newvv = vv;
				if(newvv >= Config.Nc ) newvv -= Config.Nc; 
				if(v < newvv && !Config.part[newvv].already) overlapPear(); 

				vv = Config.link[vv];
			}
			NCell.push_back(k);
			Config.usedCell[k]=true;
			
			if( MovedParticle.cellN != MovedParticle.cell && !inside){
				k = MovedParticle.cellN;

				vv = Config.head[k];
			    
				while(vv != -1 && !inside){
					newvv = vv;
					if(newvv >= Config.Nc ) newvv -= Config.Nc; 
					if(v < newvv && !Config.part[newvv].already) overlapPear(); 

					vv = Config.link[vv];
				}
				NCell.push_back(k);
				Config.usedCell[k]=true;
			}

			if(!inside){
				if(dsx < 0){
					sx[0] = Config.W[0]-1;
					sx[1] = sx[0] -1;
					sx[2] = 0;
				}else{
					sx[0] = dsx*Config.w[0];
					if(sx[0]>=Config.W[0]) sx[0]-=Config.W[0];

					sx[1] = sx[0]-1;
					if(sx[1] < 0 ) sx[1]+=Config.W[0];

					sx[2] = sx[0]+1;
					if(sx[2]>=Config.W[0]) sx[2]-=Config.W[0];
				}

				if(dsy < 0){
					sy[0] = Config.W[1]-1;
					sy[1] = sy[0] -1;
					sy[2] = 0;
				}else{
					sy[0] = dsy*Config.w[1];
					if(sy[0]>=Config.W[1]) sy[0]-=Config.W[1];

					sy[1] = sy[0]-1;
					if(sy[1] < 0 ) sy[1]+=Config.W[1];

					sy[2] = sy[0]+1;
					if(sy[2]>=Config.W[1]) sy[2]-=Config.W[1];
				}

				if(dsz < 0){
					sz[0] = Config.W[2]-1;
					sz[1] = sz[0] -1;
					sz[2] = 0;
				}else{
					sz[0] = dsz*Config.w[2];
					if(sz[0]>=Config.W[2]) sz[0]-=Config.W[2];

					sz[1] = sz[0]-1;
					if(sz[1] < 0 ) sz[1]+=Config.W[2];

					sz[2] = sz[0]+1;
					if(sz[2]>=Config.W[2]) sz[2]-=Config.W[2];
				}

				for (int ix=0; ix < 3 && !inside; ix++){

					for (int iy=0; iy < 3 && !inside; iy++){

						for (int iz=0; iz <= 3 && !inside; iz++){
							k = sx[ix] + Config.W[0]*(sy[iy] + Config.W[1]*sz[iz]);
							if( Config.usedCell[k] ) continue;
			
							vv = Config.head[k];
					    
							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= Config.Nc ) newvv -= Config.Nc; 
								if(v < newvv && !Config.part[newvv].already) overlapPear(); 
								vv = Config.link[vv];
							}
							NCell.push_back(k);
							Config.usedCell[k]=true;
						}
					}
				}
			}


		
			if(!inside){
				first=false;
				if(dsxN < 0){
					sxN[0] = Config.W[0]-1;
					sxN[1] = sxN[0] -1;
					sxN[2] = 0;
				}else{
					sxN[0] = dsxN*Config.w[0];
					if(sxN[0]>=Config.W[0]) sxN[0]-=Config.W[0];

					sxN[1] = sxN[0]-1;
					if(sxN[1] < 0 ) sxN[1]+=Config.W[0];

					sxN[2] = sxN[0]+1;
					if(sxN[2]>=Config.W[0]) sxN[2]-=Config.W[0];
				}

				if(dsyN < 0){
					syN[0] = Config.W[1]-1;
					syN[1] = syN[0] -1;
					syN[2] = 0;
				}else{
					syN[0] = dsyN*Config.w[1];
					if(syN[0]>=Config.W[1]) syN[0]-=Config.W[1];

					syN[1] = syN[0]-1;
					if(syN[1] < 0 ) syN[1]+=Config.W[1];

					syN[2] = syN[0]+1;
					if(syN[2]>=Config.W[1]) syN[2]-=Config.W[1];
				}

				if(dszN < 0){
					szN[0] = Config.W[2]-1;
					szN[1] = szN[0] -1;
					szN[2] = 0;
				}else{
					szN[0] = dszN*Config.w[2];
					if(szN[0]>=Config.W[2]) szN[0]-=Config.W[2];

					szN[1] = szN[0]-1;
					if(szN[1] < 0 ) szN[1]+=Config.W[2];

					szN[2] = szN[0]+1;
					if(szN[2]>=Config.W[2]) szN[2]-=Config.W[2];
				}
			
				if(sx[1] != sxN[1] || sy[1] != syN[1] || sz[1] != szN[1] ){

					std::swap(dsx,dsxN);
					std::swap(dsy,dsyN);
					std::swap(dsz,dszN);

					for (int ix=0; ix < 3 && !inside; ix++){

						for (int iy=0; iy < 3 && !inside; iy++){

							for (int iz=0; iz <= 3 && !inside; iz++){
								k = sxN[ix] + Config.W[0]*(syN[iy] + Config.W[1]*szN[iz]);
								if( Config.usedCell[k] ) continue;
				
								vv = Config.head[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.Nc ) newvv -= Config.Nc; 
									if(v < newvv && !Config.part[newvv].already) overlapPear(); 
									vv = Config.link[vv];
								}
							}
						}
					}
				}
			}
			for( vv = 0; vv < NPart.size(); vv++) Config.part[NPart[vv]].already=false;
			for( vv = 0; vv < NCell.size(); vv++) Config.usedCell[NCell[vv]]=false;
		}else{
			for( vv = 0; vv < Config.Nc && !inside; vv++ ){
				if(v==vv) continue;
				if(overlapP(Config.part[v],Config.part[vv], Config.l)) inside=true; 
			}

			for( vv = Config.Nc; vv < Config.part.size() && !inside; vv++ ){
				if(overlapPSPH(Config.part[v],Config.part[vv], Config.l)) inside=true;
			}
		}
	}else{
		for( vv = 0; vv < Config.Nc && !inside; vv++ ){
			if(overlapPSPH(Config.part[vv],Config.part[v], Config.l)) inside=true; 
		}

		if(!inside){

			int k = Config.part[v].cell;

			vv = Config.head[k];
	    
			while(vv != -1 && !inside){
				if(v != vv ) if(overlapSPH(Config.part[vv],Config.part[v], Config.l)) inside=true; 

				vv = Config.link[vv];
			}
		

			sx[1] = Config.part[v].pos[0]*Config.w[0];
			if(sx[1]>=Config.W[0]-1){ 
				sx[1] = Config.W[0]-1;
				sx[0] = Config.W[0]-2;
				sx[2] = 0;
			}else{
				sx[2]=sx[1]+1;
				if(sx[1]==0) sx[0] = Config.W[0]-1;
				else sx[0] = sx[1] -1;
			}

			sy[1] = Config.part[v].pos[1]*Config.w[1];
			if(sy[1]>=Config.W[1]-1){ 
				sy[1] = Config.W[1]-1;
				sy[0] = Config.W[1]-2;
				sy[2] = 0;
			}else{
				sy[2]=sy[1]+1;
				if(sy[1]==0) sy[0] = Config.W[1]-1;
				else sy[0] = sy[1] -1;
			}

			sz[1] = Config.part[v].pos[2]*Config.w[2];
			if(sz[1]>=Config.W[2]-1){ 
				sz[1] = Config.W[2]-1;
				sz[0] = Config.W[2]-2;
				sz[2] = 0;
			}else{
				sz[2]=sz[1]+1;
				if(sz[1]==0) sz[0] = Config.W[2]-1;
				else sz[0] = sz[1] -1;
			}

			for (int ix=0; ix < 3 && !inside; ix++){

				for (int iy=0; iy < 3 && !inside; iy++){

					for (int iz=0; iz <= 3 && !inside; iz++){
						k = sx[ix] + Config.W[0]*(sy[iy] + Config.W[1]*sz[iz]);
						if( k == Config.part[v].cell ) continue;
		
						vv = Config.head[k];
				    
						while(vv != -1 && !inside){
							if(overlapSPH(Config.part[vv],Config.part[v], Config.l)) inside=true; 
							vv = Config.link[vv];
						}
					}
				}
			}
		}
        }
    }

    if(!inside){ 
	if(Compressing){
		Config.Vbox = Vn;
		Config.rhoN = N/Config.Vbox; 
		Config.rhoV = Config.rhoN*Config.Vsys; 
		vproc = 0.5*(3*vproc-1);

		if(vproc < 0.9 ) vproc=0.9;
		std::cout << "Compession successful! New rho = " << Config.rhoV << " and new vproc=" << vproc << std::endl;
		std::cout << "(lx,ly,lz) = (" << Config.l[0] << "," << Config.l[1] << "," << Config.l[2] << ")"<< std::endl;

		Config.write("Save/Config.dat",1);

		if( Config.rhoV - rhoInit < 0.01 ){
			Config.write("Results/Config" + toString(Config.rhoV) + ".dat",0);
			rhoInit = Config.rhoV;
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
        Config.l[0] *= ln[0];
        Config.l[1] *= ln[1];
        Config.l[2] *= ln[2];

	Config.l_2[0] = Config.l[0]*0.5;
	Config.l_2[1] = Config.l[1]*0.5;
	Config.l_2[2] = Config.l[2]*0.5;

        for( int i=0; i < Config.part.size(); i++){
                Config.part[i].pos[0] *= ln[0];
                Config.part[i].pos[1] *= ln[1];
        	Config.part[i].pos[2] *= ln[2];
        }

	if(Compressing){
		vproc = 0.5*(1+vproc);
		if(vproc > 0.9999){ 
		    vproc = 0.9999;
		    initCompress = true;
		}
		std::cout << "Compession unsuccessful! New vproc=" << vproc << std::endl;
		std::cout << "(lx,ly,lz) = (" << Config.l[0] << "," << Config.l[1] << "," << Config.l[2] << ")"<< std::endl;
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
		RenewList();
	}
    }
}
