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
                		phi = 2*M_PI*uni(gen);
                        	theta = acos(2*uni(gen) - 1);				
				cosphi = uni(gen)*pos_lambda;

				costheta = cos(theta)*sin(phi)*cosphi;
				sintheta = sin(theta)*sin(phi)*cosphi;
				cosphi = cos(phi)*cosphi;
                		MovedParticle.pos[0] = MovedParticle.pos[0] + costheta;
                		MovedParticle.pos[1] = MovedParticle.pos[1] + sintheta;
                		MovedParticle.pos[2] = MovedParticle.pos[2] + cosphi;

				MovedParticle.pos_msd[0] += (MovedParticle.pos[0]-Config.part[randP].pos[0]);
				MovedParticle.pos_msd[1] += (MovedParticle.pos[1]-Config.part[randP].pos[1]);
				MovedParticle.pos_msd[2] += (MovedParticle.pos[2]-Config.part[randP].pos[2]);

            		}else{
                		phi = (uni(gen)-0.5)*ori_lambda;
                		cosphi = cos(phi);
                		sinphi = sin(phi);

                		Case = 3*uni(gen);
                		switch(Case){
                			case 0:
                    			MovedParticle.ori[0] = cosphi*Config.part[randP].ori[0] + sinphi*Config.part[randP].ori[1];
                    			MovedParticle.ori[1] = -sinphi*Config.part[randP].ori[0] + cosphi*Config.part[randP].ori[1];
                    			break;

                			case 1:
                    			MovedParticle.ori[2] = cosphi*Config.part[randP].ori[2] + sinphi*Config.part[randP].ori[0];
                    			MovedParticle.ori[0] = -sinphi*Config.part[randP].ori[2] + cosphi*Config.part[randP].ori[0];
                    			break;

                			case 2:
                    			MovedParticle.ori[1] = cosphi*Config.part[randP].ori[1] + sinphi*Config.part[randP].ori[2];
                    			MovedParticle.ori[2] = -sinphi*Config.part[randP].ori[1] + cosphi*Config.part[randP].ori[2];
                    			break;
                		}

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


			NCell.resize(0);
			NPart.resize(0);

			k = MovedParticle.cell;

		    	vv = Config.headP[k];
		    
		    	while(vv != -1 && !inside){
				newvv = vv;
				if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
				if(randP != newvv && !Config.part[newvv].already){ 
					NPart.push_back(newvv);
					Config.part[newvv].already=true;	

					R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
					R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
					R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
					if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
					else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
					if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
					else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
					if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
					else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
					overlapPear(); 
				}

				vv = Config.linkP[vv];
		    	}
			NCell.push_back(k);
			Config.usedCell[k]=true;

			if(Config.Ns != 0 ){

				newvv = Config.headPS[k];
			    
				while(newvv != -1 && !inside){
					R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
					R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
					R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
					if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
					else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
					if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
					else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
					if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
					else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
					
					overlapPearSphere(); 
					newvv = Config.linkP[newvv];
				}
			}

			if( MovedParticle.cellN != MovedParticle.cell && !inside){

				k = MovedParticle.cellN;

				vv = Config.headP[k];
			    
				while(vv != -1 && !inside){
					newvv = vv;
					if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
					if(randP != newvv && !Config.part[newvv].already){ 
						NPart.push_back(newvv);
						Config.part[newvv].already=true;
						R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
						R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
						R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
						if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
						else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
						if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
						else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
						if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
						else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
						overlapPear(); 
					}

					vv = Config.linkP[vv];
				}
				NCell.push_back(k);
				Config.usedCell[k]=true;

				if(Config.Ns != 0 ){

					newvv = Config.headPS[k];
				    
					while(newvv != -1 && !inside){
						R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
						R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
						R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
						if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
						else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
						if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
						else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
						if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
						else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
						overlapPearSphere(); 
						newvv = Config.linkP[newvv];
					}
				}
			}


			if(!inside){
				if( MovedParticle.pos[0] > Config.l[0] ) MovedParticle.pos[0] -= Config.l[0];
				else if( MovedParticle.pos[0] < 0 ) MovedParticle.pos[0] += Config.l[0];
			
				if( MovedParticle.pos[1] > Config.l[1] ) MovedParticle.pos[1] -= Config.l[1];
				else if( MovedParticle.pos[1] < 0 ) MovedParticle.pos[1] += Config.l[1];
	
				if( MovedParticle.pos[2] > Config.l[2] ) MovedParticle.pos[2] -= Config.l[2];
				else if( MovedParticle.pos[2] < 0 ) MovedParticle.pos[2] += Config.l[2];
	

				dist_ori[0] = distN*MovedParticle.ori[0];
				dist_ori[1] = distN*MovedParticle.ori[1];
				dist_ori[2] = distN*MovedParticle.ori[2];

				dsx = MovedParticle.pos[0]-dist_ori[0];
				dsy = MovedParticle.pos[1]-dist_ori[1];
				dsz = MovedParticle.pos[2]-dist_ori[2];

				in_test = false;
				if(dsx < 0){
					sx[0] = Config.WP[0]-1;
					sx[1] = sx[0] -1;
					sx[2] = 0;
					in_test = true;
				}else{
					if(dsx > Config.l[0]){
						sx[0] = 0;
						sx[1] = Config.WP[0]-1;
						sx[2] = 1;
						in_test = true;
					}else{
						sx[0] = dsx*Config.wP[0];

						if( sx[0] >= Config.WP[0]-2 || sx[0] <= 1){ 
							in_test = true;

							sx[1] = sx[0]-1;
							if(sx[1] < 0 ) sx[1] = Config.WP[0]-1;

							sx[2] = sx[0]+1;
							if(sx[2]> Config.WP[0]-1) sx[2] = 0;
						}else{
							sx[1] = sx[0]-1;
							sx[2] = sx[0]+1;
						}
					}
				}
			
				

				if(dsy < 0){
					sy[0] = Config.WP[1]-1;
					sy[1] = sy[0] -1;
					sy[2] = 0;
					in_test = true;
				}else{
					if(dsy > Config.l[1]){
						sy[0] = 0;
						sy[1] = Config.WP[1]-1;
						sy[2] = 1;
						in_test = true;
					}else{
						sy[0] = dsy*Config.wP[1];
						if( sy[0] >= Config.WP[1]-2 || sy[0] <= 1){ 
							in_test = true;

							sy[1] = sy[0]-1;
							if(sy[1] < 0 ) sy[1] = Config.WP[1]-1;

							sy[2] = sy[0]+1;
							if(sy[2]> Config.WP[1]-1) sy[2] = 0;
						}else{
							sy[1] = sy[0]-1;
							sy[2] = sy[0]+1;
						}
					}
				}

				if(dsz < 0){
					sz[0] = Config.WP[2]-1;
					sz[1] = sz[0] -1;
					sz[2] = 0;
					in_test = true;
				}else{
					if(dsz > Config.l[2]){
						sz[0] = 0;
						sz[1] = Config.WP[2]-1;
						sz[2] = 1;
						in_test = true;
					}else{
						sz[0] = dsz*Config.wP[2];
						if( sz[0] >= Config.WP[2]-2 || sz[0] <= 1){ 
							in_test = true;

							sz[1] = sz[0]-1;
							if(sz[1] < 0 ) sz[1] = Config.WP[2]-1;

							sz[2] = sz[0]+1;
							if(sz[2]> Config.WP[2]-1) sz[2] = 0;
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
								k = sx[ix] + Config.WP[0]*(sy[iy]+ Config.WP[1]*sz[iz]);
								if( Config.usedCell[k] ) continue;
					
								vv = Config.headP[k];

								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
									if(!Config.part[newvv].already){ 
										NPart.push_back(newvv);
										Config.part[newvv].already=true;	
										
										R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
										if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
										else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
										if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
										else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
										if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
										else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
										overlapPear(); 
									}
									vv = Config.linkP[vv];
								}
								NCell.push_back(k);
								Config.usedCell[k]=true;

								if(Config.Ns != 0 ){

									newvv = Config.headPS[k];

									lb[0]= 0;
									lb[1]= 0;
									lb[2]= 0;

									R[0] = sx[ix]/Config.wP[0] - MovedParticle.pos[0];
									if(R[0] > Config.l_2[0]) lb[0] = -Config.l[0];
									else if (R[0] < -Config.l_2[0])  lb[0] = Config.l[0];

									R[1] = sy[iy]/Config.wP[1] - MovedParticle.pos[1];
									if(R[1] > Config.l_2[1]) lb[1] = -Config.l[1];
									else if (R[1] < -Config.l_2[1]) lb[1] = Config.l[1];
							    
									R[2] = sy[iy]/Config.wP[2] - MovedParticle.pos[2];
									if(R[2] > Config.l_2[2]) lb[2] = -Config.l[2];
									else if (R[2] < -Config.l_2[2]) lb[2] = Config.l[2];
								    
									while(newvv != -1 && !inside){
										R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0] + lb[0];
										R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1] + lb[1];
										R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2] + lb[2];
										if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0]) std::cout << "PSb("  << sx[ix] << ") (" << sx[0] << ") " << R[0] << " " << Config.l[0] << " (" << MovedParticle.pos[0] << "," <<  Config.part[newvv].pos[0] << "," << Config.WP[0] << "," << Config.wP[0] << ")" << std::endl;
										if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1]) std::cout << "PSb(" << sy[iy] << ") (" << sy[1] << ") " << R[1] << " " << Config.l[1] << " (" << MovedParticle.pos[1] << "," <<  Config.part[newvv].pos[1] << "," << sy[iy]*Config.wP[1] << "," << Config.wP[0] << ")" << std::endl;
										if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2]) std::cout << "PSb(" << sz[iz] << ") (" << sz[2] << ") " << R[2] << " " << Config.l[2] << " (" << MovedParticle.pos[2] << "," <<  Config.part[newvv].pos[2] << "," << sz[iz]*Config.wP[2] << "," << Config.wP[0] << ")" << std::endl;
										overlapPearSphere(); 
										newvv = Config.linkP[newvv];
									}
								}
							}
						}
					}
				}else{
					for (int ix=0; ix < 3 && !inside; ix++){
						for (int iy=0; iy < 3 && !inside; iy++){
							for (int iz=0; iz < 3 && !inside; iz++){
								k = sx[ix] + Config.WP[0]*(sy[iy]+ Config.WP[1]*sz[iz]);
								if( Config.usedCell[k] ) continue;
					
								vv = Config.headP[k];

								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
									if(!Config.part[newvv].already){ 
										NPart.push_back(newvv);
										Config.part[newvv].already=true;	
										
										R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];

										if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "PPb(" << sx[ix] << "," << sy[iy] << ") (" << sx[0] << "," << sy[0] << ") " << R[0] << " " << Config.l_2[0] << " (" << bx[ix] << ",0," << lb[0] << " )" << std::endl;
										if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "PPb(" << sx[ix] << "," << sy[iy] << ") (" << sx[0] << "," << sy[0] << ") " << R[1] << " " << Config.l_2[1] << " (" << by[iy] << ",1," << lb[1] << " )" <<std::endl;
										if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "PPb(" << sx[ix] << "," << sz[iz] << ") (" << sx[0] << "," << sz[0] << ") " << R[2] << " " << Config.l_2[2] << " (" << bz[iz] << ",1," << lb[2] << " )" <<std::endl;
										overlapPear(); 
									}
									vv = Config.linkP[vv];
								}
								NCell.push_back(k);
								Config.usedCell[k]=true;

								if(Config.Ns != 0 ){

									newvv = Config.headPS[k];

									while(newvv != -1 && !inside){
										R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
										if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0]) std::cout << "PSb("  << sx[ix] << ") (" << sx[0] << ") " << R[0] << " " << Config.l[0] << " (" << MovedParticle.pos[0] << "," <<  Config.part[newvv].pos[0] << "," << Config.WP[0] << "," << Config.wP[0] << ")" << std::endl;
										if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1]) std::cout << "PSb(" << sy[iy] << ") (" << sy[1] << ") " << R[1] << " " << Config.l[1] << " (" << MovedParticle.pos[1] << "," <<  Config.part[newvv].pos[1] << "," << sy[iy]*Config.wP[1] << "," << Config.wP[0] << ")" << std::endl;
										if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2]) std::cout << "PSb(" << sz[iz] << ") (" << sz[2] << ") " << R[2] << " " << Config.l[2] << " (" << MovedParticle.pos[2] << "," <<  Config.part[newvv].pos[2] << "," << sz[iz]*Config.wP[2] << "," << Config.wP[0] << ")" << std::endl;
										overlapPearSphere(); 
										newvv = Config.linkP[newvv];
									}
								}
							}
						}
					}
				}
			}
		
			if(!inside){

				dsxN = MovedParticle.pos[0]+dist_ori[0];

				dsyN = MovedParticle.pos[1]+dist_ori[1];

				dszN = MovedParticle.pos[2]+dist_ori[2];

				in_test=false;
				if(dsxN < 0){
					sxN[0] = Config.WP[0]-1;
					in_test = true;
				}else{
					if(dsxN > Config.l[0]){
						sxN[0] = 0;
						in_test = true;
					}else{
						sxN[0] = dsxN*Config.wP[0];
						if( sxN[0] >= Config.WP[0]-2 || sxN[0] <= 1) in_test = true;
					}
				}

				if(dsyN < 0){
					syN[0] = Config.WP[1]-1;
					in_test = true;
				}else{
					if(dsyN > Config.l[1]){
						syN[0] = 0;
						in_test = true;
					}else{
						syN[0] = dsyN*Config.wP[1];
						if( syN[0] >= Config.WP[1]-2 || syN[0] <= 1) in_test = true;
					}
				}

				if(dszN < 0){
					szN[0] = Config.WP[2]-1;
					in_test = true;
				}else{
					if(dszN > Config.l[2]){
						szN[0] = 0;
						in_test = true;
					}else{
						szN[0] = dszN*Config.wP[2];
						if( szN[0] >= Config.WP[2]-2 || szN[0] <= 1) in_test = true;
					}
				}

				if(sx[0] != sxN[0] || sy[0] != syN[0] || sz[0] != szN[0] ){

					sxN[1] = sxN[0]-1;
					if(sxN[1] < 0 ) sxN[1] = Config.WP[0]-1;

					sxN[2] = sxN[0]+1;
					if(sxN[2]> Config.WP[0]-1) sxN[2] = 0;

					syN[1] = syN[0]-1;
					if(syN[1] < 0 ) syN[1] = Config.WP[1]-1;

					syN[2] = syN[0]+1;
					if(syN[2]> Config.WP[1]-1) syN[2] = 0;
					
					szN[1] = szN[0]-1;
					if(szN[1] < 0 ) szN[1] = Config.WP[2]-1;

					szN[2] = szN[0]+1;
					if(szN[2]> Config.WP[2]-1) szN[2] = 0;
					
					if(in_test){
						for (int ix=0; ix < 3 && !inside; ix++){

							for (int iy=0; iy < 3 && !inside; iy++){
								for (int iz=0; iz < 3 && !inside; iz++){

									k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
									if( Config.usedCell[k] ) continue;
					
									vv = Config.headP[k];
							    
									while(vv != -1 && !inside){
										newvv = vv;
										if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
										if(!Config.part[newvv].already){ 
											NPart.push_back(newvv);
											Config.part[newvv].already=true;	

											R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
											R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
											R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
											if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
											else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
											if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
											else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
											if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
											else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
											overlapPear(); 
										}
										vv = Config.linkP[vv];
									}
									if(Config.Ns != 0 ){

										newvv = Config.headPS[k];

										lb[0]= 0;
										lb[1]= 0;
										lb[2]= 0;

										R[0] = sxN[ix]/Config.wP[0] - MovedParticle.pos[0];
										if(R[0] > Config.l_2[0]) lb[0] = -Config.l[0];
										else if (R[0] < -Config.l_2[0])  lb[0] = Config.l[0];

										R[1] = syN[iy]/Config.wP[1] - MovedParticle.pos[1];
										if(R[1] > Config.l_2[1]) lb[1] = -Config.l[1];
										else if (R[1] < -Config.l_2[1])  lb[1] = Config.l[1];

										R[2] = syN[iy]/Config.wP[2] - MovedParticle.pos[2];
										if(R[2] > Config.l_2[2]) lb[2] = -Config.l[2];
										else if (R[2] < -Config.l_2[2])  lb[2] = Config.l[2];
									    
										while(newvv != -1 && !inside){
											R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0] + lb[0];
											R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1] + lb[1];
											R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2] + lb[2];
											if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "PSt(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[0] << " " << Config.l_2[0] << " (" << bx[ix] << ",0," << lb[0] << " )" << randP << " " << newvv << std::endl;
											if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "PSt(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[1] << " " << Config.l_2[1] << " (" << by[iy] << ",1," << lb[1] << " )" << randP << " " << newvv << std::endl;
											if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "PSt(" << sxN[ix] << "," << szN[iz] << ") (" << sxN[0] << "," << szN[0] << ") " << R[2] << " " << Config.l_2[2] << " (" << bz[iz] << ",1," << lb[2] << " )" << randP << " " << newvv << std::endl;
											overlapPearSphere(); 
											newvv = Config.linkP[newvv];
										}
									}
								}
							}
						}
					}else{
						for (int ix=0; ix < 3 && !inside; ix++){

							for (int iy=0; iy < 3 && !inside; iy++){
								for (int iz=0; iz < 3 && !inside; iz++){

									k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
									if( Config.usedCell[k] ) continue;
					
									vv = Config.headP[k];
							    
									while(vv != -1 && !inside){
										newvv = vv;
										if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
										if(!Config.part[newvv].already){ 
											NPart.push_back(newvv);
											Config.part[newvv].already=true;	

											R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
											R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
											R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
											if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "PPt(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[0] << " " << Config.l_2[0] << " (" << bx[ix] << ",0," << lb[0] << " )" << std::endl;
											if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "PPt(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[1] << " " << Config.l_2[1] << " (" << by[iy] << ",1," << lb[1] << " )" <<std::endl;
											if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "PPt(" << sxN[ix] << "," << szN[iz] << ") (" << sxN[0] << "," << szN[0] << ") " << R[2] << " " << Config.l_2[2] << " (" << bz[iz] << ",1," << lb[2] << " )" <<std::endl;
											overlapPear(); 
										}
										vv = Config.linkP[vv];
									}
									if(Config.Ns != 0 ){

										newvv = Config.headPS[k];

										while(newvv != -1 && !inside){
											R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
											R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
											R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
											if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "PSt(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[0] << " " << Config.l_2[0] << " (" << bx[ix] << ",0," << lb[0] << " )" << randP << " " << newvv << std::endl;
											if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "PSt(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[1] << " " << Config.l_2[1] << " (" << by[iy] << ",1," << lb[1] << " )" << randP << " " << newvv << std::endl;
											if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "PSt(" << sxN[ix] << "," << szN[iz] << ") (" << sxN[0] << "," << szN[0] << ") " << R[2] << " " << Config.l_2[2] << " (" << bz[iz] << ",1," << lb[2] << " )" << randP << " " << newvv << std::endl;
											overlapPearSphere(); 
											newvv = Config.linkP[newvv];
										}
									}
								}
							}
						}
					}
				}
			}

			for( vv = 0; vv < NPart.size(); vv++) Config.part[NPart[vv]].already=false;
			for( vv = 0; vv < NCell.size(); vv++) Config.usedCell[NCell[vv]]=false;

        	}else{
			phi = 2*M_PI*uni(gen);
			theta = acos(2*uni(gen) - 1);				
			cosphi = uni(gen)*pos_lambda;

			costheta = cos(theta)*sin(phi)*cosphi;
			sintheta = sin(theta)*sin(phi)*cosphi;
			cosphi = cos(phi)*cosphi;
			MovedParticle.pos[0] = MovedParticle.pos[0] + costheta;
			MovedParticle.pos[1] = MovedParticle.pos[1] + sintheta;
			MovedParticle.pos[2] = MovedParticle.pos[2] + cosphi;

			MovedParticle.pos_msd[0] += (MovedParticle.pos[0]-Config.part[randP].pos[0]);
			MovedParticle.pos_msd[1] += (MovedParticle.pos[1]-Config.part[randP].pos[1]);
			MovedParticle.pos_msd[2] += (MovedParticle.pos[2]-Config.part[randP].pos[2]);

			if( MovedParticle.pos[0] > Config.l[0] ) MovedParticle.pos[0] -= Config.l[0];
			else if( MovedParticle.pos[0] < 0 )MovedParticle.pos[0] += Config.l[0];
			if( MovedParticle.pos[1] > Config.l[1] ) MovedParticle.pos[1] -= Config.l[1];
			else if( MovedParticle.pos[1] < 0 ) MovedParticle.pos[1] += Config.l[1];
			if( MovedParticle.pos[2] > Config.l[2] ) MovedParticle.pos[2] -= Config.l[2];
			else if( MovedParticle.pos[2] < 0 ) MovedParticle.pos[2] += Config.l[2];
	
			sx[0] = MovedParticle.pos[0]*Config.wS[0];
			sy[0] = MovedParticle.pos[1]*Config.wS[1];
			sz[0] = MovedParticle.pos[2]*Config.wS[2];

			kN = sx[0] + Config.WS[0]*(sy[0]+Config.WS[1]*sz[0]);
		
			if(Config.headS[kN] != -1 && Config.headS[kN] != randP) inside = true;

			if(!inside ){
				NPart.resize(0);
				
				k = MovedParticle.cellN;

				vv = Config.headP[k];
		    
				while(vv != -1 && !inside){
					newvv = vv;
					if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
					if(!Config.part[newvv].already){ 
						NPart.push_back(newvv);
						Config.part[newvv].already=true;	

						R[0] = MovedParticle.pos[0] - Config.part[newvv].pos[0];
						R[1] = MovedParticle.pos[1] - Config.part[newvv].pos[1];
						R[2] = MovedParticle.pos[2] - Config.part[newvv].pos[2];
						if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
						else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
						if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
						else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
						if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
						else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
						overlapSpherePear(); 
					}

					vv = Config.linkP[vv];
				}

				if(!inside){
					in_test=false;
				
					sxN[0] = MovedParticle.pos[0]*Config.wP[0];
					if( sxN[0] >= Config.WP[0]-2 || sxN[0] <= 1){ 
						in_test = true;

						sxN[1] = sxN[0]-1;
						if(sxN[1] < 0 ) sxN[1] = Config.WP[0]-1;

						sxN[2] = sxN[0]+1;
						if(sxN[2]> Config.WP[0]-1) sxN[2] = 0;
					}else{
						sxN[1] = sxN[0]-1;

						sxN[2] = sxN[0]+1;
					}

					syN[0] = MovedParticle.pos[1]*Config.wP[1];
					if( syN[0] >= Config.WP[1]-2 || syN[0] <= 1){ 
						in_test = true;

						syN[1] = syN[0]-1;
						if(syN[1] < 0 ) syN[1] = Config.WP[1]-1;

						syN[2] = syN[0]+1;
						if(syN[2]> Config.WP[1]-1) syN[2] = 0;
					}else{
						syN[1] = syN[0]-1;

						syN[2] = syN[0]+1;
					}

					szN[0] = MovedParticle.pos[2]*Config.wP[2];
					if( szN[0] >= Config.WP[2]-2 || szN[0] <= 1){ 
						in_test = true;

						szN[1] = szN[0]-1;
						if(szN[1] < 0 ) szN[1] = Config.WP[2]-1;

						szN[2] = szN[0]+1;
						if(szN[2]> Config.WP[2]-1) szN[2] = 0;
					}else{
						szN[1] = szN[0]-1;

						szN[2] = szN[0]+1;
					}

					if(in_test){
						for (int ix=0; ix < 3 && !inside; ix++){

							for (int iy=0; iy < 3 && !inside; iy++){
								for (int iz=0; iz < 3 && !inside; iz++){

									k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
									if( k == MovedParticle.cellN ) continue;
						
									vv = Config.headP[k];
								    
									while(vv != -1 && !inside){
										newvv = vv;
										if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
										if(!Config.part[newvv].already){ 
											NPart.push_back(newvv);
											Config.part[newvv].already=true;	

											R[0] = MovedParticle.pos[0] - Config.part[newvv].pos[0];
											R[1] = MovedParticle.pos[1] - Config.part[newvv].pos[1];
											R[2] = MovedParticle.pos[2] - Config.part[newvv].pos[2];
											if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
											else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
											if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
											else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
											if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
											else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];

											overlapSpherePear(); 
										}
										vv = Config.linkP[vv];
									}
								}
							}
						}
					}else{
						for (int ix=0; ix < 3 && !inside; ix++){
							for (int iy=0; iy < 3 && !inside; iy++){
								for (int iz=0; iz < 3 && !inside; iz++){

									k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
									if( k == MovedParticle.cellN ) continue;
						
									vv = Config.headP[k];
								    
									while(vv != -1 && !inside){
										newvv = vv;
										if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
										if(!Config.part[newvv].already){ 
											NPart.push_back(newvv);
											Config.part[newvv].already=true;	

											R[0] = MovedParticle.pos[0] - Config.part[newvv].pos[0];
											R[1] = MovedParticle.pos[1] - Config.part[newvv].pos[1];
											R[2] = MovedParticle.pos[2] - Config.part[newvv].pos[2];
											if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "SP(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[0] << " " << Config.l_2[0] << " (" << bx[ix] << ",0," << lb[0] << "," << MovedParticle.pos[0] << "," << Config.part[newvv].pos[0] << ")" << std::endl;
											if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "SP(" << sxN[ix] << "," << syN[iy] << ") (" << sxN[0] << "," << syN[0] << ") " << R[1] << " " << Config.l_2[1] << " (" << by[iy] << ",1," << lb[1] << "," << MovedParticle.pos[0] << "," << Config.part[newvv].pos[1] << ")" <<std::endl;
											if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "SP(" << sxN[ix] << "," << szN[iz] << ") (" << sxN[0] << "," << szN[0] << ") " << R[2] << " " << Config.l_2[2] << " (" << bz[iz] << ",1," << lb[2] << "," << MovedParticle.pos[0] << "," << Config.part[newvv].pos[2] << ")" <<std::endl;

											overlapSpherePear(); 
										}
										vv = Config.linkP[vv];
									}
								}
							}
						}
					}
				}

				for( vv = 0; vv < NPart.size(); vv++) Config.part[NPart[vv]].already=false;
			}

	    		if(!inside){
				in_test=false;

				if( sx[0] >= Config.WS[0]-2 || sx[0] <= 1){ 
					in_test = true;
					sx[1] = sx[0]-1;
					if(sx[1] < 0 ){ 
						sx[1] = Config.WS[0]-1;
						bx[1]=-1;
						bx[3]=-1;
					}

					sx[2] = sx[0]+1;
					if(sx[2]> Config.WS[0]-1){ 
						sx[2] = 0;
						bx[2]=1;
						bx[4]=1;
					}

					sx[3] = sx[1]-1;
					if(sx[3] < 0 ){ 
						sx[3] = Config.WS[0]-1;
						bx[3]=-1;
					}

					sx[4] = sx[2]+1;
					if(sx[4]> Config.WS[0]-1){ 
						sx[4] = 0;
						bx[4]=1;
					}
				}else{
					sx[1] = sx[0]-1;
					sx[2] = sx[0]+1;
					sx[3] = sx[1]-1;
					sx[4] = sx[2]+1;
				}

				if( sy[0] >= Config.WS[1]-2 || sy[0] <= 1){ 
					in_test = true;

					sy[1] = sy[0]-1;
					if(sy[1] < 0 ){ 
						sy[1] = Config.WS[1]-1;
						by[1]=-1;
						by[3]=-1;
					}

					sy[2] = sy[0]+1;
					if(sy[2]> Config.WS[1]-1){ 
						sy[2] = 0;
						by[2]=1;
						by[4]=1;
					}

					sy[3] = sy[1]-1;
					if(sy[3] < 0 ){ 
						sy[3] = Config.WS[1]-1;
						by[3]=-1;
					}

					sy[4] = sy[2]+1;
					if(sy[4]> Config.WS[1]-1){ 
						sy[4] = 0;
						by[4]=1;
					}
				}else{
					sy[1] = sy[0]-1;
					sy[2] = sy[0]+1;
					sy[3] = sy[1]-1;
					sy[4] = sy[2]+1;
				}

				if( sz[0] >= Config.WS[2]-2 || sz[0] <= 1){ 
					in_test = true;

					sz[1] = sz[0]-1;
					if(sz[1] < 0 ){ 
						sz[1] = Config.WS[2]-1;
						bz[1]=-1;
						bz[3]=-1;
					}

					sz[2] = sz[0]+1;
					if(sz[2]> Config.WS[2]-1){ 
						sz[2] = 0;
						bz[2]=1;
						bz[4]=1;
					}

					sz[3] = sz[1]-1;
					if(sz[3] < 0 ){ 
						sz[3] = Config.WS[2]-1;
						bz[3]=-1;
					}

					sz[4] = sz[2]+1;
					if(sz[4]> Config.WS[2]-1){ 
						sz[4] = 0;
						bz[4]=1;
					}
				}else{
					sz[1] = sz[0]-1;
					sz[2] = sz[0]+1;
					sz[3] = sz[1]-1;
					sz[4] = sz[2]+1;
				}

				if(in_test){
					for (int ix=0; ix < 5 && !inside; ix++){

						for (int iy=0; iy < 5 && !inside; iy++){
							for (int iz=0; iz < 5 && !inside; iz++){

								k = sx[ix] + Config.WS[0]*(sy[iy]+Config.WS[1]*sz[iz]);
								if( k == MovedParticle.cell || k == kN) continue;
				
								vv = Config.headS[k];

								lb[0]= bx[ix]*Config.l[0];
								lb[1]= by[iy]*Config.l[1];
								lb[2]= bz[iz]*Config.l[2];
						    
								if(vv != -1){
									R[0] = MovedParticle.pos[0] - Config.part[vv].pos[0] - lb[0];
									R[1] = MovedParticle.pos[1] - Config.part[vv].pos[1] - lb[1];
									R[2] = MovedParticle.pos[2] - Config.part[vv].pos[2] - lb[2];
									if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "1SS(" << sx[ix] << "," << sy[iy] << ") (" << sx[0] << "," << sy[0] << ") (" << bx[ix] << ",0," << lb[0] << ") (" << MovedParticle.pos[0] << "," << MovedParticle.pos[1] << ") (" << Config.part[vv].pos[0] << "," << Config.part[vv].pos[1] << ") " << k << " " << Config.part[vv].cell << std::endl;
									if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "1SS(" << sx[ix] << "," << sy[iy] << ") (" << sx[0] << "," << sy[0] << ") (" << bx[ix] << ",0," << lb[0] << ") (" << MovedParticle.pos[0] << "," << MovedParticle.pos[1] << ") (" << Config.part[vv].pos[0] << "," << Config.part[vv].pos[1] << ") " << k << " " << Config.part[vv].cell << std::endl;
									if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "1SS(" << sx[ix] << "," << sz[iz] << ") (" << sx[0] << "," << sz[0] << ") (" << bx[ix] << ",0," << lb[0] << ") (" << MovedParticle.pos[0] << "," << MovedParticle.pos[2] << ") (" << Config.part[vv].pos[0] << "," << Config.part[vv].pos[2] << ") " << k << " " << Config.part[vv].cell << std::endl;
									overlapSphere(); 
								}
							}
						}
					}
					bx[0] =0;
					bx[1] =0;
					bx[2] =0;
					bx[3] =0;
					bx[4] =0;

					by[0] =0;
					by[1] =0;
					by[2] =0;
					by[3] =0;
					by[4] =0;

					bz[0] =0;
					bz[1] =0;
					bz[2] =0;
					bz[3] =0;
					bz[4] =0;
				}else{
					for (int ix=0; ix < 5 && !inside; ix++){

						for (int iy=0; iy < 5 && !inside; iy++){
							for (int iz=0; iz < 5 && !inside; iz++){

								k = sx[ix] + Config.WS[0]*(sy[iy]+Config.WS[1]*sz[iz]);
								if( k == MovedParticle.cell || k == kN) continue;
				
								vv = Config.headS[k];

								if(vv != -1){
									R[0] = MovedParticle.pos[0] - Config.part[vv].pos[0];
									R[1] = MovedParticle.pos[1] - Config.part[vv].pos[1];
									R[2] = MovedParticle.pos[2] - Config.part[vv].pos[2];
									if (R[0] > Config.l_2[0] || R[0] < -Config.l_2[0])  std::cout << "2SS(" << sx[ix] << "," << sy[iy] << ") (" << sx[0] << "," << sy[0] << ") (" << bx[ix] << ",0," << lb[0] << ") (" << MovedParticle.pos[0] << "," << MovedParticle.pos[1] << ") (" << Config.part[vv].pos[0] << "," << Config.part[vv].pos[1] << ") " << k << " " << Config.part[vv].cell << std::endl;
									if (R[1] > Config.l_2[1] || R[1] < -Config.l_2[1])  std::cout << "2SS(" << sx[ix] << "," << sy[iy] << ") (" << sx[0] << "," << sy[0] << ") (" << bx[ix] << ",0," << lb[0] << ") (" << MovedParticle.pos[0] << "," << MovedParticle.pos[1] << ") (" << Config.part[vv].pos[0] << "," << Config.part[vv].pos[1] << ") " << k << " " << Config.part[vv].cell << std::endl;
									if (R[2] > Config.l_2[2] || R[2] < -Config.l_2[2])  std::cout << "2SS(" << sx[ix] << "," << sz[iz] << ") (" << sx[0] << "," << sz[0] << ") (" << bx[ix] << ",0," << lb[0] << ") (" << MovedParticle.pos[0] << "," << MovedParticle.pos[2] << ") (" << Config.part[vv].pos[0] << "," << Config.part[vv].pos[2] << ") " << k << " " << Config.part[vv].cell << std::endl;
									overlapSphere(); 
								}
							}
						}
					}
				}
			}
		}


		if(!inside){
	            	acc++;

        		if( randP < Config.Nc ){
				k = sx[0] + Config.WP[0]*(sy[0]+Config.WP[1]*sz[0]); 

				if(MovedParticle.cell != k){
					if(Config.headP[MovedParticle.cell] != randP ){
						int vvv = Config.headP[MovedParticle.cell];
						vv = Config.linkP[vvv];
						while(vv!=randP){ 
							vvv = vv;
							vv = Config.linkP[vvv];
						}
				
						Config.linkP[vvv] = Config.linkP[randP];
					}else{
						Config.headP[MovedParticle.cell] = Config.linkP[randP];
					}
		    
					Config.linkP[randP] = Config.headP[k];
					Config.headP[k] = randP;
					MovedParticle.cell = k;
				}

				k = sxN[0] + Config.WP[0]*(syN[0]+Config.WP[1]*szN[0]); 

				if(MovedParticle.cellN != k){
					if(Config.headP[MovedParticle.cellN] != randP+Config.part.size() ){
						int vvv = Config.headP[MovedParticle.cellN];
						vv = Config.linkP[vvv];
						while(vv!=randP+Config.part.size()){ 
							vvv = vv;
							vv = Config.linkP[vvv];
						}
				
						Config.linkP[vvv] = Config.linkP[randP+Config.part.size()];
					}else{
						Config.headP[MovedParticle.cellN] = Config.linkP[randP+Config.part.size()];
					}
		    
					Config.linkP[randP+Config.part.size()] = Config.headP[k];
					Config.headP[k] = randP+Config.part.size();
					MovedParticle.cellN = k;
				}


			}else{ 
				k = sx[0] + Config.WS[0]*(sy[0]+Config.WS[1]*sz[0]); 

				if(MovedParticle.cell != k){
					Config.headS[MovedParticle.cell] = -1;
					Config.headS[k] = randP;
					MovedParticle.cell = k;
				}
			   
				k = sxN[0] + Config.WP[0]*(syN[0]+Config.WP[1]*szN[0]); 

				if(MovedParticle.cellN != k){
					if(Config.headPS[MovedParticle.cellN] != randP ){
						int vvv = Config.headPS[MovedParticle.cellN];
						vv = Config.linkP[vvv];
						while(vv!=randP){ 
							vvv = vv;
							vv = Config.linkP[vvv];
						}
				
						Config.linkP[vvv] = Config.linkP[randP];
					}else{
						Config.headPS[MovedParticle.cellN] = Config.linkP[randP];
					}
		    
					Config.linkP[randP] = Config.headPS[k];
					Config.headPS[k] = randP;
					MovedParticle.cellN = k;
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
	MovedParticle=Config.part[v];
        if(v < Config.Nc){
		pear_mesh.UpdateTrans(id[0], MovedParticle.trans);

		NCell.resize(0);
		NPart.resize(0);

		k = MovedParticle.cell;

		vv = Config.headP[k];
    
		while(vv != -1 && !inside){
			newvv = vv;
			if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
			if(v < newvv && !Config.part[newvv].already){ 
				NPart.push_back(newvv);
				Config.part[newvv].already=true;

				R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
				R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
				R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
				if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
				else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
				if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
				else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
				if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
				else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
				overlapPear(); 
			}

			vv = Config.linkP[vv];
		}
		NCell.push_back(k);
		Config.usedCell[k]=true;

		if( MovedParticle.cellN != MovedParticle.cell && !inside){
			k = MovedParticle.cellN;

			vv = Config.headP[k];
		    
			while(vv != -1 && !inside){
				newvv = vv;
				if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
				if(v < newvv && !Config.part[newvv].already){
					NPart.push_back(newvv);
					Config.part[newvv].already=true;
					R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
					R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
					R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
					if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
					else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
					if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
					else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
					if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
					else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
					overlapPear(); 
				}
				vv = Config.linkP[vv];
			}
			NCell.push_back(k);
			Config.usedCell[k]=true;
		}

		if(!inside){

			dist_ori[0] = distN*MovedParticle.ori[0];
			dist_ori[1] = distN*MovedParticle.ori[1];
			dist_ori[2] = distN*MovedParticle.ori[2];

			dsx = MovedParticle.pos[0]-dist_ori[0];
			dsy = MovedParticle.pos[1]-dist_ori[1];
			dsz = MovedParticle.pos[2]-dist_ori[2];

			in_test = false;
			if(dsx < 0){
				sx[0] = Config.WP[0]-1;
				sx[1] = sx[0] -1;
				sx[2] = 0;
				in_test = true;
			}else{
				if(dsx > Config.l[0]){
					sx[0] = 0;
					sx[1] = Config.WP[0]-1;
					sx[2] = 1;
					in_test = true;
				}else{
					sx[0] = dsx*Config.wP[0];

					if( sx[0] >= Config.WP[0]-2 || sx[0] <= 1){ 
						in_test = true;

						sx[1] = sx[0]-1;
						if(sx[1] < 0 ) sx[1] = Config.WP[0]-1;

						sx[2] = sx[0]+1;
						if(sx[2]> Config.WP[0]-1) sx[2] = 0;
					}else{
						sx[1] = sx[0]-1;
						sx[2] = sx[0]+1;
					}
				}
			}
		
			

			if(dsy < 0){
				sy[0] = Config.WP[1]-1;
				sy[1] = sy[0] -1;
				sy[2] = 0;
				in_test = true;
			}else{
				if(dsy > Config.l[1]){
					sy[0] = 0;
					sy[1] = Config.WP[1]-1;
					sy[2] = 1;
					in_test = true;
				}else{
					sy[0] = dsy*Config.wP[1];
					if( sy[0] >= Config.WP[1]-2 || sy[0] <= 1){ 
						in_test = true;

						sy[1] = sy[0]-1;
						if(sy[1] < 0 ) sy[1] = Config.WP[1]-1;

						sy[2] = sy[0]+1;
						if(sy[2]> Config.WP[1]-1) sy[2] = 0;
					}else{
						sy[1] = sy[0]-1;
						sy[2] = sy[0]+1;
					}
				}
			}

			if(dsz < 0){
				sz[0] = Config.WP[2]-1;
				sz[1] = sz[0] -1;
				sz[2] = 0;
				in_test = true;
			}else{
				if(dsz > Config.l[2]){
					sz[0] = 0;
					sz[1] = Config.WP[2]-1;
					sz[2] = 1;
					in_test = true;
				}else{
					sz[0] = dsz*Config.wP[2];
					if( sz[0] >= Config.WP[2]-2 || sz[0] <= 1){ 
						in_test = true;

						sz[1] = sz[0]-1;
						if(sz[1] < 0 ) sz[1] = Config.WP[2]-1;

						sz[2] = sz[0]+1;
						if(sz[2]> Config.WP[2]-1) sz[2] = 0;
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
							k = sx[ix] + Config.WP[0]*(sy[iy]+Config.WP[1]*sz[iz]);
							if( Config.usedCell[k] ) continue;
			
							vv = Config.headP[k];

							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
								if(v < newvv && !Config.part[newvv].already){
									NPart.push_back(newvv);
									Config.part[newvv].already=true;

									R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
									R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
									R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
									if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
									else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
									if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
									else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
									if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
									else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
									overlapPear(); 
								}
								vv = Config.linkP[vv];
							}
							NCell.push_back(k);
							Config.usedCell[k]=true;
						}
					}
				}
			}else{
				for (int ix=0; ix < 3 && !inside; ix++){
					for (int iy=0; iy < 3 && !inside; iy++){
						for (int iz=0; iz < 3 && !inside; iz++){
							k = sx[ix] + Config.WP[0]*(sy[iy]+Config.WP[1]*sz[iz]);
							if( Config.usedCell[k] ) continue;
			
							vv = Config.headP[k];

							while(vv != -1 && !inside){
								newvv = vv;
								if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
								if(v < newvv && !Config.part[newvv].already){
									NPart.push_back(newvv);
									Config.part[newvv].already=true;

									R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
									R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
									R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
									overlapPear(); 
								}
								vv = Config.linkP[vv];
							}
							NCell.push_back(k);
							Config.usedCell[k]=true;
						}
					}
				}
			}
		}

		if(!inside){
			dsxN = MovedParticle.pos[0]+dist_ori[0];

			dsyN = MovedParticle.pos[1]+dist_ori[1];

			dszN = MovedParticle.pos[2]+dist_ori[2];

			in_test=false;
			if(dsxN < 0){
				sxN[0] = Config.WP[0]-1;
				in_test = true;
			}else{
				if(dsxN > Config.l[0]){
					sxN[0] = 0;
					in_test = true;
				}else{
					sxN[0] = dsxN*Config.wP[0];
					if( sxN[0] >= Config.WP[0]-2 || sxN[0] <= 1) in_test = true;
				}
			}

			if(dsyN < 0){
				syN[0] = Config.WP[1]-1;
				in_test = true;
			}else{
				if(dsyN > Config.l[1]){
					syN[0] = 0;
					in_test = true;
				}else{
					syN[0] = dsyN*Config.wP[1];
					if( syN[0] >= Config.WP[1]-2 || syN[0] <= 1) in_test = true;
				}
			}

			if(dszN < 0){
				szN[0] = Config.WP[2]-1;
				in_test = true;
			}else{
				if(dszN > Config.l[2]){
					szN[0] = 0;
					in_test = true;
				}else{
					szN[0] = dszN*Config.wP[2];
					if( szN[0] >= Config.WP[2]-2 || szN[0] <= 1) in_test = true;
				}
			}

			if(sx[0] != sxN[0] || sy[0] != syN[0] || sz[0] != szN[0]){

				sxN[1] = sxN[0]-1;
				if(sxN[1] < 0 ) sxN[1] = Config.WP[0]-1;

				sxN[2] = sxN[0]+1;
				if(sxN[2]> Config.WP[0]-1) sxN[2] = 0;

				syN[1] = syN[0]-1;
				if(syN[1] < 0 ) syN[1] = Config.WP[1]-1;

				syN[2] = syN[0]+1;
				if(syN[2]> Config.WP[1]-1) syN[2] = 0;
				
				szN[1] = szN[0]-1;
				if(szN[1] < 0 ) szN[1] = Config.WP[2]-1;

				szN[2] = szN[0]+1;
				if(szN[2]> Config.WP[2]-1) szN[2] = 0;

				if(in_test){
					for (int ix=0; ix < 3 && !inside; ix++){

						for (int iy=0; iy < 3 && !inside; iy++){

							for (int iz=0; iz < 3 && !inside; iz++){

								k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
								if( Config.usedCell[k] ) continue;
				
								vv = Config.headP[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
									if(v < newvv && !Config.part[newvv].already){
										NPart.push_back(newvv);
										Config.part[newvv].already=true;

										R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
										if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
										else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
										if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
										else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
										if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
										else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
										overlapPear(); 
									}
									vv = Config.linkP[vv];
								}
							}
						}
					}
				}else{
					for (int ix=0; ix < 3 && !inside; ix++){

						for (int iy=0; iy < 3 && !inside; iy++){

							for (int iz=0; iz < 3 && !inside; iz++){

								k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
								if( Config.usedCell[k] ) continue;
				
								vv = Config.headP[k];
						    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
									if(v < newvv && !Config.part[newvv].already){
										NPart.push_back(newvv);
										Config.part[newvv].already=true;

										R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0];
										R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1];
										R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2];
										overlapPear(); 
									}
									vv = Config.linkP[vv];
								}
							}
						}
					}

				}
			}
		}
		for( vv = 0; vv < NPart.size(); vv++) Config.part[NPart[vv]].already=false;
		for( vv = 0; vv < NCell.size(); vv++) Config.usedCell[NCell[vv]]=false;

	}else{

		sx[0] = MovedParticle.pos[0]*Config.wS[0];
		sy[0] = MovedParticle.pos[1]*Config.wS[1];
		sz[0] = MovedParticle.pos[2]*Config.wS[2];

		kN = sx[0] + Config.WS[0]*(sy[0]+Config.WS[1]*sz[0]);
	
		if(Config.headS[kN] != -1 && Config.headS[kN] != randP) inside = true;

		if(!inside ){
			NPart.resize(0);

			k = MovedParticle.cellN;

			vv = Config.headP[k];
	    
			while(vv != -1 && !inside){
				newvv = vv;
				if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
				if(!Config.part[newvv].already){ 
					NPart.push_back(newvv);
					Config.part[newvv].already=true;	

					R[0] = MovedParticle.pos[0] - Config.part[newvv].pos[0];
					R[1] = MovedParticle.pos[1] - Config.part[newvv].pos[1];
					R[2] = MovedParticle.pos[2] - Config.part[newvv].pos[2];
					if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
					else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
					if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
					else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
					if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
					else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
					overlapSpherePear(); 
				}

				vv = Config.linkP[vv];
			}

			if(!inside){
				in_test=false;
			
				sxN[0] = MovedParticle.pos[0]*Config.wP[0];
				if( sxN[0] >= Config.WP[0]-2 || sxN[0] <= 1){ 
					in_test = true;

					sxN[1] = sxN[0]-1;
					if(sxN[1] < 0 ) sxN[1] = Config.WP[0]-1;

					sxN[2] = sxN[0]+1;
					if(sxN[2]> Config.WP[0]-1) sxN[2] = 0;
				}else{
					sxN[1] = sxN[0]-1;

					sxN[2] = sxN[0]+1;
				}

				syN[0] = MovedParticle.pos[1]*Config.wP[1];
				if( syN[0] >= Config.WP[1]-2 || syN[0] <= 1){ 
					in_test = true;

					syN[1] = syN[0]-1;
					if(syN[1] < 0 ) syN[1] = Config.WP[1]-1;

					syN[2] = syN[0]+1;
					if(syN[2]> Config.WP[1]-1) syN[2] = 0;
				}else{
					syN[1] = syN[0]-1;

					syN[2] = syN[0]+1;
				}

				szN[0] = MovedParticle.pos[2]*Config.wP[2];
				if( szN[0] >= Config.WP[2]-2 || szN[0] <= 1){ 
					in_test = true;

					szN[1] = szN[0]-1;
					if(szN[1] < 0 ) szN[1] = Config.WP[2]-1;

					szN[2] = szN[0]+1;
					if(szN[2]> Config.WP[2]-1) szN[2] = 0;
				}else{
					szN[1] = szN[0]-1;

					szN[2] = szN[0]+1;
				}


				if(in_test){
					for (int ix=0; ix < 3 && !inside; ix++){

						for (int iy=0; iy < 3 && !inside; iy++){

							for (int iz=0; iz < 3 && !inside; iz++){
								k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
								if( k == MovedParticle.cellN ) continue;
					
								vv = Config.headP[k];
							    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
									if(!Config.part[newvv].already){ 
										NPart.push_back(newvv);
										Config.part[newvv].already=true;	

										R[0] = MovedParticle.pos[0] - Config.part[newvv].pos[0];
										R[1] = MovedParticle.pos[1] - Config.part[newvv].pos[1];
										R[2] = MovedParticle.pos[2] - Config.part[newvv].pos[2];
										if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
										else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];
										if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
										else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];
										if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
										else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];
										overlapSpherePear(); 
									}
									vv = Config.linkP[vv];
								}
							}
						}
					}
				}else{
					for (int ix=0; ix < 3 && !inside; ix++){

						for (int iy=0; iy < 3 && !inside; iy++){

							for (int iz=0; iz < 3 && !inside; iz++){
								k = sxN[ix] + Config.WP[0]*(syN[iy]+Config.WP[1]*szN[iz]);
								if( k == MovedParticle.cellN ) continue;
					
								vv = Config.headP[k];
							    
								while(vv != -1 && !inside){
									newvv = vv;
									if(newvv >= Config.part.size() ) newvv -= Config.part.size(); 
									if(!Config.part[newvv].already){ 
										NPart.push_back(newvv);
										Config.part[newvv].already=true;	

										R[0] = MovedParticle.pos[0] - Config.part[newvv].pos[0];
										R[1] = MovedParticle.pos[1] - Config.part[newvv].pos[1];
										R[2] = MovedParticle.pos[2] - Config.part[newvv].pos[2];
										overlapSpherePear(); 
									}
									vv = Config.linkP[vv];
								}
							}
						}
					}
				}
			}

			for( vv = 0; vv < NPart.size(); vv++) Config.part[NPart[vv]].already=false;
		}


		if(!inside){
			in_test=false;

			if( sx[0] >= Config.WS[0]-2 || sx[0] <= 1){ 
				in_test = true;
				sx[1] = sx[0]-1;
				if(sx[1] < 0 ){ 
					sx[1] = Config.WS[0]-1;
					bx[1]=-1;
					bx[3]=-1;
				}

				sx[2] = sx[0]+1;
				if(sx[2]> Config.WS[0]-1){ 
					sx[2] = 0;
					bx[2]=1;
					bx[4]=1;
				}

				sx[3] = sx[1]-1;
				if(sx[3] < 0 ){ 
					sx[3] = Config.WS[0]-1;
					bx[3]=-1;
				}

				sx[4] = sx[2]+1;
				if(sx[4]> Config.WS[0]-1){ 
					sx[4] = 0;
					bx[4]=1;
				}
			}else{
				sx[1] = sx[0]-1;
				sx[2] = sx[0]+1;
				sx[3] = sx[1]-1;
				sx[4] = sx[2]+1;
			}

			if( sy[0] >= Config.WS[1]-2 || sy[0] <= 1){ 
				in_test = true;

				sy[1] = sy[0]-1;
				if(sy[1] < 0 ){ 
					sy[1] = Config.WS[1]-1;
					by[1]=-1;
					by[3]=-1;
				}

				sy[2] = sy[0]+1;
				if(sy[2]> Config.WS[1]-1){ 
					sy[2] = 0;
					by[2]=1;
					by[4]=1;
				}

				sy[3] = sy[1]-1;
				if(sy[3] < 0 ){ 
					sy[3] = Config.WS[1]-1;
					by[3]=-1;
				}

				sy[4] = sy[2]+1;
				if(sy[4]> Config.WS[1]-1){ 
					sy[4] = 0;
					by[4]=1;
				}
			}else{
				sy[1] = sy[0]-1;
				sy[2] = sy[0]+1;
				sy[3] = sy[1]-1;
				sy[4] = sy[2]+1;
			}

			if( sz[0] >= Config.WS[2]-2 || sz[0] <= 1){ 
				in_test = true;

				sz[1] = sz[0]-1;
				if(sz[1] < 0 ){ 
					sz[1] = Config.WS[2]-1;
					bz[1]=-1;
					bz[3]=-1;
				}

				sz[2] = sz[0]+1;
				if(sz[2]> Config.WS[2]-1){ 
					sz[2] = 0;
					bz[2]=1;
					bz[4]=1;
				}

				sz[3] = sz[1]-1;
				if(sz[3] < 0 ){ 
					sz[3] = Config.WS[2]-1;
					bz[3]=-1;
				}

				sz[4] = sz[2]+1;
				if(sz[4]> Config.WS[2]-1){ 
					sz[4] = 0;
					bz[4]=1;
				}
			}else{
				sz[1] = sz[0]-1;
				sz[2] = sz[0]+1;
				sz[3] = sz[1]-1;
				sz[4] = sz[2]+1;
			}

			if(in_test){
				for (int ix=0; ix < 5 && !inside; ix++){

					for (int iy=0; iy < 5 && !inside; iy++){

						for (int iz=0; iz < 5 && !inside; iz++){
							k = sx[ix] + Config.WS[0]*sy[iy];
							if( k == Config.part[v].cell ) continue;
			
							vv = Config.headS[k];
					    
							lb[0]= bx[ix]*Config.l[0];
							lb[1]= by[iy]*Config.l[1];
							lb[2]= bz[iz]*Config.l[2];
						    
							if(vv != -1 && v < vv){
								R[0] = MovedParticle.pos[0] - Config.part[vv].pos[0] - lb[0];
								R[1] = MovedParticle.pos[1] - Config.part[vv].pos[1] - lb[1];
								R[2] = MovedParticle.pos[2] - Config.part[vv].pos[2] - lb[2];
								overlapSphere(); 
							}
						}
					}
				}
				bx[0] =0;
				bx[1] =0;
				bx[2] =0;
				bx[3] =0;
				bx[4] =0;

				by[0] =0;
				by[1] =0;
				by[2] =0;
				by[3] =0;
				by[4] =0;

				bz[0] =0;
				bz[1] =0;
				bz[2] =0;
				bz[3] =0;
				bz[4] =0;
			}else{
				for (int ix=0; ix < 5 && !inside; ix++){

					for (int iy=0; iy < 5 && !inside; iy++){

						for (int iz=0; iz < 5 && !inside; iz++){
							k = sx[ix] + Config.WS[0]*sy[iy];
							if( k == Config.part[v].cell ) continue;
			
							vv = Config.headS[k];
					    
							if(vv != -1 && v < vv){
								R[0] = MovedParticle.pos[0] - Config.part[vv].pos[0];
								R[1] = MovedParticle.pos[1] - Config.part[vv].pos[1];
								R[2] = MovedParticle.pos[2] - Config.part[vv].pos[2];
								overlapSphere(); 
							}
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
		std::cout << "(lx,ly) = (" << Config.l[0] << "," << Config.l[1] << ")"<< std::endl;

		Config.write(savefile,1);

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
	RenewList();

	if(Compressing){
		vproc = 0.5*(1+vproc);
		if(vproc > 0.9999){ 
		    vproc = 0.9999;
		    initCompress = true;
		}
		std::cout << "Compession unsuccessful! New vproc=" << vproc << std::endl;
		std::cout << "(lx,ly,lz) = (" << Config.l[0] << "," << Config.l[1] << ")"<< std::endl;
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
