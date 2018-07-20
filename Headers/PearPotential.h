void overlapSphere(){

	Rsq = R[0]*R[0] + R[1]*R[1]+R[2]*R[2];

	if(Rsq < rcut_SPH ) inside=true;
}

void overlapSpherePear (){
	
	Rsq = R[0]*R[0] + R[1]*R[1]+R[2]*R[2];

    	if(Rsq < rcut_PSPH){
        
        	rw = part[newvv].ori[0]*R[0] + part[newvv].ori[1]*R[1] + part[newvv].ori[2]*R[2];
       
		rwT = Rsq - rw*rw;

		if( rwT < rcut_PSPH_T ){
			if( rwT < 1e-4 ) inside = true;
			else{
				rwT = sqrt(rwT);
                            	wT[0] = (R[0] - rw*part[newvv].ori[0])/rwT;
                            	wT[1] = (R[1] - rw*part[newvv].ori[1])/rwT;
                            	wT[2] = (R[2] - rw*part[newvv].ori[2])/rwT;

				in_test = false;

				for( int i=0; i < bezier_x.size(); i++){
					if(rw < bezier_z[i] && !in_test){
						in_test=true;
						if(rwT < bezier_x[i]){ 
							inside=true;
							break;
						}
					}

					surf_point[0] = bezier_z[i]*part[newvv].ori[0] + bezier_x[i]*wT[0] - R[0];
					surf_point[1] = bezier_z[i]*part[newvv].ori[1] + bezier_x[i]*wT[1] - R[1];
					surf_point[2] = bezier_z[i]*part[newvv].ori[2] + bezier_x[i]*wT[2] - R[2];
					Rsq = surf_point[0]*surf_point[0] + surf_point[1]*surf_point[1] + surf_point[2]*surf_point[2];

					if(Rsq <  rsphere*rsphere){
						inside=true;
						break;
					}
				}
			}
		}
    	}
}

void overlapPearSphere (){
	
	Rsq = R[0]*R[0] + R[1]*R[1]+R[2]*R[2];

    	if(Rsq < rcut_PSPH){
        
        	rw = MovedParticle.ori[0]*R[0] + MovedParticle.ori[1]*R[1] + MovedParticle.ori[2]*R[2];
       
		rwT = Rsq - rw*rw;

		if( rwT < rcut_PSPH_T ){
			if( rwT < 1e-4 ) inside = true;
			else{
				rwT = sqrt(rwT);
                            	wT[0] = (R[0] - rw*MovedParticle.ori[0])/rwT;
                            	wT[1] = (R[1] - rw*MovedParticle.ori[1])/rwT;
                            	wT[2] = (R[2] - rw*MovedParticle.ori[2])/rwT;
    
				in_test = false;

				for( int i=0; i < bezier_x.size(); i++){

					if(rw < bezier_z[i] && !in_test){
						in_test=true;
						if(rwT < bezier_x[i]){ 
							inside=true;
							break;
						}
					}

					surf_point[0] = bezier_z[i]*MovedParticle.ori[0] + bezier_x[i]*wT[0] - R[0];
					surf_point[1] = bezier_z[i]*MovedParticle.ori[1] + bezier_x[i]*wT[1] - R[1];
					surf_point[2] = bezier_z[i]*MovedParticle.ori[2] + bezier_x[i]*wT[2] - R[2];

					Rsq = surf_point[0]*surf_point[0] + surf_point[1]*surf_point[1] + surf_point[2]*surf_point[2];

					if(Rsq <  rsphere*rsphere){
						inside=true;
						break;
					}
				}
			}
		}
    	}
}

void overlapPear (){
	
	Rsq = R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
    
	if(Rsq < rcut_P){
        	rw = MovedParticle.ori[0]*R[0] + MovedParticle.ori[1]*R[1] + MovedParticle.ori[2]*R[2];
       
		rwT = Rsq - rw*rw;
		if( rwT < rcut_P_T ){
			rw1 = part[newvv].ori[0]*R[0] + part[newvv].ori[1]*R[1] + part[newvv].ori[2]*R[2];
	       
			rwT = Rsq - rw1*rw1;

			if( rwT < rcut_P_T ){

				ww = MovedParticle.ori[0]*part[newvv].ori[0] + MovedParticle.ori[1]*part[newvv].ori[1] + MovedParticle.ori[2]*part[newvv].ori[2];
				www = 1 - ww*ww;
				if( www < 1e-4){ 
					if(rwT < rcut_P_II){
						part[newvv].trans[0][3] = R[0];	
						part[newvv].trans[1][3] = R[1];	
						part[newvv].trans[2][3] = R[2];	


						pear_mesh.UpdateTrans(id[1], part[newvv].trans);

						VCReport report;

						pear_mesh.Collide( &report );

						part[newvv].trans[0][3] = 0;	
						part[newvv].trans[1][3] = 0;	
						part[newvv].trans[2][3] = 0;	

						if(report.numObjPairs() > 0) inside = true;
					}
				}else{
					xlambda = (rw-ww*rw1)/www;
					xmu = xlambda*ww-rw1;

					if( fabs(xlambda) > max_z || fabs(xmu) > max_z ){

						if(fabs(xlambda) > fabs(xmu) ){
							if(xlambda < 0 ) xlambda = -max_z;
							else xlambda = max_z;
							
							xmu = xlambda*ww-rw1;
							if(xmu < -max_z) xmu = -max_z;
							else if(xmu > max_z) xmu = max_z; 
						}else{
							if(xmu < 0 ) xmu = -max_z;
							else xmu = max_z;
							
							xlambda = xmu*ww+rw;
							if(xlambda < -max_z) xlambda = -max_z;
							else if(xlambda > max_z) xlambda = max_z; 
						}

					}
/*
					xlambda = rw-ww*rw1;
					xmu = -rw1+ww*rw;

					wwww = www*www;
					rwT = Rsq*wwww + xlambda*xlambda + xmu*xmu - 2*xlambda*xmu*ww + 2*xmu*rw1*www-2*xlambda*rw*www;

					if(rwT < rcut_P_II*wwww){
*/
					
					rwT = Rsq + xlambda*xlambda +xmu*xmu + 2*(xmu*rw1-xlambda*(rw+xmu*ww));
					if(rwT < rcut_P_II){
						part[newvv].trans[0][3] = R[0];	
						part[newvv].trans[1][3] = R[1];	
						part[newvv].trans[2][3] = R[2];	


						pear_mesh.UpdateTrans(id[1], part[newvv].trans);

						VCReport report;

						pear_mesh.Collide( &report );

						part[newvv].trans[0][3] = 0;	
						part[newvv].trans[1][3] = 0;	
						part[newvv].trans[2][3] = 0;	

						if(report.numObjPairs() > 0) inside = true;
					}
				}
			}
		}
	}
}


bool overlapSPH ( class teilchen sphere1, class teilchen sphere2, std::vector<double> l){

    std::vector<double> r = DeltaR(sphere1.pos,sphere2.pos,l);
	double rsq = scal_p(r,r);

    if(rsq > rcut_SPH ) return false;
    else return true;
}




bool overlapPSPH ( class teilchen pear, class teilchen sphere, std::vector<double> l){
	
	std::vector<double> r (3);
	double r2;

	r[0] = sphere.pos[0] - pear.pos[0]; // boundary correction
	if (r[0] > l[0]*0.5)  r[0] -= l[0];
	else if (r[0] < -l[0]*0.5)  r[0] += l[0];

	r[1] = sphere.pos[1] - pear.pos[1]; // boundary correction
	if (r[1] > l[1]*0.5)  r[1] -= l[1];
	else if (r[1] < -l[1]*0.5)  r[1] += l[1];


	r[2] = sphere.pos[2] - pear.pos[2]; // boundary correction
	if (r[2] > l[2]*0.5)  r[2] -= l[2];
	else if (r[2] < -l[2]*0.5)  r[2] += l[2];

 	r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    

    	if(r2 > rcut_PSPH) return false;
    	else{
        
        	double rw = pear.ori[0]*r[0] + pear.ori[1]*r[1] + pear.ori[2]*r[2];
       
		double rwT = r2 - rw*rw;

		if( rwT > rcut_PSPH_T ) return false;
		else{
			if( rwT < 1e-4 ) return true;
			else{
				rwT = sqrt(rwT);
				std::vector<double> wT (3);
                            	wT[0] = (r[0] - rw*pear.ori[0])/rwT;
                            	wT[1] = (r[1] - rw*pear.ori[1])/rwT;
                            	wT[2] = (r[2] - rw*pear.ori[2])/rwT;

                            	//double min_val = 10000;

				bool in_test = false;
				bool wrong_side = false;

                            	for( int i=0; i < bezier_x.size(); i++){

					if(rw < bezier_z[i] && !in_test){
						in_test=true;
						if(rwT < bezier_x[i]){ 
							wrong_side=true;
							break;
						}
					}

                                	std::vector<double> surf_point (3);
                                	surf_point[0] = bezier_z[i]*pear.ori[0] + bezier_x[i]*wT[0] - r[0];
                                	surf_point[1] = bezier_z[i]*pear.ori[1] + bezier_x[i]*wT[1] - r[1];
                                	surf_point[2] = bezier_z[i]*pear.ori[2] + bezier_x[i]*wT[2] - r[2];

                                	r2 = surf_point[0]*surf_point[0] + surf_point[1]*surf_point[1]+ surf_point[2]*surf_point[2];

                                	if(r2 <  rsphere*rsphere ){ 
						wrong_side = true;
						break;
					}


                                	//if(min_val > r2) min_val = r2;
                                	//else if(in_test) break;
                            	}	
                            	return wrong_side;
			}

		}
    	}
}

bool overlapP ( class teilchen pear1, class teilchen pear2, std::vector<double> l){
	
	std::vector<double> r (3);
	double r2;
	
	r[0] = pear2.pos[0] - pear1.pos[0]; // boundary correction
	if (r[0] > l[0]*0.5){  
		r[0] -= l[0];
	}else{ 
		if (r[0] < -l[0]*0.5){  
			r[0] += l[0];
		}
	}

	r[1] = pear2.pos[1] - pear1.pos[1]; // boundary correction
	if (r[1] > l[1]*0.5){  
		r[1] -= l[1];
	}else{ 
		if (r[1] < -l[1]*0.5){ 
			r[1] += l[1];
		}
	}

	r[2] = pear2.pos[2] - pear1.pos[2]; // boundary correction
	if (r[2] > l[2]*0.5){  
		r[2] -= l[2];
	}else{ 
		if (r[2] < -l[2]*0.5){ 
			r[2] += l[2];
		}
	}

	r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    
    	if(r2 > rcut_P) return false;
    	else{

		pear2.trans[0][3] = r[0];	
		pear2.trans[1][3] = r[1];	
		pear2.trans[2][3] = r[2];	
		pear_mesh.UpdateTrans(id[1], pear2.trans);
		pear_mesh.UpdateTrans(id[0], pear1.trans);

		VCReport report;

		pear_mesh.Collide( &report );
		
		pear2.trans[0][3] = 0;	
		pear2.trans[1][3] = 0;	
		pear2.trans[2][3] = 0;	

		if(report.numObjPairs() > 0) return true;
		else return false;
    	}
}
