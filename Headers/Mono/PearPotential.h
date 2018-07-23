void overlapPear (){
	
	Rsq = R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
    
	if(Rsq < rcut_P){
        	rw = MovedParticle.trans[0][2]*R[0] + MovedParticle.trans[1][2]*R[1] + MovedParticle.trans[2][2]*R[2];
       
		rwT = Rsq - rw*rw;
		if( rwT < rcut_P_T ){
			rw1 = part[newvv].trans[0][2]*R[0] + part[newvv].trans[1][2]*R[1] + part[newvv].trans[2][2]*R[2];
	       
			rwT = Rsq - rw1*rw1;

			if( rwT < rcut_P_T ){

				ww = MovedParticle.trans[0][2]*part[newvv].trans[0][2] + MovedParticle.trans[1][2]*part[newvv].trans[1][2] + MovedParticle.trans[2][2]*part[newvv].trans[2][2];
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
