void overlapPear (){
	
	Rsq = Rx*Rx+Ry*Ry+Rz*Rz;
    
	if(Rsq < rcut_P){
        	rw = MovedParticle.trans[0][2]*Rx + MovedParticle.trans[1][2]*Ry + MovedParticle.trans[2][2]*Rz;
       
		rwT = Rsq - rw*rw;
		if( rwT < rcut_P_T ){
			rw1 = part[newvv].trans[0][2]*Rx + part[newvv].trans[1][2]*Ry + part[newvv].trans[2][2]*Rz;
	       
			rwT = Rsq - rw1*rw1;

			if( rwT < rcut_P_T ){

				ww = MovedParticle.trans[0][2]*part[newvv].trans[0][2] + MovedParticle.trans[1][2]*part[newvv].trans[1][2] + MovedParticle.trans[2][2]*part[newvv].trans[2][2];
				www = 1 - ww*ww;
				if( www < 1e-4){ 
					if(rwT < rcut_P_II){
						part[newvv].trans[0][3] = Rx;	
						part[newvv].trans[1][3] = Ry;	
						part[newvv].trans[2][3] = Rz;	


						pear_mesh.UpdateTrans(id[1], part[newvv].trans);

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
						part[newvv].trans[0][3] = Rx;	
						part[newvv].trans[1][3] = Ry;	
						part[newvv].trans[2][3] = Rz;	


						pear_mesh.UpdateTrans(id[1], part[newvv].trans);

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

bool overlapP ( double rx, double ry,  double rz, class teilchen pear1, class teilchen pear2){
	
	double r2;
	r2 = rx*rx + ry*ry + rz*rz;
    
    	if(r2 > rcut_P) return false;
    	else{

		pear2.trans[0][3] = rx;	
		pear2.trans[1][3] = ry;	
		pear2.trans[2][3] = rz;	
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
