void overlapSphere(){

	std::vector<double> R (3);
	double Rsq;

	R[0] = MovedParticle.pos[0] - Config.part[vv].pos[0]; // boundary correction
	if (R[0] > Config.l_2[0])  R[0] -= Config.l[0];
	else if (R[0] < -Config.l_2[0])  R[0] += Config.l[0];

	R[1] = MovedParticle.pos[1] - Config.part[vv].pos[1]; // boundary correction
	if (R[1] > Config.l_2[1])  R[1] -= Config.l[1];
	else if (R[1] < -Config.l_2[1])  R[1] += Config.l[1];

	R[2] = MovedParticle.pos[2] - Config.part[vv].pos[2]; // boundary correction
	if (R[2] > Config.l_2[2])  R[2] -= Config.l[2];
	else if (R[2] < -Config.l_2[2])  R[2] += Config.l[2];

	Rsq = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

	if(Rsq < rcut_SPH ) inside=true;
}

void overlapPear (){
	
	std::vector<double> R (3);
	NPart.push_back(newvv);
	Config.part[newvv].already=true;

	R[0] = Config.part[newvv].pos[0] - MovedParticle.pos[0]; // boundary correction
	if (R[0] > Config.l_2[0]) R[0] -= Config.l[0];
	else if (R[0] < -Config.l_2[0]) R[0] += Config.l[0];

	R[1] = Config.part[newvv].pos[1] - MovedParticle.pos[1]; // boundary correction
	if (R[1] > Config.l_2[1]) R[1] -= Config.l[1];
	else if (R[1] < -Config.l_2[1]) R[1] += Config.l[1];

	R[2] = Config.part[newvv].pos[2] - MovedParticle.pos[2]; // boundary correction
	if (R[2] > Config.l_2[2]) R[2] -= Config.l[2];
	else if (R[2] < -Config.l_2[2])	R[2] += Config.l[2];

	double Rsq = R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
    
	if(Rsq < rcut_P){

		Config.part[newvv].trans[0][3] = R[0];	
		Config.part[newvv].trans[1][3] = R[1];	
		Config.part[newvv].trans[2][3] = R[2];	


		pear_mesh.UpdateTrans(id[1], Config.part[newvv].trans);

		VCReport report;

		pear_mesh.Collide( &report );

		Config.part[newvv].trans[0][3] = 0;	
		Config.part[newvv].trans[1][3] = 0;	
		Config.part[newvv].trans[2][3] = 0;	

		if(report.numObjPairs() > 0) inside = true;

	}
}


bool overlapSPH ( class teilchen sphere1, class teilchen sphere2, std::vector<double> l){

    std::vector<double> R = DeltaR(sphere1.pos,sphere2.pos,l);
	double Rsq = scal_p(R,R);

    if(Rsq > rcut_SPH ) return false;
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
				rwT = 1/sqrt(rwT);
				std::vector<double> wT (3);
                            	wT[0] = (r[0] - rw*pear.ori[0])*rwT;
                            	wT[1] = (r[1] - rw*pear.ori[1])*rwT;
                            	wT[2] = (r[2] - rw*pear.ori[2])*rwT;

                            	double min_val = 10000;

                            	for( int i=0; i < bezier_x.size(); i++){
                                	std::vector<double> surf_point (3);
                                	surf_point[0] = bezier_z[i]*pear.ori[0] + bezier_x[i]*wT[0] - r[0];
                                	surf_point[1] = bezier_z[i]*pear.ori[1] + bezier_x[i]*wT[1] - r[1];
                                	surf_point[2] = bezier_z[i]*pear.ori[2] + bezier_x[i]*wT[2] - r[2];

                                	r2 = surf_point[0]*surf_point[0] + surf_point[1]*surf_point[1] +surf_point[2]*surf_point[2];

                                	if(r2 <  rsphere*rsphere ) break;

                                	if(min_val > r2) min_val = r2;
                                	else break;
                            	}
                            	if(r2 <  rsphere*rsphere ) return true;
				else return false;
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
