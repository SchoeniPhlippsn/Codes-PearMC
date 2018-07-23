void RenewList(){ // Initialise Neighbour List

	WP[0] = (int)(l[0]/rlistP); 
	WP[1] = (int)(l[1]/rlistP); 
	WP[2] = (int)(l[2]/rlistP); 

	std::cout << WP[0] << "|" << WP[1] << "|" << WP[2] << std::endl;

	wP[0] = WP[0]/l[0];
	wP[1] = WP[1]/l[1];
	wP[2] = WP[2]/l[1];
	
	headP.resize(WP[0]*WP[1]*WP[2]); 
	usedCell.resize(WP[0]*WP[1]*WP[2],false); 
	for( int i=0; i<headP.size(); i++) headP[i] = -1;

	linkP.resize(Nc+part.size()); 
	for( int i=0; i<linkP.size(); i++) linkP[i] = -1;

	s_n.resize(WP[0]);

	s_nn[0] = 0;
	s_nn[1] = 1;
	s_nn[2] = WP[0]-1;
	s_nn[3] = 2;
	s_nn[4] = WP[0]-2;
	s_n[0]=s_nn;

	s_nn[0] = 1;
	s_nn[1] = 2;
	s_nn[2] = 0;
	s_nn[3] = 3;
	s_nn[4] = WP[0]-1;
	s_n[1]=s_nn;
	for( int i=2; i<s_n.size()-2; i++){ 
		s_nn[0] = i;
		s_nn[1] = i+1;
		s_nn[3] = i+2;
		s_nn[2] = i-1;
		s_nn[4] = i-2;

		s_n[i]=s_nn;
	}

	s_nn[0] = WP[0]-2;
	s_nn[1] = WP[0]-1;
	s_nn[2] = WP[0]-3;
	s_nn[3] = 0;
	s_nn[4] = WP[0]-4;

	s_n[WP[0]-2]=s_nn;

	s_nn[0] = WP[0]-1;
	s_nn[1] = 0;
	s_nn[2] = WP[0]-2;
	s_nn[3] = 1;
	s_nn[4] = WP[0]-3;

	s_n[WP[0]-1]=s_nn;


	for( int i=0; i<Nc; i++){
		dsx = part[i].pos[0]-part[i].dist_ori[0];
		if(dsx < 0) part[i].s[0]= WP[0]-1;
		else{
			part[i].s[0] = dsx*wP[0];
			if(part[i].s[0]>WP[0]-1) part[i].s[0] = 0; 
		}

		dsy = part[i].pos[1]-part[i].dist_ori[1];
		if(dsy < 0) part[i].s[1]= WP[1]-1;
		else{
			part[i].s[1] = dsy*wP[1];
			if(part[i].s[1]>WP[1]-1) part[i].s[1] = 0; 
		}

		dsz = part[i].pos[2]-part[i].dist_ori[2];
		if(dsz < 0) part[i].s[2]= WP[2]-1;
		else{
			part[i].s[2] = dsz*wP[2];
			if(part[i].s[2]>WP[2]-1) part[i].s[2] = 0; 
		}

		k = part[i].s[0] + WP[0]*(part[i].s[1]+WP[1]*part[i].s[2]); 

		part[i].cell = k;
		if(headP[k]==-1) headP[k] = i;
		else{
		   linkP[i] = headP[k];
		   headP[k] = i; 
		}

		dsx = part[i].pos[0]+part[i].dist_ori[0];
		if(dsx < 0) part[i].sN[0]= WP[0]-1;
		else{
			part[i].sN[0] = dsx*wP[0];
			if(part[i].sN[0]>WP[0]-1) part[i].sN[0] = 0;
		}

		dsy = part[i].pos[1]+part[i].dist_ori[1];
		if(dsy < 0) part[i].sN[1]= WP[1]-1;
		else{
			part[i].sN[1] = dsy*wP[1];
			if(part[i].sN[1]>WP[1]-1) part[i].sN[1] = 0; 
		}

		dsz = part[i].pos[2]+part[i].dist_ori[2];
		if(dsz < 0) part[i].sN[2]= WP[2]-1;
		else{
			part[i].sN[2] = dsz*wP[2];
			if(part[i].sN[2]>WP[2]-1) part[i].sN[2] = 0; 
		}

		k = part[i].sN[0] + WP[0]*(part[i].sN[1]+WP[1]*part[i].sN[2]); 

		part[i].cellN = k;
		i += part.size();
		if(headP[k]==-1) headP[k] = i;
		else{
		   linkP[i] = headP[k];
		   headP[k] = i; 
		}
		i -= part.size();
	}
};
