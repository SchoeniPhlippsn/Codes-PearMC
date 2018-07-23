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

	for( int i=0; i<Nc; i++){
		dsx = part[i].pos[0]-part[i].dist_ori[0];
		if(dsx < 0) sx[1]= WP[0]-1;
		else{
			sx[1] = dsx*wP[0];
			if(sx[1]>WP[0]-1) sx[1] = 0; 
		}

		dsy = part[i].pos[1]-part[i].dist_ori[1];
		if(dsy < 0) sy[1]= WP[1]-1;
		else{
			sy[1] = dsy*wP[1];
			if(sy[1]>WP[1]-1) sy[1] = 0; 
		}

		dsz = part[i].pos[2]-part[i].dist_ori[2];
		if(dsz < 0) sz[1]= WP[2]-1;
		else{
			sz[1] = dsz*wP[2];
			if(sz[1]>WP[2]-1) sz[1] = 0; 
		}

		k = sx[1] + WP[0]*(sy[1]+WP[1]*sz[1]); 

		part[i].cell = k;
		if(headP[k]==-1) headP[k] = i;
		else{
		   linkP[i] = headP[k];
		   headP[k] = i; 
		}

		dsx = part[i].pos[0]+part[i].dist_ori[0];
		if(dsx < 0) sxN[1]= WP[0]-1;
		else{
			sxN[1] = dsx*wP[0];
			if(sxN[1]>WP[0]-1) sxN[1] = 0;
		}

		dsy = part[i].pos[1]+part[i].dist_ori[1];
		if(dsy < 0) syN[1]= WP[1]-1;
		else{
			syN[1] = dsy*wP[1];
			if(syN[1]>WP[1]-1) syN[1] = 0; 
		}

		dsz = part[i].pos[2]+part[i].dist_ori[2];
		if(dsz < 0) szN[1]= WP[2]-1;
		else{
			szN[1] = dsz*wP[2];
			if(szN[1]>WP[2]-1) szN[1] = 0; 
		}

		k = sxN[1] + WP[0]*(syN[1]+WP[1]*szN[1]); 

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
