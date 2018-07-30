void RenewList(){ // Initialise Neighbour List

	WPx = (int)(lx/rlistP); 
	WPy = (int)(ly/rlistP); 
	WPz = (int)(lz/rlistP); 

	std::cout << WPx << "|" << WPy << "|" << WPz << std::endl;

	wPx = WPx/lx;
	wPy = WPy/ly;
	wPz = WPz/lz;
	
	headP.resize(WPx*WPy*WPz); 
	usedCell.resize(WPx*WPy*WPz,false); 
	for( int i=0; i<headP.size(); i++) headP[i] = -1;

	linkP.resize(2*Nc); 
	for( int i=0; i<linkP.size(); i++) linkP[i] = -1;

	s_n.resize(WPx);

	s_nn[0] = 0;
	s_nn[1] = 1;
	s_nn[2] = WPx-1;
	s_nn[3] = 2;
	s_nn[4] = WPx-2;
	s_n[0]=s_nn;

	s_nn[0] = 1;
	s_nn[1] = 2;
	s_nn[2] = 0;
	s_nn[3] = 3;
	s_nn[4] = WPx-1;
	s_n[1]=s_nn;
	for( int i=2; i<s_n.size()-2; i++){ 
		s_nn[0] = i;
		s_nn[1] = i+1;
		s_nn[3] = i+2;
		s_nn[2] = i-1;
		s_nn[4] = i-2;

		s_n[i]=s_nn;
	}

	s_nn[0] = WPx-2;
	s_nn[1] = WPx-1;
	s_nn[2] = WPx-3;
	s_nn[3] = 0;
	s_nn[4] = WPx-4;

	s_n[WPx-2]=s_nn;

	s_nn[0] = WPx-1;
	s_nn[1] = 0;
	s_nn[2] = WPx-2;
	s_nn[3] = 1;
	s_nn[4] = WPx-3;

	s_n[WPx-1]=s_nn;


	for( int i=0; i<Nc; i++){

		dsx = posx[i]-dist_orix[i];
		if(dsx < 0) sx[i]= WPx-1;
		else{
			sx[i] = dsx*wPx;
			if(sx[i]>WPx-1) sx[i] = 0; 
		}

		dsy = posy[i]-dist_oriy[i];
		if(dsy < 0) sy[i]= WPy-1;
		else{
			sy[i] = dsy*wPy;
			if(sy[i]>WPy-1) sy[i] = 0; 
		}

		dsz = posz[i]-dist_oriz[i];
		if(dsz < 0) sz[i]= WPz-1;
		else{
			sz[i] = dsz*wPz;
			if(sz[i]>WPz-1) sz[i] = 0; 
		}

		k = sx[i] + WPx*(sy[i]+WPy*sz[i]); 

		cell[i] = k;
		if(headP[k]==-1) headP[k] = i;
		else{
		   linkP[i] = headP[k];
		   headP[k] = i; 
		}

		dsx = posx[i]+dist_orix[i];
		if(dsx < 0) sNx[i] = WPx-1;
		else{
			sNx[i] = dsx*wPx;
			if(sNx[i]>WPx-1) sNx[i] = 0;
		}

		dsy = posy[i]+dist_oriy[i];
		if(dsy < 0) sNy[i] = WPy-1;
		else{
			sNy[i] = dsy*wPy;
			if(sNy[i]>WPy-1) sNy[i] = 0;
		}

		dsz = posz[i]+dist_oriz[i];
		if(dsz < 0) sNz[i] = WPz-1;
		else{
			sNz[i] = dsz*wPz;
			if(sNz[i]>WPz-1) sNz[i] = 0;
		}

		k = sNx[i] + WPx*(sNy[i]+WPy*sNz[i]); 

		cellN[i] = k;
		i+=Nc;
		if(headP[k]==-1) headP[k] = i;
		else{
		   linkP[i] = headP[k];
		   headP[k] = i; 
		}
		i-=Nc;
	}
};
