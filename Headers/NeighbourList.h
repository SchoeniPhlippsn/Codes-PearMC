void RenewList(){ // Initialise Neighbour List

	WP[0] = (int)(l[0]/rlistP); 
	WP[1] = (int)(l[1]/rlistP); 
	WP[2] = (int)(l[2]/rlistP); 

	wP[0] = WP[0]/l[0];
	wP[1] = WP[1]/l[1];
	wP[2] = WP[2]/l[1];
	
	headP.resize(WP[0]*WP[1]*WP[2]); 
	usedCell.resize(WP[0]*WP[1]*WP[2],false); 
	for( int i=0; i<headP.size(); i++) headP[i] = -1;

	linkP.resize(Nc+part.size()); 
	for( int i=0; i<linkP.size(); i++) linkP[i] = -1;

	for( int i=0; i<Nc; i++){
		dsx = part[i].pos[0]-distN*part[i].ori[0];
		if(dsx < 0) sx[1]= WP[0]-1;
		else{
			sx[1] = dsx*wP[0];
			if(sx[1]>WP[0]-1) sx[1] = 0; 
		}

		dsy = part[i].pos[1]-distN*part[i].ori[1];
		if(dsy < 0) sy[1]= WP[1]-1;
		else{
			sy[1] = dsy*wP[1];
			if(sy[1]>WP[1]-1) sy[1] = 0; 
		}

		dsz = part[i].pos[2]-distN*part[i].ori[2];
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

		dsx = part[i].pos[0]+distN*part[i].ori[0];
		if(dsx < 0) sxN[1]= WP[0]-1;
		else{
			sxN[1] = dsx*wP[0];
			if(sxN[1]>WP[0]-1) sxN[1] = 0;
		}

		dsy = part[i].pos[1]+distN*part[i].ori[1];
		if(dsy < 0) syN[1]= WP[1]-1;
		else{
			syN[1] = dsy*wP[1];
			if(syN[1]>WP[1]-1) syN[1] = 0; 
		}

		dsz = part[i].pos[2]+distN*part[i].ori[2];
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

	if(Ns!=0){
		WS[0] = (int)(l[0]/rlistS); 
		WS[1] = (int)(l[1]/rlistS); 
		WS[2] = (int)(l[2]/rlistS); 
		
		wS[0] = WS[0]/l[0];
		wS[1] = WS[1]/l[1];
		wS[2] = WS[2]/l[2];
		
		headS.resize(WS[0]*WS[1]*WS[2]); 
		for( int i=0; i<headS.size(); i++) headS[i] = -1;

		for( int i=Nc; i<part.size(); i++){ 
		    
			if( part[i].pos[0] < 0 ) sx[1] = 0;
			else{
				sx[1] = part[i].pos[0]*wS[0];
				if(sx[1]>WS[0]-1) sx[1] = WS[0]-1;
			}
			if( part[i].pos[1] < 0 ) sy[1] = 0;
			else{
				sy[1] = part[i].pos[1]*wS[1]; 
				if(sy[1]>WS[1]-1) sy[1] = WS[1]-1;
			}
			if( part[i].pos[2] < 0 ) sz[1] = 0;
			else{
				sz[1] = part[i].pos[2]*wS[2];
				if(sz[1]>WS[2]-1) sz[1] = WS[2]-1;
			}
		   
			k = sx[1] + WS[0]*(sy[1]+WS[1]*sz[1]); 

			part[i].cell = k;
			headS[k] = i;

		}

		headPS.resize(WP[0]*WP[1]*WP[2]); 
		for( int i=0; i<headPS.size(); i++) headPS[i] = -1;

		for( int i=Nc; i<part.size(); i++){ 
		    
			if( part[i].pos[0] < 0 ) sx[1] = 0;
			else{
				sx[1] = part[i].pos[0]*wP[0];
				if(sx[1]>WP[0]-1) sx[1] = WP[0]-1; 
			}
			if( part[i].pos[1] < 0 ) sy[1] = 0;
			else{
				sy[1] = part[i].pos[1]*wP[1]; 
				if(sy[1]>WP[1]-1) sy[1] = WP[1]-1;
			}
			if( part[i].pos[2] < 0 ) sz[1] = 0;
			else{
				sz[1] = part[i].pos[2]*wP[2];
				if(sz[1]>WP[2]-1) sz[1] = WP[2]-1; 
			}
		   
			k = sx[1] + WP[0]*(sy[1]+WP[1]*sz[1]); 

			part[i].cellN = k;
			if(headPS[k]==-1) headPS[k] = i;
			else{
			   linkP[i] = headPS[k];
			   headPS[k] = i; 
			}
		}
	} 
};
