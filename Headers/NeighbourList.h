void RenewList(){ // Initialise Neighbour List

	Config.WP[0] = (int)(Config.l[0]/rlistP); 
	Config.WP[1] = (int)(Config.l[1]/rlistP); 
	Config.WP[2] = (int)(Config.l[2]/rlistP); 

	Config.wP[0] = Config.WP[0]/Config.l[0];
	Config.wP[1] = Config.WP[1]/Config.l[1];
	Config.wP[2] = Config.WP[2]/Config.l[1];
	
	Config.headP.resize(Config.WP[0]*Config.WP[1]*Config.WP[2]); 
	Config.usedCell.resize(Config.WP[0]*Config.WP[1]*Config.WP[2],false); 
	for( int i=0; i<Config.headP.size(); i++) Config.headP[i] = -1;

	Config.linkP.resize(Config.Nc+Config.part.size()); 
	for( int i=0; i<Config.linkP.size(); i++) Config.linkP[i] = -1;

	for( int i=0; i<Config.Nc; i++){
		dsx = Config.part[i].pos[0]-distN*Config.part[i].ori[0];
		if(dsx < 0) sx[1]= Config.WP[0]-1;
		else{
			sx[1] = dsx*Config.wP[0];
			if(sx[1]>Config.WP[0]-1) sx[1] = 0; 
		}

		dsy = Config.part[i].pos[1]-distN*Config.part[i].ori[1];
		if(dsy < 0) sy[1]= Config.WP[1]-1;
		else{
			sy[1] = dsy*Config.wP[1];
			if(sy[1]>Config.WP[1]-1) sy[1] = 0; 
		}

		dsz = Config.part[i].pos[2]-distN*Config.part[i].ori[2];
		if(dsz < 0) sz[1]= Config.WP[2]-1;
		else{
			sz[1] = dsz*Config.wP[2];
			if(sz[1]>Config.WP[2]-1) sz[1] = 0; 
		}

		k = sx[1] + Config.WP[0]*(sy[1]+Config.WP[1]*sz[1]); 

		Config.part[i].cell = k;
		if(Config.headP[k]==-1) Config.headP[k] = i;
		else{
		   Config.linkP[i] = Config.headP[k];
		   Config.headP[k] = i; 
		}

		dsx = Config.part[i].pos[0]+distN*Config.part[i].ori[0];
		if(dsx < 0) sxN[1]= Config.WP[0]-1;
		else{
			sxN[1] = dsx*Config.wP[0];
			if(sxN[1]>Config.WP[0]-1) sxN[1] = 0;
		}

		dsy = Config.part[i].pos[1]+distN*Config.part[i].ori[1];
		if(dsy < 0) syN[1]= Config.WP[1]-1;
		else{
			syN[1] = dsy*Config.wP[1];
			if(syN[1]>Config.WP[1]-1) syN[1] = 0; 
		}

		dsz = Config.part[i].pos[2]+distN*Config.part[i].ori[2];
		if(dsz < 0) szN[1]= Config.WP[2]-1;
		else{
			szN[1] = dsz*Config.wP[2];
			if(szN[1]>Config.WP[2]-1) szN[1] = 0; 
		}

		k = sxN[1] + Config.WP[0]*(syN[1]+Config.WP[1]*szN[1]); 

		Config.part[i].cellN = k;
		i += Config.part.size();
		if(Config.headP[k]==-1) Config.headP[k] = i;
		else{
		   Config.linkP[i] = Config.headP[k];
		   Config.headP[k] = i; 
		}
		i -= Config.part.size();
	}

	if(Config.Ns!=0){
		Config.WS[0] = (int)(Config.l[0]/rlistS); 
		Config.WS[1] = (int)(Config.l[1]/rlistS); 
		Config.WS[2] = (int)(Config.l[2]/rlistS); 
		
		Config.wS[0] = Config.WS[0]/Config.l[0];
		Config.wS[1] = Config.WS[1]/Config.l[1];
		Config.wS[2] = Config.WS[2]/Config.l[2];
		
		Config.headS.resize(Config.WS[0]*Config.WS[1]*Config.WS[2]); 
		for( int i=0; i<Config.headS.size(); i++) Config.headS[i] = -1;

		for( int i=Config.Nc; i<Config.part.size(); i++){ 
		    
			if( Config.part[i].pos[0] < 0 ) sx[1] = 0;
			else{
				sx[1] = Config.part[i].pos[0]*Config.wS[0];
				if(sx[1]>Config.WS[0]-1) sx[1] = Config.WS[0]-1;
			}
			if( Config.part[i].pos[1] < 0 ) sy[1] = 0;
			else{
				sy[1] = Config.part[i].pos[1]*Config.wS[1]; 
				if(sy[1]>Config.WS[1]-1) sy[1] = Config.WS[1]-1;
			}
			if( Config.part[i].pos[2] < 0 ) sz[1] = 0;
			else{
				sz[1] = Config.part[i].pos[2]*Config.wS[2];
				if(sz[1]>Config.WS[2]-1) sz[1] = Config.WS[2]-1;
			}
		   
			k = sx[1] + Config.WS[0]*(sy[1]+Config.WS[1]*sz[1]); 

			Config.part[i].cell = k;
			Config.headS[k] = i;

		}

		Config.headPS.resize(Config.WP[0]*Config.WP[1]*Config.WP[2]); 
		for( int i=0; i<Config.headPS.size(); i++) Config.headPS[i] = -1;

		for( int i=Config.Nc; i<Config.part.size(); i++){ 
		    
			if( Config.part[i].pos[0] < 0 ) sx[1] = 0;
			else{
				sx[1] = Config.part[i].pos[0]*Config.wP[0];
				if(sx[1]>Config.WP[0]-1) sx[1] = Config.WP[0]-1; 
			}
			if( Config.part[i].pos[1] < 0 ) sy[1] = 0;
			else{
				sy[1] = Config.part[i].pos[1]*Config.wP[1]; 
				if(sy[1]>Config.WP[1]-1) sy[1] = Config.WP[1]-1;
			}
			if( Config.part[i].pos[2] < 0 ) sz[1] = 0;
			else{
				sz[1] = Config.part[i].pos[2]*Config.wP[2];
				if(sz[1]>Config.WP[2]-1) sz[1] = Config.WP[2]-1; 
			}
		   
			k = sx[1] + Config.WP[0]*(sy[1]+Config.WP[1]*sz[1]); 

			Config.part[i].cellN = k;
			if(Config.headPS[k]==-1) Config.headPS[k] = i;
			else{
			   Config.linkP[i] = Config.headPS[k];
			   Config.headPS[k] = i; 
			}
		}
	} 
};
