void RenewList(){ // Initialise Neighbour List

    Config.W[0] = (int)(Config.l[0]/rlist); 
    Config.W[1] = (int)(Config.l[1]/rlist); 
    Config.W[2] = (int)(Config.l[2]/rlist); 
    
    Config.w[0] = Config.W[0]/Config.l[0];
    Config.w[1] = Config.W[1]/Config.l[1];
    Config.w[2] = Config.W[2]/Config.l[2];
	
    Config.head.resize(Config.W[0]*Config.W[1]*Config.W[2]); 
    Config.usedCell.resize(Config.W[0]*Config.W[1]*Config.W[2],false); 

	for( int i=0; i<Config.head.size(); i++) Config.head[i] = -1;

	if( Config.Ns == 0){
		Config.link.resize(2*Config.part.size()); 
		for( int i=0; i<Config.link.size(); i++) Config.link[i] = -1;

		for( int i=0; i<Config.part.size(); i++){ 
		    
			
			dsx = Config.part[i].pos[0]-Config.part[i].dist*Config.part[i].ori[0];
			if(dsx < 0) sx[1]= Config.W[0]-1;
			else{
				sx[1] = dsx*Config.w[0];
				if(sx[1]>Config.W[0]-1) sx[1] = 0; 
			}

			dsy = Config.part[i].pos[1]-Config.part[i].dist*Config.part[i].ori[1];
			if(dsy < 0) sy[1]= Config.W[1]-1;
			else{
				sy[1] = dsy*Config.w[1];
				if(sy[1]>Config.W[1]-1) sy[1] = 0; 
			}

			dsz = Config.part[i].pos[2]-Config.part[i].dist*Config.part[i].ori[2];
			if(dsz < 0) sz[1]= Config.W[2]-1;
			else{
				sz[1] = dsz*Config.w[2];
				if(sz[1]>Config.W[2]-1) sz[1] = 0; 
			}

			int k = sx[1] + Config.W[0]*(sy[1] + Config.W[1]*sz[1]); 
			
			Config.part[i].cell = k;
			if(Config.head[k]==-1) Config.head[k] = i;
			else{
			   Config.link[i] = Config.head[k];
			   Config.head[k] = i; 
			}

			dsx = Config.part[i].pos[0]+Config.part[i].dist*Config.part[i].ori[0];
			if(dsx < 0) sxN[1]= Config.W[0]-1;
			else{
				sxN[1] = dsx*Config.w[0];
				if(sxN[1]>Config.W[0]-1) sxN[1] = 0; 
			}

			dsy = Config.part[i].pos[1]+Config.part[i].dist*Config.part[i].ori[1];
			if(dsy < 0) syN[1]= Config.W[1]-1;
			else{
				syN[1] = dsy*Config.w[1];
				if(syN[1]>Config.W[1]-1) syN[1] = 0; 
			}

			dsz = Config.part[i].pos[2]+Config.part[i].dist*Config.part[i].ori[2];
			if(dsz < 0) szN[1]= Config.W[2]-1;
			else{
				szN[1] = dsz*Config.w[2];
				if(szN[1]>Config.W[2]-1) szN[1] = 0; 
			}

			k = sxN[1] + Config.W[0]*(syN[1] + Config.W[1]*szN[1]); 
			
			Config.part[i].cellN = k;
			i += Config.Nc;
			if(Config.head[k]==-1) Config.head[k] = i;
			else{
			   Config.link[i] = Config.head[k];
			   Config.head[k] = i; 
			}
			i -= Config.Nc;
		}
	}else{
		Config.link.resize(Config.part.size()); 
		for( int i=0; i<Config.link.size(); i++) Config.link[i] = -1;

		for( int i=Config.Nc; i<Config.part.size(); i++){ 
		    
			sx[1] = Config.part[i].pos[0]*Config.w[0];
			if(sx[1]==Config.W[0]) sx[1] = Config.W[0] - 1; 
			sy[1] = Config.part[i].pos[1]*Config.w[1]; 
			if(sy[1]==Config.W[1]) sy[1] = Config.W[1] - 1; 
			sz[1] = Config.part[i].pos[2]*Config.w[2]; 
			if(sz[1]==Config.W[2]) sz[1] = Config.W[2] - 1; 
		   
			int k = sx[1] + Config.W[0]*(sy[1] + Config.W[1]*sz[1]); 

			Config.part[i].cell = k;
			if(Config.head[k]==-1) Config.head[k] = i;
			else{
			   Config.link[i] = Config.head[k];
			   Config.head[k] = i; 
			}
		}
	} 
};
