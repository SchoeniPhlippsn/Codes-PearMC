void Compression(){

    if(fabs(rho0-rhoV) > 1e-4){

        std::cout <<  rhoV  << " " << rho0 << std::endl;
        step=0;
	rhoInit = rhoV;
	Compressing = true;	
        while (rho0>rhoV){
            acc = 0;   
         //   if(initCompress) exit(0);
            if(initCompress) Move_step();
            acceptance = static_cast<double>(acc)/N;
            Compression_step(); 
            std::cout << acceptance << std::endl;
        }
	
        std::cout <<  rhoV  << " " << rho0 << std::endl;
        
        rhoV = rho0;
        rhoN = rhoV/Vsys;
        Vbox = N/rhoN;
        
        lxn = pow(Vbox, 1.0/3.0)/lx;
        lyn = pow(Vbox, 1.0/3.0)/ly;
        lzn = pow(Vbox, 1.0/3.0)/lz;

        lx *= lxn;
        ly *= lyn;
        lz *= lzn;

        for( int i=0; i < Nc; i++){
                posx[i] *= lxn;
                posy[i] *= lyn;
                posz[i] *= lzn;

                pos_msdx[i] = 0;
                pos_msdy[i] = 0;
                pos_msdz[i] = 0;
        }

        RenewList();
       
	writeSave();
//	exit(0);
    }
    Compressing = false;
    vproc = 0.01;

    if( pos_lambda < 0){ 
	    pos_lambda= 0.02;
	    ori_lambda= 0.02;

	    zahl = 0;
	    zahl1 = 0;
	
	    while( zahl < 10 ){ 
		acc = 0;   
		for( int i = 0; i < 10; i++) Trans_step(); 
		acceptance = static_cast<double>(acc)/(10*N);
		zahl++;
		if (acceptance > 0.55){
		    pos_lambda += 0.01;
		    if( pos_lambda > maxpos) pos_lambda = maxpos;
		    zahl=0;
		}
		if (acceptance < 0.45){
		    pos_lambda -= 0.001;
		    if( pos_lambda < 0.005 ) pos_lambda = 0.005;
		    zahl=0;
		}
		
		std::cout << acceptance << "\tpos_lambda = "<< pos_lambda  << std::endl;
	    }

	    zahl = 0;
	    zahl1 = 0;
	    while( zahl < 10 ){ 
		acc = 0;   
		for( int i = 0; i < 10; i++) Rot_step(); 
		acceptance = static_cast<double>(acc)/(10*N);
		zahl++;
		if (acceptance > 0.55){
		    ori_lambda += 0.01;
		    if( ori_lambda > maxpos) ori_lambda = maxpos;
		    zahl=0;
		}
		if (acceptance < 0.45){
		    ori_lambda -= 0.001;
		    if( ori_lambda < 0.005 ) ori_lambda = 0.005;
		    zahl=0;
		}
		
		std::cout << acceptance << "\tori_lambda = " << ori_lambda << std::endl;
	    }
    }
    zahl=0;
}
