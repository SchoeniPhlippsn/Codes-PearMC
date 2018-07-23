void Compression(){

    if(fabs(rho0-rhoV) > 1e-4){
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
        
        ln[0] = pow(Vbox, 1.0/3.0)/l[0];
        ln[1] = pow(Vbox, 1.0/3.0)/l[1];
        ln[2] = pow(Vbox, 1.0/3.0)/l[2];

        l[0] *= ln[0];
        l[1] *= ln[1];
        l[2] *= ln[2];

        for( int i=0; i < part.size(); i++){
                part[i].pos[0] *= ln[0];
                part[i].pos[1] *= ln[1];
                part[i].pos[2] *= ln[2];

                part[i].pos_msd[0] = 0;
                part[i].pos_msd[1] = 0;
                part[i].pos_msd[2] = 0;
        }

        RenewList();
       
	write(savefile,1);
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
		for( int i = 0; i < 10; i++) Move_step(); 
		acceptance = static_cast<double>(acc)/(10*N);
		zahl++;
		if (acceptance > 0.55){
		    pos_lambda += 0.01;
		    if( pos_lambda > maxpos) pos_lambda = maxpos;
		    ori_lambda += 0.01;
		    if( ori_lambda > maxpos) ori_lambda = maxpos;
		    zahl=0;
		}
		if (acceptance < 0.45){
		    pos_lambda -= 0.001;
		    if( pos_lambda < 0.005 ) pos_lambda = 0.005;
		    ori_lambda -= 0.001;
		    if( ori_lambda < 0.005 ) ori_lambda = 0.005;
		    zahl=0;
		}
		
		std::cout << acceptance << "\tpos_lambda = "<< pos_lambda << "\tand\tori_lambda = " << ori_lambda << std::endl;
	    }
    }
    zahl=0;
}
