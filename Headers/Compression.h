void Compression(){

    if(fabs(rho0-Config.rhoV) > 1e-4){
        Config.step=0;
	rhoInit = Config.rhoV;
	Compressing = true;	
        while (rho0>Config.rhoV){
            acc = 0;   
         //   if(initCompress) exit(0);
            if(initCompress) Move_step();
            acceptance = static_cast<double>(acc)/N;
            Compression_step(); 
            std::cout << acceptance << std::endl;
        }
	
        std::cout <<  Config.rhoV  << " " << rho0 << std::endl;
        
        Config.rhoV = rho0;
        Config.rhoN = Config.rhoV/Config.Vsys;
        Config.Vbox = N/Config.rhoN;
        
        ln[0] = pow(Config.Vbox, 1.0/3.0)/Config.l[0];
        ln[1] = pow(Config.Vbox, 1.0/3.0)/Config.l[1];
        ln[2] = pow(Config.Vbox, 1.0/3.0)/Config.l[2];

        Config.l[0] *= ln[0];
        Config.l[1] *= ln[1];
        Config.l[2] *= ln[2];

        for( int i=0; i < Config.part.size(); i++){
                Config.part[i].pos[0] *= ln[0];
                Config.part[i].pos[1] *= ln[1];
                Config.part[i].pos[2] *= ln[2];

                Config.part[i].pos_msd[0] = 0;
                Config.part[i].pos_msd[1] = 0;
                Config.part[i].pos_msd[2] = 0;
        }

        RenewList();
       
	Config.write(savefile,1);
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
