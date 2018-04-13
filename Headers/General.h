double kth;
double aspect;
const double a1=0.5*sqrt(0.5);
double a2;
double a3;
double b1; 
const double two_sixth = pow(2.0,1.0/6.0); 
double distN;

//hardWalls main
double rlist = 2*(b1-distN)+0.02;
double rlist_2;
double rcut_SPH;
double rcut_PSPH;
double rcut_PSPH_T;
double rcut_P;
double rcut_PII;
double rcut_PT;
double maxpos;
long N; //# of particles
long seed;
long finstep;
long compstep;
long savestep;
long logstep;
long backupstep;
long acc;
double rho0, rhoInit; //density
double Vratio;
double Vn;
double powVn;
double rsphere;
double pos_lambda;
double ori_lambda;
double vproc;
double dr, r;
std::vector<double> rij;
int randP;
int zahl, zahl1;
double acceptance;
double sqrtz;
double phi, cosphi, sinphi;
int Case;

int vv, newvv;

double dsx;
double dsy;
double dsz;

double dsxN;
double dsyN;
double dszN;

std::vector<int> sx (3);
std::vector<int> sy (3);
std::vector<int> sz (3);

std::vector<int> sxN (3);
std::vector<int> syN (3);
std::vector<int> szN (3);

std::vector<double> ln;
std::vector<double> dist_ori (3) ;

bool inside;
bool Compressing;
bool WallMove;
bool firstTest;
int restart;
bool initCompress;

std::string file_name;
std::string file_pre;
std::string logfile;
std::string RHO;
std::string NC;
std::string NS;
std::string VR;
std::string KTH;
std::string ASP;

VCollide pear_mesh;
int id[2];
std::vector<int> NCell;
std::vector<int> NPart;

class teilchen{
    public:
        
    std::vector<double> pos;       	// position of the particle
    std::vector<double> ori;      	// orientation of the particle
    std::vector<double> pos_msd;       	// position of the particle
    double trans[4][4];
    bool already;
    int cell;
    int cellN;
    double dist;

        private:
};


class system{
public:
    
    std::vector<teilchen> part;       
    std::vector<double> l;
    std::vector<double> l_2;
    double rhoN;
    double rhoV;
    double Vbox;
    double Vsys;
    int Nc;
    int Ns;
    double msd;
    std::vector<int> head;
    std::vector<int> link;
    std::vector<int> W;
    std::vector<double> w;
    std::vector<bool> usedCell;
    int step;
    void write(std::string, bool);
    void read(std::string);
    void writeLog(std::string);
	                                                     
	private:
    int v,vv;
    double sqrtz,a11,a12,a13,a21,a22,a23,a31,a32,a33;
    std::vector<double> RR;
};

class system Config;
class teilchen MovedParticle;

int num_tri;
std::vector<double> bezier_x;
std::vector<double> bezier_z;
double max_x, max_z;

#include "Functions.h"
#include "PearPotential.h"
#include "System.h"
//hardWalls main
