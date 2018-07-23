double kth;
double aspect;
const double a1=0.5*sqrt(0.5);
double a2;
double a3;
double b1; 
const double two_sixth = pow(2.0,1.0/6.0); 
const double two_pi = 2*M_PI;
double distN;

//hardWalls main
double rlistP;
double rcut_P;
double rcut_P_T;
double rcut_P_II;
double maxpos;
long N; //# of particles
long seed;
long finstep;
long compstep;
long savestep;
long logstep;
long backupstep;
long acc;
double rho0, rhoInit;
double norm;
double Vratio;
double Vn;
double powVn;
double pos_lambda;
double ori_lambda;
double vproc;
double dr, r;
std::vector<double> rij;
int randP, ky, kz;
int zahl, zahl1;
double acceptance;
double sqrtz;
double phi, cosphi, sinphi;
double theta, costheta, sintheta;
int Case;
int k, kN;

int vv, newvv;

double dsx;
double dsy;
double dsz;

double dsxN;
double dsyN;
double dszN;

std::vector<std::vector<int> > s_n;
std::vector<int> s_nn (5);

std::vector<double> ln;

bool inside, in_test;
bool Compressing;
bool WallMove;
bool firstTest;
int restart;
bool initCompress;

std::string file_name;
std::string file_pre;
std::string logfile;
std::string savefile;
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
    std::vector<double> dist_ori;      	// orientation of the particle
    std::vector<double> pos_msd;       	// position of the particle
    double trans[4][4];
    int cell;
    int s[3];
    int cellN;
    int sN[3];
    bool already;

        private:
};


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
std::vector<int> headP;
std::vector<int> linkP;
std::vector<int> WP;
std::vector<double> wP;
std::vector<bool> usedCell;
int step;
void write(std::string, bool);
void read(std::string);
void writeLog(std::string);
	                                                     
class teilchen MovedParticle;

int num_tri;
double max_x, max_z;

std::vector<double> R (3);
double Rsq, rw, rwT, rw1, ww, www, wwww;
double xlambda, xmu; 


#include "Functions.h"
#include "PearPotential.h"
#include "System.h"
//hardWalls main
