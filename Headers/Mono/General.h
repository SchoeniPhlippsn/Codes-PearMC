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
long savestep, sstep;
long logstep, lstep;
long backupstep, bstep;
long acc;
double rho0, rhoInit;
double norm;
double Vratio;
double Vn;
double powVn;
double pos_lambda;
double ori_lambda;
double rand_lambda;
double vproc;
double r;
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

double lxn,lyn,lzn;

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
VCReport report;
std::vector<int> NCell;
std::vector<int> NPart;

class teilchen{
    public:
        
    double trans[4][4];
        private:
};

double newposx;       	// position of the particle
double newposy;       	// position of the particle
double newposz;       	// position of the particle
double newdist_orix;      	// orientation of the particle
double newdist_oriy;      	// orientation of the particle
double newdist_oriz;      	// orientation of the particle
double newpos_msdx;       	// position of the particle
double newpos_msdy;       	// position of the particle
double newpos_msdz;       	// position of the particle
int newcell;
int newsx;
int newsy;
int newsz;
int newcellN;
int newsNx;
int newsNy;
int newsNz;

double *posx;       	// position of the particle
double *posy;       	// position of the particle
double *posz;       	// position of the particle
double *dist_orix;      	// orientation of the particle
double *dist_oriy;      	// orientation of the particle
double *dist_oriz;      	// orientation of the particle
double *pos_msdx;       	// position of the particle
double *pos_msdy;       	// position of the particle
double *pos_msdz;       	// position of the particle
int *cell;
int *sx;
int *sy;
int *sz;
int *cellN;
int *sNx;
int *sNy;
int *sNz;
bool *already;


std::ifstream iFile;
std::ofstream oFile;
std::fstream fFile;


std::vector<teilchen> part;       
double lx,ly,lz;
double lx_2,ly_2,lz_2;
double rhoN;
double rhoV;
double Vbox;
double Vsys;
int Nc;
int Ns;
double msd;
std::vector<int> headP;
std::vector<int> linkP;
int WPx;
int WPy;
int WPz;
double wPx;
double wPy;
double wPz;
std::vector<bool> usedCell;
int step;
	                                                     
class teilchen MovedParticle;

int num_tri;
double max_x, max_z;

double Rx,Ry,Rz;
double Rsq, rw, rwT, rw1, ww, www, wwww;
double xlambda, xmu; 


template<class T>
inline std::string toString(const T& t){
	std::ostringstream os;
	os << t;
  	return os.str();
}

template<class T>
inline T fromString(const std::string& s){
	T t;
	std::istringstream is(s);
	is >> t;
	return t;
}
#include "PearPotential.h"
#include "System.h"
//hardWalls main
