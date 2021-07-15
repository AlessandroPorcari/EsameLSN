/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

// configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblocks, seed;
double delta;
bool previous;

// block average
double blocco_pot, blocco_kin, blocco_etot, blocco_temp;
double sum_pot, sum_kin, sum_etot, sum_temp;
double sum2_pot, sum2_kin, sum2_etot, sum2_temp;
double error_pot, error_kin, error_etot, error_temp;

// function g(r)
const int nbins = 100;
double gofr[nbins];
double blocco_gofr[nbins];
double sum_gofr[nbins], sum2_gofr[nbins];
double error_gofr;

//functions
void Input();
void Move(void);
void ConfFinal(void);
void ConfPrevious(void);
void ConfXYZ(int);
void Measure(void);
void GetSums(void);
void Reset(void);
void Uncertainty(int);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
