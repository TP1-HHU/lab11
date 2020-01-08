#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u1, double* const u0, const double dx,const double dt, const double xmin,
                const int N);

void analytic_solution(const double dx, const double tend, const double N,
                                const double xmin);
//---------------------------------------
int main(){

  const double tEnd = 0.15 ;


  const int N  = 64;
  const double xmin = 0;
  const double xmax = 1;
  const double dx = (xmax-xmin)/N ;
  double dt = 0.001;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* u2 = new double[N];
  double* h;
  stringstream strm;

  initialize(u1,u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N);

  analytic_solution(dx,tEnd,N,xmin);

  cout << "Nk = " << Nk << endl;

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

        // <--- some things here

        t += dt;
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
  }

  cout << "t = " << t << endl;

  delete[] u0;
  delete[] u1;
  delete[] u2;
  return 0;
}
//-----------------------------------------------
void analytic_solution(const double dx, const double tend, const double N,
                      const double xmin)
{
  ofstream out("ana");

  for(int i=0; i<N; i++)
  {
    double x = xmin + i*dx;
    double u = sin(2*M_PI*x);
    double xi = x + u*tend ;
    out << xi << "\t" << u << endl;
  }

  out.close();
}
//-----------------------------------------------

// <--- Step function

//-----------------------------------------------
void initialize(double* const u1, double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{

// <--- generate initial conditions

}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << endl;
   }
   out.close();
}
