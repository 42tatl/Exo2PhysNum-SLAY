#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> // Se usi fun
using namespace std;

class Exercice2
{

private:
  double t, dt, tFin;
  double m, g, L;
  double Omega, kappa;
  double theta, thetadot;
  double B0, B1;
  double mu;
  double om_0, om_1, nu, Ig;

  int N_excit, nsteps_per, Nperiod;
  int sampling;
  int last;
  ofstream *outputFile;

//Test pour section de Poincare
void printOut(bool force)
  {
      const double T = (2*M_PI)/Omega; // Period
      if (std::abs(std::fmod(t, T)) < dt)
    {
      double emec = 1.0/24.0*m*pow(L,2)*pow(thetadot,2) -mu*cos(theta)*B0; // TODO: Evaluer l'energie mecanique
      double pnc  = -mu*thetadot*B1*sin(Omega*t)*sin(theta) - kappa*pow(thetadot,2); // TODO: Evaluer la puissance des forces non conservatives

      *outputFile << t << " " << theta << " " << thetadot << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }


    // TODO define angular acceleration functions and separate contributions
    // for a[0] = function of (x,t), a[1] = function of (v)
  valarray<double> acceleration(double x, double v, double t_)
  {
    valarray<double> acc = valarray<double>(2);

    acc[0] = (- 12 * mu / (m * pow(L,2)) )* (B0 + B1 * sin(Omega*t))*sin(theta); // angular acceleration depending on x and t only
    acc[1] = -12 * kappa / (m * pow(L,2)) * thetadot; // angular acceleration depending on v only

    return acc;
  }


    void step()
  {
    // TODO: implement the extended Verlet scheme Section 2.7.4
    valarray<double> a = acceleration(theta, thetadot, t);

    theta = theta + thetadot*dt + 0.5*(a[0]+a[1])*pow(dt,2); //x_j+1 (2.128)
    double thetadot_plus0_5 = thetadot+0.5*(a[0]+a[1])*dt; //v_j+1/2 (2.129)
    double a_plus1 = (- 12 * mu / (m * pow(L,2)) )* (B0 + B1 * sin(Omega*(t+dt)))*sin(theta); //a_1(x_j+1,t_j+1) in (2.130)
    double a_plus0_5 = -12 * kappa / (m * pow(L,2)) * thetadot_plus0_5; //a_2(v_j+0.5,t_j+0.5) in (2.130)

    thetadot = thetadot + 0.5*(a[0]+a_plus1)*dt + a_plus0_5*dt; //expression 2.130 
  }


public:
  Exercice2(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("confige.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    Omega    = configFile.get<double>("Omega");     // frequency of oscillating magnetic field
    kappa    = configFile.get<double>("kappa");     // coefficient for friction
    m        = configFile.get<double>("m");         // mass
    L        = configFile.get<double>("L");         // length
    B1       = configFile.get<double>("B1");     //  B1 part of magnetic fields (oscillating part amplitude)
    B0       = configFile.get<double>("B0");     //  B0 part of magnetic fields (static part amplitude)
    mu       = configFile.get<double>("mu");     //  magnetic moment 
    theta    = configFile.get<double>("theta0");    // initial condition in theta
    thetadot = configFile.get<double>("thetadot0"); // initial condition in thetadot
    sampling = configFile.get<int>("sampling");     // number of time steps between two writings on file
    N_excit  = configFile.get<int>("N_excit");      // number of periods of excitation
    Nperiod  = configFile.get<int>("Nperiod");      // number of periods of oscillation of the eigenmode
    nsteps_per= configFile.get<int>("nsteps");      // number of time step per period

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);

    // define auxiliary variables if you need/want
    double omega0 = sqrt(12*mu*B0 / (m*pow(L,2))); //ADDED
    // TODO: implement the expression for tFin, dt
    double T = 2*pi/Omega; //Period

    if(N_excit>0){
      // simulate N_excit periods of excitation
      tFin = N_excit*T;
      dt   = T/nsteps_per;
    }
    else{
      // simulate Nperiod periods of the eigenmode
      tFin = Nperiod*T;
      dt   = T/nsteps_per;
    } 
    cout << "final time is "<<"  "<< tFin << endl; 
    cout << "omega0 is "<<"  "<< omega0 << endl; //ADDED

  }
  

  ~Exercice2()
  {
    outputFile->close();
    delete outputFile;
  };

    void run()
  {
    t = 0.;
    last = 0;
    printOut(true);

    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};

int main(int argc, char* argv[])
{
  Exercice2 engine(argc, argv);
  engine.run();

  return 0;
}
