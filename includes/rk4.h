#ifndef _RUNGE_KUTTA
#define _RUNGE_KUTTA

/*
RUNGE KUTTA 4 CLASS

This class allows time stepping. For flexibility each RK step is taken as follows:
  initRK();
  for(s = 0; s<rk.getNumSteps(); s++)
    1. stepUi(dt)
    2. impose BCs over the current Ui which can be obtained through currentU()
    3. compute current F (i.e., Fi = F(Ui) )
    5. setFi(Fi)
  finalizeRK(dt);

*/

#include "DataStructs.h"
#include "FluxFunctions.h"
#include "RHSoperator.h"

template<class T>
class RungeKutta4 
{
  private:
    int nSteps;
    int currentStep;

    T *coeffsA, *coeffsB;

    // reference to solution (Un)
    DataStruct<T> &Un;

    // intermediate solution
    DataStruct<T> Ui;

    // RHS 
    DataStruct<T> *fi;

  public:

    // default constructor
    RungeKutta4(DataStruct<T> &_Un);

    // default destructor
    ~RungeKutta4();

    int getNumSteps();

    // initialize the RK
    void initRK();

    // finalizes the RK (updates Un)
    void finalizeRK(const T dt);

    /*
    For the step to work properly, the user must provide the appropriate F for the current Ui.
    This is done this way becase the user might want to modify the Ui or Fi so that 
    Boundary conditions can be imposed
    */
   void stepUi(T dt);
   void setFi(DataStruct<T> &_F);

   // current Ui
   DataStruct<T>* currentU();
};

template <class T>
class EulerRungeKutta4 {
private:
    int N;
    // Arrays para guardar el estado al inicio del paso
    DataStruct<T> rho_n, rho_u_n, rho_E_n;
    // Arrays para guardar los estados temporales de los sub-pasos
    DataStruct<T> rho_tmp, rho_u_tmp, rho_E_tmp;
    // Acumuladores de la media ponderada final
    DataStruct<T> sum_rho, sum_rho_u, sum_rho_E;
    
    EulerRHSoperator<T>* rhs_op;

public:
    EulerRungeKutta4(int size, EulerRHSoperator<T>* op) : 
        N(size), 
        rho_n(size), rho_u_n(size), rho_E_n(size),
        rho_tmp(size), rho_u_tmp(size), rho_E_tmp(size),
        sum_rho(size), sum_rho_u(size), sum_rho_E(size),
        rhs_op(op) {}

    void takeStep(DataStruct<T>& rho, DataStruct<T>& rho_u, DataStruct<T>& rho_E, T dt) {
        DataStruct<T> R_rho(N), R_rhou(N), R_rhoE(N);

        // Guardamos el estado base u^n y reseteamos acumuladores
        for(int i=0; i<N; i++) {
            rho_n[i] = rho[i]; rho_u_n[i] = rho_u[i]; rho_E_n[i] = rho_E[i];
            sum_rho[i] = 0.0; sum_rho_u[i] = 0.0; sum_rho_E[i] = 0.0;
        }

        // Paso 1 (k1)
        rhs_op->eval(rho_n, rho_u_n, rho_E_n, R_rho, R_rhou, R_rhoE);
        for(int i=0; i<N; i++) {
            rho_tmp[i]   = rho_n[i]   + 0.5 * dt * R_rho[i];
            rho_u_tmp[i] = rho_u_n[i] + 0.5 * dt * R_rhou[i];
            rho_E_tmp[i] = rho_E_n[i] + 0.5 * dt * R_rhoE[i];
            sum_rho[i]   += R_rho[i]; sum_rho_u[i] += R_rhou[i]; sum_rho_E[i] += R_rhoE[i];
        }

        // Paso 2 (k2)
        rhs_op->eval(rho_tmp, rho_u_tmp, rho_E_tmp, R_rho, R_rhou, R_rhoE);
        for(int i=0; i<N; i++) {
            rho_tmp[i]   = rho_n[i]   + 0.5 * dt * R_rho[i];
            rho_u_tmp[i] = rho_u_n[i] + 0.5 * dt * R_rhou[i];
            rho_E_tmp[i] = rho_E_n[i] + 0.5 * dt * R_rhoE[i];
            sum_rho[i]   += 2.0 * R_rho[i]; sum_rho_u[i] += 2.0 * R_rhou[i]; sum_rho_E[i] += 2.0 * R_rhoE[i];
        }

        // Paso 3 (k3)
        rhs_op->eval(rho_tmp, rho_u_tmp, rho_E_tmp, R_rho, R_rhou, R_rhoE);
        for(int i=0; i<N; i++) {
            rho_tmp[i]   = rho_n[i]   + dt * R_rho[i];
            rho_u_tmp[i] = rho_u_n[i] + dt * R_rhou[i];
            rho_E_tmp[i] = rho_E_n[i] + dt * R_rhoE[i];
            sum_rho[i]   += 2.0 * R_rho[i]; sum_rho_u[i] += 2.0 * R_rhou[i]; sum_rho_E[i] += 2.0 * R_rhoE[i];
        }

        // Paso 4 (k4) y actualización final u^{n+1}
        rhs_op->eval(rho_tmp, rho_u_tmp, rho_E_tmp, R_rho, R_rhou, R_rhoE);
        for(int i=0; i<N; i++) {
            rho[i]   = rho_n[i]   + (dt / 6.0) * (sum_rho[i]   + R_rho[i]);
            rho_u[i] = rho_u_n[i] + (dt / 6.0) * (sum_rho_u[i] + R_rhou[i]);
            rho_E[i] = rho_E_n[i] + (dt / 6.0) * (sum_rho_E[i] + R_rhoE[i]);
        }
    }
};

#endif // _RUNGE_KUTTA 