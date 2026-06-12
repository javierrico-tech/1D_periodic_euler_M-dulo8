#ifndef _FLUX_FUNCTIONS
#define _FLUX_FUNCTIONS

#include "DataStructs.h"

// base class for flux functions
template<class T>
class FluxFunction
{
  public:
    FluxFunction();

    // in U; out F
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F) = 0;

    // same as above but at the node level
    virtual T computeFlux(const T &Ui) = 0;
};

template<class T>
class LinearFlux : public FluxFunction<T>
{
  private:
    T c;

  public:
    LinearFlux();

    // in U; out F
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F);

    // same as above but at the node level
    virtual T computeFlux(const T &Ui);
};

template <class T>
class EulerFlux {
private:
    T gamma = 1.4; // Constante de los gases
public:
    // Esta función recibe nuestras 3 variables conservadas y devuelve los 3 flujos calculados
    void eval(DataStruct<T>& rho, DataStruct<T>& rho_u, DataStruct<T>& rho_E, 
              DataStruct<T>& f_rho, DataStruct<T>& f_rho_u, DataStruct<T>& f_rho_E) {
        
        int N = rho.getSize();
        
        for (int i = 0; i < N; i++) {
            // 1. Extraemos las variables físicas intermedias (velocidad y presión)
            T u = rho_u[i] / rho[i];
            T p = (gamma - 1.0) * (rho_E[i] - 0.5 * rho[i] * u * u);

            // 2. Calculamos los flujos correspondientes a cada ecuación de Euler respetando los signos
            f_rho[i]   = rho_u[i];
            f_rho_u[i] = rho_u[i] * u + p;
            f_rho_E[i] = u * (rho_E[i] + p);
        }
    }
};
#endif // _FLUX_FUNCTIONS