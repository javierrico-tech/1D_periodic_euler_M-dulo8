#ifndef _RHS_OPERATOR
#define _RHS_OPERATOR

#include "DataStructs.h"
#include "FluxFunctions.h"

template<class T>
class RHSOperator
{
public:
  RHSOperator();
  ~RHSOperator();

  virtual void eval() = 0;
  virtual void eval(DataStruct<T> &Uin) = 0;

  virtual DataStruct<T>& ref2RHS() = 0;
};

template<class T>
class Central1D : public RHSOperator<T>
{
private:

  // structure containing the RHS values
  DataStruct<T> RHS;
  
  // reference to current solution
  DataStruct<T> &U;

  // reference to mesh 
  // TODO: change to a mesh structure
  DataStruct<T> &mesh;

  // reference to flux function
  FluxFunction<T> &F;

  void evalRHS(DataStruct<T> &Uin);

public:
  Central1D(DataStruct<T> &_U, DataStruct<T> &_mesh, FluxFunction<T> &_F);
  ~Central1D();

  virtual void eval();
  virtual void eval(DataStruct<T> &Uin);

  virtual DataStruct<T>& ref2RHS();

};

template <class T>
class EulerRHSoperator {
private:
    T dx;
    EulerFlux<T> flux_func;
    
    // Almacenamiento temporal para los flujos evaluados
    DataStruct<T> f_rho;
    DataStruct<T> f_rho_u;
    DataStruct<T> f_rho_E;

public:
    // Constructor: Inicializa la distancia dx y el tamaño de los vectores temporales
    EulerRHSoperator(T dx_in, int N) : dx(dx_in), f_rho(N), f_rho_u(N), f_rho_E(N) {}

    // Evalúa las derivadas espaciales (RHS = Right Hand Side de la ecuación)
    void eval(DataStruct<T>& rho, DataStruct<T>& rho_u, DataStruct<T>& rho_E,
              DataStruct<T>& rhs_rho, DataStruct<T>& rhs_rho_u, DataStruct<T>& rhs_rho_E) {
        
        int N = rho.getSize();
        
        // 1. Calculamos los flujos físicos en todos los puntos
        flux_func.eval(rho, rho_u, rho_E, f_rho, f_rho_u, f_rho_E);

        // 2. Calculamos las diferencias centradas para cada punto
        for (int i = 0; i < N; i++) {
            // Condiciones de contorno periódicas (anillo cerrado)
            int i_plus  = (i + 1) % N;       // Vecino de la derecha
            int i_minus = (i - 1 + N) % N;   // Vecino de la izquierda
            
            // La ecuación es dU/dt + df/dx = 0  ->  dU/dt = -df/dx
            // Por tanto, RHS = - (f_{i+1} - f_{i-1}) / (2 * dx)
            rhs_rho[i]   = - (f_rho[i_plus] - f_rho[i_minus]) / (2.0 * dx);
            rhs_rho_u[i] = - (f_rho_u[i_plus] - f_rho_u[i_minus]) / (2.0 * dx);
            rhs_rho_E[i] = - (f_rho_E[i_plus] - f_rho_E[i_minus]) / (2.0 * dx);
        }
    }
};

#endif // _RHS_OPERATOR