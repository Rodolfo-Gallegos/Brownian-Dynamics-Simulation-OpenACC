#pragma once        

/*      #pragma once es una directiva de preprocesador de C++ que se utiliza para indicar que un archivo de encabezado (header file) 
debe incluirse una sola vez en un programa, incluso si se incluye varias veces. Esto se utiliza para evitar problemas de dependencias cíclicas 
            entre archivos de encabezado y para asegurar que un archivo de encabezado se incluya solo una vez en un programa.               */

//////////////////////////////////////////
//       Declaración de funciones       //
//////////////////////////////////////////
//  Eduardo Rodolfo Rodríguez Gallegos  //
//     Facultad de Ciencias - UASLP     //
//////////////////////////////////////////

#include <unordered_set>    // Biblioteca para manejar conjuntos
#include <algorithm>        // Decimal a fracción.
#include <openacc.h>        // Pragmas OpenACC
#include <iostream>         // cout
#include <curand.h>         // Números aleatorios con CUDA
#include <stdlib.h>         // Standard. Para gestión de memoria dinámica (malloc y free) y números pseudo-aleatorios.
#include <iomanip>          // setprecision().
#include <stdio.h>          // Standard. 
#include <fstream>          // Creación/Lectura/Escritura de archivos de texto
#include <sstream>
#include <thread>           // Procesadores disponibles
#include <random>           // Números pseudo-aleatorios en C++
#include <string>
#include <fenv.h>           // Checa divisiones por cero, operaciones inválidas y sobrecargas.
#include <cmath>            // Standard. Operaciones matemáticas.
#include <time.h>           // Standard.
#include <omp.h>            // Pragmas de OpenMP.

//#include </usr/local/cuda-11.2/targets/x86_64-linux/include/curand.h>       // Números aleatorios con CUDA en Maxwell
//#include </usr/local/cuda/targets/x86_64-linux/include/curand.h>            // Números aleatorios con CUDA en Mangore


using namespace std;


/*  Esfera de Thomson  */

void Posiciones_iniciales(int N, double *__restrict x, double *__restrict y,  double *__restrict z, int &np); 

void Proyecta_esfera(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double R);

void Separa_puntos(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double R);

double Potencial(int N, double *__restrict x, double *__restrict y,  double *__restrict z);

void Gradiente_potencial(int N, double *__restrict x, double *__restrict y,  double *__restrict z, 
                            double *__restrict px, double *__restrict py,  double *__restrict pz);

void Descenso_gradiente(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double &U, double R, int itDG, int &con);

void Programa_Thomson(int N, double l_B, double D, double dt, double A, double *__restrict xm, double *__restrict ym,  double *__restrict zm);

void Guarda_posiciones(int NT, double *__restrict xT, double *__restrict yT,  double *__restrict zT, int itDG, double U);


/*  Macroion  */

void Iones_iniciales(int N, double *__restrict x, double *__restrict y, double *__restrict z, double R, double Rc, double r_ion, int &np);

void Separa_iones(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double R, double Rc, double r_ion);


double Potencial_M(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double *__restrict v, double l_B, 
                          double *__restrict xT, double *__restrict yT, double *__restrict zT, double *__restrict vT, int NT);

void Programa_Macroion(int N, double l_B, double D, double dt, double A, double R, double Rc, int tam_dr, double pi, long int it, double r_ion, 
                              double *__restrict xT, double *__restrict yT,  double *__restrict zT, int NT, int qM, double sigma);
    
void Valencias(double *__restrict v, double *__restrict vT, int N, int NT, int qM);


void LeeArchivoPosiciones(double *__restrict xT, double *__restrict yT,  double *__restrict zT, string filename);


void Crea_archivos(int N, double *__restrict x, double *__restrict y, double *__restrict z, double *__restrict v, int NT, double *__restrict xT, 
                double *__restrict yT, double *__restrict zT, double *__restrict rho_rp, int tam_dr,  double *__restrict rho_rn, double *__restrict vT, 
                double dr, double im, double pi, double *__restrict P_r, double R, double Rc, double dt, long int it, double *__restrict V_r, double sigma);
