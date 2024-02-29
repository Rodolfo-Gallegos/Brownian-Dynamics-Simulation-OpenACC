//////////////////////////////////////////
//          Programa principal          //
//////////////////////////////////////////
//  Eduardo Rodolfo Rodríguez Gallegos  //
//     Facultad de Ciencias - UASLP     //
//////////////////////////////////////////
//    Última modificación:  18/04/23    //
//////////////////////////////////////////
//               OpenACC                //
//////////////////////////////////////////

#include "funciones.h"

//#include "Esfera_Thomson.cpp"
//#include "Macroion.cpp"

/* TODO: 
    Misma configuración de Thomson
*/


int main(){
    time_t start, finish;
    start = time(nullptr);                                      // Empieza a contar el reloj

    /* Declaración de constantes */

    //const double e_0 = 1.6021766e-19;                         // Carga fundamental
    //const double k_B = 1.380649e-23;                          // Constante de Boltzman
    //const double ep_0 = 8.85418782e-12;                       // Permitividad del vacío
    //const double T = 298.0;                                   // Temperatura ambiente
    //const double ep = 78.5;                                   // Permitividad del agua (Checar en un canal nanométrico)

    const double pi = 3.14159265359;                            // PI
    const double A = 1.0e-10;                                   // Tamaño de un Ángstrom
    const double D = 1.0e-12;                                   // Coeficiente de difusión

    const double dt = 1.0e-11;                                  // Cambio en el tiempo Macroion
    const double dt_T = 1.0e-12;                              // Cambio en el tiempo Thomson

    //const double l_B = (e_0*e_0)/(4*pi*ep_0*ep*k_B*T*A);      // Longitud de Bjerrum en unidades reducidas
    const double l_B = 7.1431944;                               // Longitud de Bjerrum en unidades reducidas

    const double R = 25.0;                                      // Radio del macroion en unidades reducidas
    const double Rc = 105.9482749;                              // Radio del cascarón en unidades reducidas
    const double r_ion = 2.5;                                   // Radio de los iones en unidades reducidas

    const double sigma = 5.0;                                   // Sigma RC en unidades reducidas

    
    /* Si NT, N o qM NO son pares, se tiene que editar la función Valencias() de Macroion.cpp */
    
    const int NT = 1700;                                        // Macroión
    const int N = 1272;                                         // Número de iones
    const int qM = -100;                                        // Carga Macroion

    // it = 3e8 for good statistics 
    const long int it = 3e6;                                    // Número de iteraciones Macroion (it > im * 10)
    const int tam_dr = (Rc-R)*16;                               // Número de capas que van del macroion al cascarón
    
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);    // Busca errores en las operaciones

    
    /*       Inicio del programa        */

    #pragma acc set device_num(0)
    
    double *x, *y, *z;                                          // Posiciones Thomson

    x = new double[NT];                                         // Reservación de memoria dinámica para los arreglos, C++
    y = new double[NT];    
    z = new double[NT];


    /* Para NT = 1 */

    //x[0] = 0.0; y[0] = 0.0; z[0] = 0.0;                       

    
    /* Para NT > 1 */

    omp_set_num_threads(16);

    Programa_Thomson(NT, l_B, D, dt_T, A, x, y, z);             // Crea un archivo con NT posiciones para el Macroion y las regresa


    
    /* Para leer archivos de posiciones (creados para visualizar con Jmol) con cualquier NT 
    
    string archivo = "Salidas/PosicionesThomson_1700NT_93000itDG.xyz";    // Ubicación del archivo con las posiciones (en caso de no usar Thomson)
    
    LeeArchivoPosiciones(x, y, z, archivo);                             // Lee el archivo de posiciones de Thomson generado previamente

    */

    
    /* Interacciones entre electrolito y macroion */

    omp_set_num_threads(1);

    //Programa_Macroion(N, l_B, D, dt, A, R, Rc, tam_dr, pi, it, r_ion, x, y, z, NT, qM, sigma);  



    /*      Liberación de memoria       */

    delete x;   delete y;   delete z;


    /*        Fin del programa          */

    finish = time(nullptr);
    double tiempo = difftime(finish, start);
    cout<<"\tTiempo total: "<<tiempo<<" s\n"<<endl;

    return 0;
}
