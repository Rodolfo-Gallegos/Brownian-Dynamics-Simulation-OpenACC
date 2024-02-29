//////////////////////////////////////////
//  Funciones del Problema de Thomson   //
//////////////////////////////////////////
//  Eduardo Rodolfo Rodríguez Gallegos  //
//     Facultad de Ciencias - UASLP     //
//////////////////////////////////////////

#include "funciones.h"


/* TODO: 
    
*/


void Programa_Thomson(int N, double l_B, double D, double dt, double A, double *__restrict xm, double *__restrict ym,  double *__restrict zm){

    time_t start, finish;
    double tiempo = 0;   start = clock();                                             // Empieza a contar el tiempo

    cout<<setprecision(10);

    double U = 0.0, R = 1.0;                                                        // Potencial y radio de la esfera

    int np = 1;                                                                     // Guarda el número de procesadores utilizados


    /*     Declaración y alojamiento de arreglos      */

    double *restrict x = new double[N];                                             // Reservación de memoria dinámica para los arreglos
    double *restrict y = new double[N];    
    double *restrict z = new double[N];

    double *restrict Fx = new double[N];                                            // Reservación de memoria dinámica para los arreglos
    double *restrict Fy = new double[N];   
    double *restrict Fz = new double[N];

    for(int i = 0; i < N; i++){                                                     // Inicializa los arreglos de las Fuerzas
        Fx[i] = 0.0;    Fy[i] = 0.0;    Fz[i] = 0.0;                                
    }

    /*     Funciones iniciales     */

    Posiciones_iniciales(N, x, y, z, np);                                           // Calcula las posiciones iniciales de los iones
    
    Proyecta_esfera(N, x, y, z, R);                                                 // Proyecta los iones a la superficie de la esfera

    Separa_puntos(N, x, y, z, R);                                                   // Separa los iones iniciales


    cout<<"\n\tProblema de Thomson \t\t U = 1/r "
    <<"\t\t Radio: "<<R<<" A"<<"\t\t dt = "<<dt<<endl;                              // Cabecera del programa


    /*     Declaración de variables y objetos     */

    int it = 1e4, itDG = 1e5, con = 0;                                              // Número de iteraciones coulombianas, iteraciones de 
    const double difussion = sqrt(2 * D*dt / (A*A));                                // descenso de gradiente y contador


    /*     Variables para OpenACC    */

    int tam = N * 3;                                                                // Tamaño del arreglo que guarda los números aleatorios

    double *restrict nrand = new double[tam];                                       // Arreglo de números aleatorios, de una distribución normal 
                                                                                    // con media 0 y desviación estandar 1
    
    #pragma acc enter data copyin(x[:N], y[:N], z[:N], Fx[:N], Fy[:N], Fz[:N], \
                                    N, R, dt, l_B, D, A, difussion, nrand[:tam])    // Copia de arreglos y variables a la GPU
    
    cudaStream_t stream; 

    curandGenerator_t gen;                                                          // Generador de números aleatorios en CUDA

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);                         // Crea el generador de CUDA

    #pragma acc host_data use_device(nrand)                             
    {
        stream = (cudaStream_t) acc_get_cuda_stream(acc_async_sync);                // Consigue el stream

        curandSetPseudoRandomGeneratorSeed(gen, time(0));                           // Asigna la semilla al generador
    }

    curandSetStream(gen, stream); 

    
    /*     Main loop     */
    
    for(int n = 1; n <= it; n++){

        #pragma acc parallel loop gang present(x[:N], y[:N], z[:N], Fx[:N], Fy[:N], Fz[:N]) independent 
        for(int i = 0; i < N; i++){

            double Fxx = 0.0, Fyy = 0.0, Fzz = 0.0;

            #pragma acc loop vector reduction(+:Fxx, Fyy, Fzz) 
            for(int j = 0; j < N; j++){
                
                if(i != j){           
                    double dx = x[i] - x[j];
                    double dy = y[i] - y[j];
                    double dz = z[i] - z[j];

                    double disC = dx*dx + dy*dy + dz*dz;                            // Distancia entre iones al cuadrado

                    double dis = sqrt(disC);                                        // Distancia entre iones

                    double Fi = l_B / (A * dis*dis*dis);                            // Fuerza electrostática ion-ion

                    Fxx += Fi * dx; 
                    Fyy += Fi * dy;
                    Fzz += Fi * dz;
                }
            }

            Fx[i] = Fxx; Fy[i] = Fyy; Fz[i] = Fzz;
        }

        
        #pragma acc host_data use_device(nrand) 
        {        
            curandGenerateNormalDouble(gen, nrand, tam, 0.0, 1.0);                  // Distribución normal con media 0 y desviación estandar 1.
        }

        #pragma acc wait


        #pragma acc parallel loop present(nrand[:tam], x[:N], y[:N], z[:N], Fx[:N], Fy[:N], Fz[:N], R) independent
        for(int i = 0; i < N; i++){     

            x[i] += difussion * nrand[3*i]   + D*dt * Fx[i]/A;
            y[i] += difussion * nrand[3*i+1] + D*dt * Fy[i]/A;
            z[i] += difussion * nrand[3*i+2] + D*dt * Fz[i]/A;

            Fx[i] = 0.0;    Fy[i] = 0.0;    Fz[i] = 0.0;                            // Reinicia los arreglos de las fuerzas

            double xi = x[i], yi = y[i], zi = z[i];

            double n = sqrt(xi*xi + yi*yi + zi*zi);                                 // Magnitud del vector i
            
            x[i] = R*xi/n;  y[i] = R*yi/n;  z[i] = R*zi/n;                          // Proyecta las posiciones nuevas a la superficie de la esfera  
        }
        
        if(n % 1000 == 0 && n <= it){

            //#pragma acc update host(x[:N], y[:N], z[:N])

            //U = Potencial(N, x, y, z);                                              // Calcula el nuevo potencial

            double UT = 0.0;

            #pragma acc enter data copyin(UT)

            #pragma acc parallel loop collapse(2) reduction(+:UT) present(x[:N], y[:N], z[:N]) independent 
            for(int i = 0; i < N; i++){                                                     
                
                for(int j = 0; j < N; j++){
                    if(i != j){
                        double dx = x[i] - x[j];
                        double dy = y[i] - y[j];
                        double dz = z[i] - z[j];

                        double disC = dx*dx + dy*dy + dz*dz;                        // Distancia entre iones al cuadradp

                        double d = sqrt(disC);                                      // Distancia entre iones

                        UT += 1 / d;                                                // U = (z1*z2)/r
                    }                                                                
                }
            }

            #pragma acc exit data copyout(UT)

            UT *= 0.5;

            cout<<"\n\tU("<<N<<") = "<<UT<<"\t\tIteraciones: "<<n<<endl;
        }
    }  

    #pragma acc exit data copyout(x[:N], y[:N], z[:N])

    #pragma acc exit data delete(x[:N], y[:N], z[:N], Fx[:N], Fy[:N], Fz[:N], N, R, dt, l_B, D, A, nrand[:tam])

    U = Potencial(N, x, y, z);                                                      // Calcula el potencial de forma serial

    cout<<"\n\n\tU_Serial("<<N<<") = "<<U<<"\t\tIteraciones: "<<it<<"\n"<<endl;     // Imprime el potencial calculado con fuerza electrostática

    
    Descenso_gradiente(N, x, y, z, U, R, itDG, con);                                // Aplica el método de descenso de gradiente


    U = Potencial(N, x, y, z);                                                      // Calcula el nuevo potencial de forma serial
                        
    cout<<"\n\n\tU_DG_Serial("<<N<<") = "<<U<<"\t\tIteraciones: "<<con<<"\n"<<endl; // Imprime el nuevo potencial


    Guarda_posiciones(N, x, y, z, con, U);                                          // Crea el archivo para visualizar con Jmol


    /*     Regresa las posiciones     */

    for(int i = 0; i < N; i++){                                                     // Regresa las posiciones finales al main
        xm[i] = x[i];
        ym[i] = y[i];
        zm[i] = z[i];
    }


    /*     Liberación de memoria reservada     */
    
    delete x;   delete y;   delete z;
    delete Fx;  delete Fy;  delete Fz; 


    finish = clock();                                                               // Termina de contar el tiempo
    tiempo = (double(finish-start)/CLOCKS_PER_SEC);                                   
    cout<<"\n\tTiempo Thomson: "<<tiempo/np<<" s\n"<<endl;
}


void Posiciones_iniciales(int N, double *__restrict x, double *__restrict y,  double *__restrict z, int &np){
    
    random_device rd{};    
    mt19937 gen{rd()};    
    uniform_real_distribution<>d{-1, 1};                                            // Genera números aleatorios entre -1 y 1

    #pragma omp parallel for
    for(int i = 0; i < N; i++){         
        x[i] = d(gen); y[i] = d(gen); z[i] = d(gen);                                // Posiciones iniciales aleatorias

        if(i == 0){
            cout<<"\n\n\tNúmero de procesadores disponibles: "
                <<thread::hardware_concurrency();
            
            np = omp_get_num_threads();                                             // Número de procesadores que se usan
            
            cout<<"\t\tNúmero de procesadores en uso: "<<np<<endl;
        }
    }
}


void Proyecta_esfera(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double R){
    
    #pragma omp parallel for
    for(int i = 0; i < N; i++){                                                     // Coloca las posiciones en la superficie de la esfera
        
        double n = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);                         // Magnitud del vector i
        
        x[i] = R*x[i]/n;  y[i] = R*y[i]/n;  z[i] = R*z[i]/n;                        // Posiciones en la superficie de la esfera   
    }
}


void Separa_puntos(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double R){
    
    random_device rd{};    mt19937 gen{rd()};    uniform_real_distribution<>dis{-1, 1};

    bool Distancia_correcta = false;

    while(!Distancia_correcta){

        Distancia_correcta = true;

        #pragma omp parallel for
        for(int i = 0; i < N; i++){
            
            for(int j = 0; j < N; j++){

                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];

                double disC = dx*dx + dy*dy + dz*dz;                                // Distancia entre iones al cuadradp

                double d = sqrt(disC);                                              // Distancia entre iones
                
                if(i != j && d < 1e-4){                                             // Si un par de iones están a una distancia muy pequeña

                    x[j] = dis(gen); y[j] = dis(gen); z[j] = dis(gen);              // Da nuevas posiciones al ion j
                    
                    double n = sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j]);              
                    
                    x[j] = R*x[j]/n;  y[j] = R*y[j]/n;  z[j] = R*z[j]/n;            // Proyecta en la esfera el ion jj

                    Distancia_correcta = false;
                }
            }
        }
    }
}


double Potencial(int N, double *__restrict x, double *__restrict y,  double *__restrict z){

    double UT = 0.0;
    
    #pragma omp parallel for reduction(+:UT)
    for(int i = 0; i < N; i++){                                                     // Calcula el potencial
        
        for(int j = 0; j < N; j++){
            
            if(i != j){
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];

                double disC = dx*dx + dy*dy + dz*dz;                                // Distancia entre iones al cuadradp

                double d = sqrt(disC);                                              // Distancia entre iones

                UT += 1 / d;                                                        // U = 1/r
            }
        }
    }
    
    UT *= 0.5;                                                                      // U/2

    return UT;
}


void Gradiente_potencial(int N, double *__restrict x, double *__restrict y,  double *__restrict z, 
                            double *__restrict px, double *__restrict py,  double *__restrict pz){
    
    #pragma omp parallel for
    for(int i = 0; i < N; i++){                                                     // Calculo gradiente de U(x, y, z)

        for(int j = 0; j < N; j++){
            
            if(i != j){
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];

                double disC = dx*dx + dy*dy + dz*dz;                                // Distancia entre iones al cuadradp

                double dis = sqrt(disC);                                            // Distancia entre iones

                px[i] += dx / (dis*dis*dis);                                        // p = 1/r²
                py[i] += dy / (dis*dis*dis);
                pz[i] += dz / (dis*dis*dis);
            }
        }
    }
}


void Descenso_gradiente(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double &UT, double R, int itDG, int &conE){

    double *restrict xa, *restrict ya, *restrict za;                                    // Posiciones anteriores
    double *restrict px, *restrict py, *restrict pz;                                    // Gradiente del potencial

    xa = new double[N];                                                                 // Reservación de memoria dinámica para los arreglos
    ya = new double[N];    
    za = new double[N];

    px = new double[N];      
    py = new double[N];    
    pz = new double[N];

    for(int i = 0; i < N; i++){         
        px[i] = 0;    py[i] = 0;    pz[i] = 0;                                          // Inicializa los arrelgos del gradiente del potencial
    }
    
    Gradiente_potencial(N, x, y, z, px, py, pz);                                        // Calcula los gradientes de acuerdo a las posiciones iniciales

    double alfa = 0.5, UA = 0.0;                                                        // Tamaño de paso, potencial anterior y auxiliar

    long int aux = 0;

    double U = Potencial(N, x, y, z);                                                   // Calcula el potencial inicial


    #pragma acc enter data copyin(x[:N], y[:N], z[:N], px[:N], py[:N], pz[:N], xa[:N], ya[:N], za[:N], \
                                    N, R, alfa, UA, U)                                  // Copia de arreglos y variables a la GPU
    
    for(int con = 1; con <= itDG; con++){                                                
      
        #pragma acc parallel loop present(x[:N], y[:N], z[:N], px[:N], py[:N], pz[:N], xa[:N], ya[:N], za[:N], R, alfa, U, UA) independent
        for(int i = 0; i < N; i++){         
            
            xa[i] = x[i];                                                               // Guarda las posiciones anteriores
            ya[i] = y[i];
            za[i] = z[i];

            px[i] *= 0.5;   py[i] *= 0.5;   pz[i] *= 0.5;                               // Parte en dos la suma de los gradientes

            x[i] = xa[i] + alfa * px[i];                                                // Asigna las nuevas posiciones
            y[i] = ya[i] + alfa * py[i];
            z[i] = za[i] + alfa * pz[i];  

            px[i] = 0;    py[i] = 0;    pz[i] = 0;                                      // Reinicia los arreglos de los gradientes

            double xi = x[i], yi = y[i], zi = z[i];

            double n = sqrt(xi*xi + yi*yi + zi*zi);                                     // Magnitud del vector i
            
            x[i] = R*xi/n;  y[i] = R*yi/n;  z[i] = R*zi/n;                              // Proyecta las posiciones nuevas a la superficie de la esfera  

            if(i == 0){
                UA = U;                                                                 // Guarda el potencial anterior
                U = 0.0;
            }
        }

        //U = Potencial(N, x, y, z);                                                    // Calcula el nuevo potencial

        #pragma acc parallel loop collapse(2) reduction(+:U) present(x[:N], y[:N], z[:N], U) independent 
        for(int i = 0; i < N; i++){                                                     
            
            for(int j = 0; j < N; j++){
                if(i != j){
                    double dx = x[i] - x[j];
                    double dy = y[i] - y[j];
                    double dz = z[i] - z[j];

                    double disC = dx*dx + dy*dy + dz*dz;                                // Distancia entre iones al cuadradp

                    double d = sqrt(disC);                                              // Distancia entre iones

                    U += 0.5 / d;                                                       // U = (z1*z2)/r
                }                                                                
            }
        }


        #pragma acc parallel present(x[:N], y[:N], z[:N], xa[:N], ya[:N], za[:N], U, UA, alfa) 
        {
            if(U > UA){
                for(int i = 0; i < N; i++){                                             // Vuelve a las posiciones anteriores
                    x[i] = xa[i];
                    y[i] = ya[i];
                    z[i] = za[i];
                }

                alfa *= 0.5;                                                            // Parte en dos el tamaño de paso
            }
        }


        //Gradiente_potencial(N, x, y, z, px, py, pz);                                  // Calcula los gradientes del potencial nuevo

        #pragma acc parallel loop gang present(x[:N], y[:N], z[:N], px[:N], py[:N], pz[:N]) independent 
        for(int i = 0; i < N; i++){

            double pxx = 0.0, pyy = 0.0, pzz = 0.0;

            #pragma acc loop vector reduction(+:pxx, pyy, pzz) 
            for(int j = 0; j < N; j++){
                
                if(i != j){           
                    double dx = x[i] - x[j];
                    double dy = y[i] - y[j];
                    double dz = z[i] - z[j];

                    double disC = dx*dx + dy*dy + dz*dz;                                // Distancia entre iones al cuadrado

                    double dis = sqrt(disC);                                            // Distancia entre iones

                    pxx += dx / (dis*dis*dis);                                          // p = 1/r²
                    pyy += dy / (dis*dis*dis);
                    pzz += dz / (dis*dis*dis);
                }
            }

            px[i] = pxx; py[i] = pyy; pz[i] = pzz;
        }


        if(con > 0 && con % 1000 == 0){

            #pragma acc update host(U, UA)

            cout<<"\n\tU_DG("<<N<<") = "<<U<<"\t\tIteraciones: "<<con<<endl;            // Imprime el nuevo potencial

            if(abs(U-UA) < 1e-10){
                
                aux = con;
                con = itDG + 1;
            }
        }
    }

    #pragma acc exit data copyout(x[:N], y[:N], z[:N]) 

    #pragma acc exit data delete(x[:N], y[:N], z[:N], px[:N], py[:N], pz[:N], xa[:N], ya[:N], za[:N], R, alfa, U, UA)


    if(aux != 0){
        conE = aux;
    }else{
        conE = itDG;
    }


    delete xa;  delete ya;  delete za;                                                  // Liberación de memoria reservada
    delete px;  delete py;  delete pz;
}


void Guarda_posiciones(int NT, double *__restrict xT, double *__restrict yT,  double *__restrict zT, int itDG, double U){

    ofstream PosiT("Salidas/PosicionesThomson_"+to_string(NT)+"NT_"+to_string(itDG)+"itDG.xyz");

    if(PosiT.is_open()){                                                            // Comprueba si el archivo Posi se abrió correctamente
                
        PosiT<<" "<<NT<<endl;                                                       // Para visualizar en Jmol

        PosiT<<"\tPosiciones Thomson \t\t U/N = "<<U<<" \t\tIteraciones de descenso de gradiente: "<<itDG<<endl;

        for(int i = 0; i < NT; i++){
            PosiT<<"H\t"<<xT[i]<<"\t"<<yT[i]<<"\t"<<zT[i]<<endl;                    // Posiciones Thomson
        }

        PosiT<<"\n\n";

    }else{
        cout<<"\n\t*****   Error al abrir el archivo de Posiciones Thomson   *****"<<endl;
    }

    PosiT.close();

}
