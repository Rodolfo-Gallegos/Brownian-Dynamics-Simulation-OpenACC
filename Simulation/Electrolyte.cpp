//////////////////////////////////////////
//       Funciones del  Macroion        //
//////////////////////////////////////////
//  Eduardo Rodolfo Rodríguez Gallegos  //
//     Facultad de Ciencias - UASLP     //
//////////////////////////////////////////

#include "funciones.h"

/* TODO: 
    
*/


void Programa_Macroion(int N, double l_B, double D, double dt, double A, double R, double Rc, int tam_dr, double pi, long int it, double r_i, 
                              double *__restrict xT, double *__restrict yT,  double *__restrict zT, int NT, int qM, double sigma){

    time_t start, mid, finish;
    long int tiempo = 0;   
    start = time(nullptr);                                                          // Empieza a contar el reloj

    cout<<setprecision(12);                                                         // Número de decimales a mostrar

    double U_Thomson = Potencial(NT, xT, yT, zT);                                   // Calcula la energía por partícula de las posiciones dadas por
                                                                                    // el problema de Thomson, para una esfera de radio 1 A.

    double Dijm = R + r_i - sigma;                                                  // Límite interior

    double Dijc = Rc - r_i + sigma;                                                 // Límite exterior

    
    /* Coloca las posiciones, dadas por el programa de Thomson, en la superficie de la esfera de radio Dijm */

    for(int i = 0; i < NT; i++){                                                    

        double n = sqrt(xT[i]*xT[i] + yT[i]*yT[i] + zT[i]*zT[i]);              

        xT[i] *= Dijm/n;  yT[i] *= Dijm/n;  zT[i] *= Dijm/n;                                 
    }

    
    /*     Declaración y alojamiento de arreglos      */

    const int O_Rc = Rc*16;                                                         // Número de capas del origen al cascarón

    double *restrict rho_rp = new double[O_Rc];                                     // Reservación de memoria dinámica para los arreglos
    double *restrict rho_rn = new double[O_Rc]; 

    double *restrict P_r = new double[tam_dr];
    double *restrict V_r = new double[tam_dr];

    double *restrict x = new double[N];                                             // Reservación de memoria dinámica para los arreglos
    double *restrict y = new double[N];    
    double *restrict z = new double[N];
    double *restrict d = new double[N];

    double *restrict Fx = new double[N];                                            // Reservación de memoria dinámica para los arreglos
    double *restrict Fy = new double[N];   
    double *restrict Fz = new double[N];

    double *restrict v   = new double[N];                                           // Reservación de memoria dinámica para los arreglos
    double *restrict vT  = new double[NT];

    
    for(int i = 0; i < O_Rc; i++){                                                  // Inicializa los arreglos para Rho(r)
        rho_rp[i] = 0.0;                   
        rho_rn[i] = 0.0;  
    }

    for(int i = 0; i < tam_dr; i++){                                                // Inicializa los arreglos para P(r) y V(r)                
        P_r[i] = 0.0;                 
        V_r[i] = 0.0;
    }

    for(int i = 0; i < N; i++){                                                     // Inicializa los arreglos de las Fuerzas
        Fx[i] = 0.0;    Fy[i] = 0.0;    Fz[i] = 0.0;                                
    }

    for(int i = 0; i < NT; i++){                                                     // Inicializa los arreglos de las Fuerzas
        vT[i] = 0.0;
    }

    int np = 1;                                                                     // Guarda el número de procesadores utilizados


    /*     Inicia los archivos que guardan el historial del potencial (UpppT) y cuando el potencial está en equilibrio (Uppp)     */

    ofstream U_ppp("Salidas/Uppp_"+to_string(N)+"N_"+to_string(NT)+
                    "NT_"+to_string(sigma)+"sigma_"+to_string(it)+"it.dat");        // Crea el archivo Uppp  

    ofstream U_pppT("Salidas/UpppT_"+to_string(N)+"N_"+to_string(NT)+
                    "NT_"+to_string(sigma)+"sigma_"+to_string(it)+"it.dat");        // Crea el archivo UpppT


    if(U_ppp.is_open()){                                                            // Abre el archivo U_ppp
        U_ppp<<"\n#\tdt\tEnergía promedio por partícula || Estables"<<endl;         // Cabecera de archivo U_ppp
    }else{
        cout<<"\n\t*****   Error al abrir el archivo de U_ppp   *****"<<endl;
    }
    
    if(U_pppT.is_open()){                                                           // Abre el archivo U_pppT
        U_pppT<<"\n#\tdt\tEnergía promedio por partícula || Total"<<endl;           // Cabecera de archivo U_pppT
    }else{
        cout<<"\n\t*****   Error al abrir el archivo de U_pppT   *****"<<endl;
    }


     /*     Funciones iniciales     */

    Iones_iniciales(N, x, y, z, R, Rc, r_i, np);                                    // Calcula las posiciones iniciales de los iones

    Valencias(v, vT, N, NT, qM);                                                    // Asigna las valencias

    Separa_iones(N, x, y, z, R, Rc, r_i);                                           // Separa los iones iniciales

    cout<<"\n\n\tEnergía por partícula de las posiciones dadas por Thomson"
        <<"para una esfera de radio 1 A: "<<U_Thomson<<endl;                        // Imprime U/NT del problema de Thomson

    cout<<"\n\n\tMacroion \tRadio: "<<R<<" A\t Radio cascarón: "<<Rc
        <<" A\t dt = "<<dt<<"\tsigma = "<<sigma<<"\tIteraciones: "<<it<<endl;       // Imprime los radios, el dt, la sigma y el número de iteraciones

    
    /*     Declaración de variables y objetos     */

    double U_prom = 0.0, U = 0.0;                                                   // Promedio de potenciales por partícula
    double Ui = Potencial_M(N, x, y, z, v, l_B, xT, yT, zT, vT, NT);                // Guarda el potencial inicial por partícula

    const double dr = (Rc - R) / tam_dr;                                            // dr -> Distancia entre capas que van del macroion al cascarón
    const double difussion = sqrt(2 * D*dt / (A*A));

    double con2 = 0;                                                                // Contadores
    int con = 0;

    long int inicio_fotos = 1e6;

    long int toma_fotos = 1e1;

    long int chequeo = 1e5;

    long int guarda_archivos = 1e8;

    
    /* Variables de las Fuerzas */
    
    double rj = r_i, Dij = r_i + rj - sigma;                                        // Constantes de núcleo repulsivo

    double eps = 1e-2, ls = Dij + pow(2, 0.166667) * sigma;         

    double lm = (Dijm + pow(2, 0.166667) * sigma);
    double lc = (Dijc - pow(2, 0.166667) * sigma);


    /*     Variables para OpenACC    */

    int tam = N * 3;                                                                // Tamaño del arreglo que guarda los números aleatorios

    double *restrict nrand = new double[tam];                                       // Arreglo de números aleatorios, de una distribución normal 
                                                                                    // con media 0 y desviación estandar 1
    
    #pragma acc enter data copyin(x[:N], y[:N], z[:N], v[:N], vT[:NT], Fx[:N], Fy[:N], Fz[:N], xT[:NT], yT[:NT], zT[:NT], d[:N], \
                                    N, NT, tam_dr, dr, Dij, Dijm, Dijc, R, Rc, l_B, D, A, qM, sigma, lm, lc, eps, \
                                    rho_rp[:O_Rc], rho_rn[:O_Rc], P_r[:tam_dr], nrand[:tam])
    
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

    for(long int n = 1; n <= it; n++){

        #pragma acc parallel loop gang present(x[:N], y[:N], z[:N], v[:N], vT[:NT], Fx[:N], Fy[:N], Fz[:N], xT[:NT], yT[:NT], zT[:NT]) independent async
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

                    double Fi = v[i] * v[j] * l_B / (A * dis*dis*dis);              // Fuerza electrostática ion-ion

                    Fxx += Fi * dx; 
                    Fyy += Fi * dy;
                    Fzz += Fi * dz;

                    if(dis < ls){                                                   // Si se acercan
                        double s = 0.0;

                        if(dis <= Dij){s = eps;}
                        else{s = dis - Dij;}

                        double Fri = 24 * (2 * pow(sigma, 12) / pow(s, 13) 
                                                - pow(sigma, 6) / pow(s, 7)) / A;   // Fuerza repulsiva
                        
                        Fxx += Fri * dx / dis;    
                        Fyy += Fri * dy / dis;
                        Fzz += Fri * dz / dis;
                    }
                }
            }

            #pragma acc loop vector reduction(+:Fxx, Fyy, Fzz) 
            for(int j = 0; j < NT; j++){
        
                double dxT = x[i] - xT[j];
                double dyT = y[i] - yT[j];
                double dzT = z[i] - zT[j];

                double disC = dxT*dxT + dyT*dyT + dzT*dzT;                          // Distancia ion-macroion al cuadrado

                double dm = sqrt(disC);                                             // Distancia ion-macroion

                double Fm = v[i] * vT[j] * l_B / (A * dm*dm*dm);                    // Fuerza electrostática ion-macroion
                
                Fxx += Fm * dxT;    
                Fyy += Fm * dyT;
                Fzz += Fm * dzT;
            }

            Fx[i] = Fxx; Fy[i] = Fyy; Fz[i] = Fzz;
        }

        
        #pragma acc parallel loop present(x[:N], y[:N], z[:N], Fx[:N], Fy[:N], Fz[:N]) independent async
        for(int i = 0; i < N; i++){
            
            double diO = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);                   // Distancia ion-origen

            if(diO < lm){
                double s = 0.0;

                if(diO <= Dijm){s = eps;}
                else{s = diO - Dijm;}

                double Frm = 24 * (2 * pow(sigma, 12) / pow(s, 13) 
                                        - pow(sigma, 6) / pow(s, 7)) / A;           // Fuerza repulsiva macroion

                Fx[i] += Frm * x[i] / diO;    
                Fy[i] += Frm * y[i] / diO;
                Fz[i] += Frm * z[i] / diO;
            }

            if(diO > lc){
                double s = 0.0;

                if(diO >= Dijc){s = eps;}
                else{s = Dijc - diO;}

                double Frc = 24 * (2 * pow(sigma, 12) / pow(s, 13) 
                                        - pow(sigma, 6) / pow(s, 7)) / A;           // Fuerza repulsiva cascarón

                Fx[i] += Frc * (-1) * x[i] / diO;    
                Fy[i] += Frc * (-1) * y[i] / diO;
                Fz[i] += Frc * (-1) * z[i] / diO;
            }
        }      

        #pragma acc host_data use_device(nrand) 
        {        
            curandGenerateNormalDouble(gen, nrand, tam, 0.0, 1.0);                  // Distribución normal con media 0 y desviación estandar 1.
        }

        #pragma acc wait

        // Calcula las nuevas posiciones (x, y, z)

        #pragma acc parallel loop present(nrand[:tam], x[:N], y[:N], z[:N], Fx[:N], Fy[:N], Fz[:N]) independent
        for(int i = 0; i < N; i++){                                                 

            x[i] = x[i] + difussion * nrand[3*i]   + D*dt * Fx[i]/A;
            y[i] = y[i] + difussion * nrand[3*i+1] + D*dt * Fy[i]/A;
            z[i] = z[i] + difussion * nrand[3*i+2] + D*dt * Fz[i]/A;

            Fx[i] = 0.0;    Fy[i] = 0.0;    Fz[i] = 0.0;                            // Reinicia los arreglos de las fuerzas
        }
        
        
        if(n % chequeo == 0){                                                       // Guarda el potencial cada 10000 iteraciones
            
            double UT = 0.0;

            double dm = Rc, dM = 0.0;                                               // Variables auxiliares

            #pragma acc enter data copyin(UT, dm, dM)

            #pragma acc parallel loop collapse(2) reduction(+:UT) present(x[:N], y[:N], z[:N], v[:N], xT[:NT], yT[:NT], zT[:NT], vT[:NT]) independent async
            for(int i = 0; i < N; i++){                                                     
                
                for(int j = 0; j < N; j++){
                    if(i != j){
                        double dx = x[i] - x[j];
                        double dy = y[i] - y[j];
                        double dz = z[i] - z[j];

                        double disC = dx*dx + dy*dy + dz*dz;                        // Distancia entre iones al cuadradp

                        double d = sqrt(disC);                                      // Distancia entre iones

                        UT += v[i] * v[j] / d;                                      // U = (z1*z2)/r
                    }                                                                
                }
            }

            #pragma acc parallel loop collapse(2) reduction(+:UT) present(x[:N], y[:N], z[:N], v[:N], xT[:NT], yT[:NT], zT[:NT], vT[:NT]) independent async
            for(int i = 0; i < N; i++){ 

                for(int j = 0; j < NT; j++){   

                    double dx = x[i] - xT[j];
                    double dy = y[i] - yT[j];
                    double dz = z[i] - zT[j];

                    double disC = dx*dx + dy*dy + dz*dz;                            // Distancia ion-macroion al cuadradp

                    double dm = sqrt(disC);                                         // Distancia ion-macroion                     

                    UT += v[i] * vT[j] / dm;                                        // U = (z1*z2)/r
                }
            }

            #pragma acc parallel loop reduction(min:dm) reduction(max:dM) present(x[:N], y[:N], z[:N]) independent
            for(int i = 0; i < N; i++){                                             // Checa que ningún ión haya salido de los límites

                double dis = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);               // Calcula la distancia ion-origen

                if(dis < dm){dm = dis;}
                if(dis > dM){dM = dis;}
            }

            #pragma acc wait

            #pragma acc exit data copyout(UT, dm, dM)

            UT *= l_B * 0.5 / (N + NT);
                
            U_pppT<<"\n\t"<<n<<"\t\t"<<UT;                                          // Escribe en el archivo U_pppT                                                        

            cout<<"\n\tn = "<<n<<"\tU("<<N<<") = "<<UT<<endl;                       // Imprime el potencial

            if(n > inicio_fotos){
                U_ppp<<"\n\t"<<n<<"\t\t"<<UT;                                       // Escribe en el archivo U_ppp
                
                U_prom += UT;

                if(n % int(1e7) == 0 && con2 != 0){
                    int horas, minutos, rest;

                    mid = time(nullptr);
                    tiempo = difftime(mid, start);

                    horas = tiempo / 3600; 
                    rest = tiempo % 3600;

                    minutos = rest / 60; 
                    rest = rest % 60;

                    cout<<"\n\n\t<U("<<N<<")> = "<<(U_prom*chequeo)/(con2*toma_fotos);    // Promedio de U/N

                    cout<<"\t\tTiempo transcurrido: "<<horas<<" horas, "<<minutos<<" minutos, "<<rest<<" segundos\n"<<endl;
                }
            }

            if(dM > Rc + r_i || dm < R - r_i){                          
                cout<<"\n\t***   Los iones salieron de los límites   ***"<<endl;
                n = it + 1;
                con++;
            }
        }


        if(n > inicio_fotos && n % toma_fotos == 0){                                // Toma una "foto" cada 10 iteraciones en el equilibrio
            
            con2++;

            #pragma acc parallel loop present(x[:N], y[:N], z[:N], d[:N], rho_rp[:O_Rc], rho_rn[:O_Rc]) independent 
            for (int i = 0; i < N; i++){
                
                d[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);                     // Distancia ion-origen

                int capa = int(d[i]/dr);                                            // Capa correspondiente a la distancia

                if(v[i] > 0){
                    #pragma acc atomic update
                    rho_rp[capa] += v[i];                                           // Guarda el número de iones positivos en la capa correspondiente 
                }else{
                    #pragma acc atomic update
                    rho_rn[capa] -= v[i];                                           // Guarda el número de iones negativos en la capa correspondiente 
                }
            }

            #pragma acc wait

            #pragma acc parallel loop gang present(x[:N], y[:N], z[:N], v[:N], d[:N], P_r[:tam_dr]) independent
            for(int i = 0; i < tam_dr; i++){
                
                double Prr = 0.0;

                #pragma acc loop vector reduction(+:Prr)
                for(int j = 0; j < N; j++){                                    
                    
                    if(d[j] <= i * dr + R){                                         // Si la distancia del ion j es menor o igual 
                        Prr += v[j];                                                // a la capa i, se suma la valencia del ion
                    }

                    if(j == 0){
                        Prr += qM;                                                  // Carga integrada / e_0
                    }
                }

                P_r[i] += Prr;
            }
        }

        
        if(n > inicio_fotos && n % guarda_archivos == 0){                           // Crea archivos cada it/5 iteraciones después del equilibrio

            #pragma acc update host(x[:N], y[:N], z[:N], rho_rp[:O_Rc], rho_rn[:O_Rc], P_r[:tam_dr])
            Crea_archivos(N, x, y, z, v, NT, xT, yT, zT, rho_rp, tam_dr, rho_rn, vT, dr, con2, pi, P_r, R, Rc, dt, n, V_r, sigma);
        }             
    }

    curandDestroyGenerator(gen);                                                    // Destruye el generador

    U_ppp.close();  U_pppT.close();                                                 // Cierra los archivos de los potenciales

    if(con == 0){

        #pragma acc exit data copyout(x[:N], y[:N], z[:N])

        U = Potencial_M(N, x, y, z, v, l_B, xT, yT, zT, vT, NT);                    // Calcula el potencial final

        cout<<"\n\n\tU("<<N<<") inicial = "<<Ui;                                    // Imprime el potencial inicial (U/N);                                             
        cout<<"\t\tU("<<N<<") final = "<<U;                                         // Imprime el potencial final

        if(con2 != 0){
            cout<<"\t\t<U("<<N<<")> = "<<(U_prom*chequeo)/(con2*toma_fotos)<<"\n"<<endl;    // Promedio de U/N
        }
    }
    
    
    /*     Liberación de memoria     */ 

    #pragma acc exit data delete(d[:N], x[:N], y[:N], z[:N], v[:N], vT[:NT], Fx[:N], Fy[:N], Fz[:N], xT[:NT], yT[:NT], zT[:NT], \
                                N, NT, tam_dr, dr, Dij, Dijm, Dijc, pi, R, Rc, l_B, D, A, qM, sigma, lm, lc, eps, \
                                rho_rp[:O_Rc], rho_rn[:O_Rc], P_r[:tam_dr], nrand[:tam])
    
    delete x;   delete y;   delete z;   delete v;                                  
    delete Fx;  delete Fy;  delete Fz;  delete vT; 
    delete P_r; delete d;   delete V_r;                         
    delete rho_rp;          delete rho_rn;    
    

    finish = time(nullptr);                                                         // Termina de contar el tiempo
    tiempo = difftime(finish, start);
    cout<<"\n\tTiempo Macroion: "<<tiempo<<" s\n\n"<<endl;                          // Tiempo de CPUs/número de CPUs
}


/*    Editar si N, NT o qM son impares    */

void Valencias(double *__restrict v, double *__restrict vT, int N, int NT, int qM){

    random_device rd;                                                               // Inicializar el generador de números aleatorios
    mt19937 gen(rd());                                                              // Generador de números aleatorios con semilla aleatoria
    uniform_int_distribution<> dis(0, NT-1);                                        // Distribución uniforme entre 0 y NT-1
    unordered_set<int> posiciones_asignadas;                                        // Conjunto de posiciones ya asignadas

    int cont = 0;                                                                   // Contador de posiciones asignadas

    /* Asigna las valencias a las posiciones dadas por un archivo, cuando la valencia es +1 o -1 para cada punto en la superficie */

    while(cont < NT/2 + qM/2){                                                      // Para NT y qM par
        
        int pos = dis(gen);                                                         // Genera una posición aleatoria del arreglo vT
        
        if(!posiciones_asignadas.count(pos)){                                       // Verifica si la posición ya ha sido asignada
            posiciones_asignadas.insert(pos);                                       // Agrega la posición al conjunto de posiciones asignadas            
            vT[pos] = 1;                                                            // Asigna valor +1 al arreglo vT en la posición dada
            cont++;                                                                 // Incrementa el contador de posiciones asignadas
        }
    }

    for(int i = 0; i < NT; i++) {                                                   // Para NT par
        if(vT[i] == 0) {
            vT[i] = -1;                                                             // Asigna valor -1 a las posiciones del arreglo vT restantes
        }
    }

    /* Asigna las valencias a las posiciones dadas por el programa de Thomson, cuando la valencia es +1 o -1 para cada punto en la superficie */

    /*for(int i = 0; i < NT; i++){                                                  // Para NT y qM par
        if(i < NT/2 + qM/2){                                                        
            vT[i] = 1;                                                              
        }else{
            vT[i] = -1;
        }
    }*/

    /*for(int i = 0; i < NT; i++){                                                  // Reparte la valencia del macroion en las NT posiciones dadas por Thomson
        vT[i] = qM / NT;                                                            
    }*/

    for(int i = 0; i < N; i++){                                                     // Asigna las valencias de los iones del electrolito
        if(i < N/2 - qM/2){                                                         // cuando la valencia es +1 o -1 y N y qM son pares.
            v[i] = 1;
        }else{
            v[i] = -1;
        }
    }
}


void Iones_iniciales(int N, double *__restrict x, double *__restrict y, double *__restrict z, double R, double Rc, double r_ion, int &np){

    random_device rd{}; mt19937 gen{rd()};                                          // rd da la semilla para el mersenne twister
    uniform_real_distribution<> d(-1.0, 1.0);                                       // [-1.0, 1.0)
    uniform_real_distribution<> r(R + r_ion, Rc - r_ion);                           // [27.5, 103.5)

    #pragma omp parallel for
    for(int i = 0; i < N; i++){         
        x[i] = d(gen); y[i] = d(gen); z[i] = d(gen);                                // Posiciones aleatorias

        double Rran = r(gen);                                                       // Radio aleatorio entre macroion y cascarón

        double n = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]); 

        x[i] = Rran*x[i]/n;  y[i] = Rran*y[i]/n;  z[i] = Rran*z[i]/n; 

        if(i == 0){
            cout<<"\n\n\tNúmero de procesadores disponibles: "<<thread::hardware_concurrency();
            np = omp_get_num_threads();
            cout<<"\t\tNúmero de procesadores en uso: "<<np<<endl;
        }

        //cout<<"\tx["<<i<<"] = "<<x[i]<<"\t\ty["<<i<<"] = "<<y[i]<<"\t\tz["<<i<<"] = "<<z[i]<<endl;
        //cout<<"\n\tx²+y²+z²["<< i <<"] = "<<sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i])<<endl;
    }
}


void Separa_iones(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double R, double Rc, double r_ion){
    random_device rd{};    mt19937 gen{rd()};    
    uniform_real_distribution<>dis{-1.0, 1.0};
    uniform_real_distribution<> r(R + r_ion, Rc - r_ion);

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
                
                if(i != j && d < 3.0 * r_ion){                                      // Si un par de iones están a una distancia muy pequeña

                    x[j] = dis(gen); y[j] = dis(gen); z[j] = dis(gen);              // Da nuevas posiciones al ion j

                    double Rran = r(gen);                                           // Radio aleatorio entre macro-ion y cascarón
                    
                    double n = sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j]);              
                    
                    x[j] = Rran*x[j]/n;  y[j] = Rran*y[j]/n;  z[j] = Rran*z[j]/n;   // Posiciones nuevas

                    Distancia_correcta = false;
                }
            }
        }
    }
}


double Potencial_M(int N, double *__restrict x, double *__restrict y,  double *__restrict z, double *__restrict v, double l_B, 
                          double *__restrict xT, double *__restrict yT, double *__restrict zT, double *__restrict vT, int NT){

    double UT = 0.0;
    
    #pragma omp parallel for reduction(+:UT)                                        // reduction evita race condition
    for(int i = 0; i < N; i++){                                                     
        
        for(int j = 0; j < N; j++){
            if(i != j){
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];

                double disC = dx*dx + dy*dy + dz*dz;                                // Distancia entre iones al cuadradp

                double d = sqrt(disC);                                              // Distancia entre iones

                UT += v[i] * v[j] / d;                                              // U = (z1*z2)/r
            }                                                                
        }

        for(int j = 0; j < NT; j++){   
            double dx = x[i] - xT[j];
            double dy = y[i] - yT[j];
            double dz = z[i] - zT[j];

            double disC = dx*dx + dy*dy + dz*dz;                                    // Distancia ion-macroion al cuadradp

            double dm = sqrt(disC);                                                 // Distancia ion-macroion                     

            UT += v[i] * vT[j] / dm;                                                // U = (z1*z2)/r
        }
    }

    UT *= l_B * 0.5 / (N + NT);                                                     // UT = l_B * (z1*z2)/ 2*r

    return UT;
}


void Crea_archivos(int N, double *__restrict x, double *__restrict y, double *__restrict z, double *__restrict v, int NT, double *__restrict xT, 
                double *__restrict yT, double *__restrict zT, double *__restrict rho_rp, int tam_dr,  double *__restrict rho_rn, double *__restrict vT, 
                double dr, double im, double pi, double *__restrict P_r, double R, double Rc, double dt, long int it, double *__restrict V_r, double sigma){


    ofstream Rho_r("Salidas/Rho_r_"+to_string(N)+"N_"+to_string(NT)+"NT_"+to_string(it)+"it_"+to_string(sigma)+"sigma.dat");

    const double N_A = 6.02214076e23, e_0 = 1.602e-19, ep_0 = 8.85418782e-12, ep = 78.5;

    if(Rho_r.is_open()){                                                            // Calcula <rho_r> y lo escribe en el archivo de Rho_r

        Rho_r<<"\n#\t<rho(r)> \tPositivos \t Negativos"<<endl;        

        for(double r = 0; r < Rc; r += dr){                  

            double dV = 4 * pi * (pow(r + dr, 3) - pow(r, 3)) / 3;                  // Volumen de las capas

            double rho_positivos = rho_rp[int(r/dr)]/(N_A * 1e-27 * im)/dV;         // rho_positivos es la densidad de iones por capa
                                                                                               
            double rho_negativos = rho_rn[int(r/dr)]/(N_A * 1e-27 * im)/dV;         // Está en unidades mol/litro

            Rho_r<<"\n\t"<<r + dr * 0.5<<"\t\t"<<rho_positivos
                 <<"\t\t"<<rho_negativos;                                           // Escribe <rho(r)> | Positivos | Negativos
        }                                                               
    }else{
        cout<<"\n\t*****   Error al abrir el archivo de Rho_r   *****"<<endl;
    }

    Rho_r.close();


    ofstream Qin("Salidas/Qin_"+to_string(N)+"N_"+to_string(NT)+"NT_"+to_string(it)+"it_"+to_string(sigma)+"sigma.dat");

    if(Qin.is_open()){                                                              // Comprueba si el archivo Qin se abrió correctamente
                
        Qin<<"\n#\tDistancia\tCarga integrada <Q(r)>"<<endl;                        // Calcula <Q(r)> y lo escribe en el archivo de Qin

        for(int i = 0; i < tam_dr; i++){                   
            Qin<<"\n\t"<< i*dr + R<<"\t\t"<< P_r[i]/im;
        }

    }else{
        cout<<"\n\t*****   Error al abrir el archivo de Qin   *****"<<endl;
    }

    Qin.close();

    
    double V = 0.0;

    ofstream V_rr("Salidas/V_"+to_string(N)+"N_"+to_string(NT)+"NT_"+to_string(it)+"it_"+to_string(sigma)+"sigma.dat");

    if(V_rr.is_open()){                                                             // Comprueba si el archivo V_r se abrió correctamente
                
        V_rr<<"\n#\tr\t\tPotencial electrostático V(r)"<<endl;                      // Calcula V(r) y lo escribe en el archivo V

        for(int i = tam_dr - 1; i >= 0; i--){    
            double r = R + i*dr + dr/2;                                             // Calcula r (distancia origen-capa)

            V += (P_r[i] * dr) / (r*r*im);                                          // Calcula P(r)/r²

            V_r[i] = (V * e_0 * 1e3 * 1e10) / (4*pi*ep_0*ep);                       // Guarda el potencial electrostático en mV
        }

        for(int i = 0; i < tam_dr; i++){                   
            V_rr<<"\n\t"<< R + i*dr<<"\t\t"<< V_r[i];
        }

    }else{
        cout<<"\n\t*****   Error al abrir el archivo de V_r   *****"<<endl;
    }

    V_rr.close();


    ofstream Posi("Salidas/Posiciones_"+to_string(N)+"N_"+to_string(NT)+"NT_"+to_string(it)+"it_"+to_string(sigma)+"sigma.xyz");

    if(Posi.is_open()){                                                             // Comprueba si el archivo Posi se abrió correctamente
                
        Posi<<" "<<N+NT<<endl;                                                      // Para visualizar en Jmol

        Posi<<"\tMacroion \t Radio: "<<R<<" A\t Radio cascarón: "<<Rc<<" A\t dt = "<<dt<<"\tIteraciones: "<<it<<endl;

        for(int i = 0; i < N; i++){
            if(v[i] > 0){
                Posi<<"N\t"<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<endl;                    // iones con valencia positiva
            }else{
                Posi<<"O\t"<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<endl;                    // iones con valencia negativa
            }
        }

        for(int i = 0; i < NT; i++){
            if(vT[i] > 0){
                Posi<<"K\t"<<xT[i]<<"\t"<<yT[i]<<"\t"<<zT[i]<<endl;                 // iones con valencia positiva del Macroion
            }else{
                Posi<<"P\t"<<xT[i]<<"\t"<<yT[i]<<"\t"<<zT[i]<<endl;                 // iones con valencia negativa del Macroion
            }
        }

        Posi<<"\n\n";

    }else{
        cout<<"\n\t*****   Error al abrir el archivo de Posiciones   *****"<<endl;
    }

    Posi.close();

}


void LeeArchivoPosiciones(double *__restrict xT, double *__restrict yT,  double *__restrict zT, string filename){
    
    ifstream archivo(filename);             // Abre el archivo

    if (!archivo.is_open()) {
        cout << "No se pudo abrir el archivo " << filename << endl;
    }

    int NT;                                 
    string line;                            
    
    getline(archivo, line);                 // Lee la primera línea del archivo que contiene el número de posiciones
    
    stringstream(line) >> NT;               // Guarda el número de posiciones en NT

    getline(archivo, line);                 // Lee la segunda línea del archivo

    for(int i = 0; i < NT; i++){
        
        getline(archivo, line);             // Lee la siguiente línea del archivo
        
        stringstream ss(line);              // Guarda la línea en el objeto ss

        string color;                       
        double x, y, z;
        
        ss >> color >> x >> y >> z;         // Separa los valores de la línea guardada

        xT[i] = x;  yT[i] = y;  zT[i] = z;  // Asgina los valores al arreglo correspondiente
    }

    archivo.close();                        // Cierra el archivo

    /*for (int i = 0; i < NT; i++) {
        cout << "H\t" << xT[i] << "\t" << yT[i] << "\t" << zT[i] << endl;           // Imprimir los valores de xT, yT y zT
    }*/

}
