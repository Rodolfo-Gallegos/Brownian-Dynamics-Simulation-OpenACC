#!/bin/bash

mkdir Salidas

# Para Boltzmann, Ampere y FJEstrada

banderas="nvc++ -g -O3 -acc=gpu -Mpreprocess -gpu=cc86 -lcurand -fstrict-aliasing -Minfo=accel -mp main.cpp Thomson.cpp Macroion.cpp -o $1_OpenACC"


# Compila el programa

eval $banderas


# Quita los reportes generados que no necesitamos

rm *.o


# Número de CPUs a usar

export OMP_NUM_THREADS=1


# Ejecuta el programa en segundo plano y guarda la salida

#nohup time ./$1_OpenACC > Salidas/salida_$1_OpenACC.txt &

time ./$1_OpenACC

#nvprof ./$1_OpenACC



# Para investigar y probar

#-gpu:managed -c99 


# Para mi laptop

#banderas="g++ -g -O3 -fstrict-aliasing -fopenacc -fopt-info-optimized-omp main.cpp Esfera_Thomson.cpp Macroion.cpp -o $1_$2CPUs"



# Genera archivo de texto que contiene las banderas utilizadas

# echo $banderas > "Salidas/Banderas_$1_$2CPUs.txt"



# Banderas para nvc++

# -g: Incluye información de depuración en el ejecutable para permitir la depuración del programa.
# -O3: activa la optimización de nivel 3, lo que permite al compilador realizar optimizaciones más complejas y agresivas.
# -acc=gpu: indica al compilador que use la GPU para ejecutar el código.
# -Mpreprocess: Preprocesamiento del código fuente.
# -Mcuda: Habilita la compilación de código CUDA.
# -alias=ansi: Ayuda a hacer que el código sea más portable entre diferentes sistemas y compiladores.
# -Minfo=accel: Proporcion información de depuración adicional sobre cómo se está utilizando la GPU.
# -mp: Habilita los pragmas de OpenMP.
# -lcudart indica al compilador que incluya la biblioteca de CUDA runtime en el proceso de compilación. (#include <cuda_runtime.h>)

# Las optimizaciones de nivel 3 a menudo implican la reorganización de la memoria y la reordenación de las operaciones, lo que puede afectar la coherencia de la caché y>


# Banderas para g++ 
#banderas="g++ -O3 -march=native -fopenmp -fopenmp-simd -ftree-vectorize -fopt-info-optimized-omp"

#banderas="g++ -g -O3 -fstrict-aliasing -fopenacc -fopt-info-optimized-omp"

# -O3: activa la optimización de nivel 3, lo que permite al compilador realizar optimizaciones más complejas y agresivas.
# -march=native: optimiza el código para la arquitectura de la máquina en la que se está ejecutando.
# -fopenmp: habilita el soporte para OpenMP.
# -fopenmp-simd: habilita la vectorización automática de bucles con OpenMP.
# -ftree-vectorize: habilita la vectorización automática de bucles.

