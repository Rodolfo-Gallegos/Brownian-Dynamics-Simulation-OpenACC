# Simulation of an electrolyte parallelized with OpenACC
Brownian dynamics simulation of an electrolyte confined around a macroion with discrete charge (Parallelized with OpenACC).
In the [thesis](Tesis_Rodolfo.pdf) you can find all the information about the simulation program, such as the theoretical framework, the tools used, the algorithms, the tests, and some runtime results as well as the behavior of the simulation. To run the simulation, you need to have [nvcc](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html) installed, which you can download [here](https://developer.nvidia.com/hpc-sdk).



Español:

Simulación de dinámica browniana de un electrolito confinado alrededor de un macroion con carga discreta (Paralelizada con OpenACC).
En la [tesis](Tesis_Rodolfo.pdf) puedes encontrar toda la información sobre el programa de simulación, como el marco teórico, las herramientas utilizadas, los algoritmos, las pruebas y algunos resultados de tiempo de ejecución así como el comportamiento de la simulación. En la [tesis](Tesis_Rodolfo.pdf) puedes encontrar toda la información sobre el programa de simulación, como el marco teórico, las herramientas utilizadas, los algoritmos, las pruebas y algunos resultados de tiempo de ejecución así como el comportamiento de la simulación.





<div style="display: flex; flex-wrap: wrap;">
    <img src="Images/15.6.NT8100.png" alt="Electrolito alrededor de un macroion que tiene su carga distribuida en 8100 sitios de la superficie de una esfera" width="500" style="margin-right: 10px;">
    <p>Electrolito alrededor de un macroion que tiene su carga distribuida en 8100 sitios de la superficie de una esfera</p>
    <img src="Images/15.0.NT1700.png" alt="Electrolito alrededor de un macroion que tiene su carga distribuida en 1700 sitios de la superficie de una esfera" width="500" style="margin-right: 10px;">
    <p>Electrolito alrededor de un macroion que tiene su carga distribuida en 1700 sitios de la superficie de una esfera</p>
    <img src="Images/14.2.NT2.png" alt="Electrolito alrededor de un macroion que tiene su carga distribuida en 2 sitios de la superficie de una esfera" width="500" style="margin-right: 10px;">
    <p>Electrolito alrededor de un macroion que tiene su carga distribuida en 2 sitios de la superficie de una esfera</p>
    <img src="Images/13.0.NT1000.png" alt="Electrolito alrededor de un macroion que tiene su carga distribuida en 1000 sitios de la superficie de una esfera" width="500" style="margin-right: 10px;">
    <p>Electrolito alrededor de un macroion que tiene su carga distribuida en 1000 sitios de la superficie de una esfera</p>
    <img src="Images/13.Tiempos_Electrolito.png" alt="Gráfico de tiempos para la simulación del electrolito alrededor de un macroion con carga discreta. OpenMP vs OpenACC vs Serial" width="500" style="margin-right: 10px;">
    <p>Gráfico de tiempos para la simulación del electrolito alrededor de un macroion con carga discreta. OpenMP vs OpenACC vs Serial</p>
    <img src="Images/11.1.Tiempos_Thomson.png" alt="Gráfico de tiempos para la simulación del problema de Thomson. OpenMP vs OpenACC vs Serial" width="500" style="margin-right: 10px;">
    <p>Gráfico de tiempos para la simulación del problema de Thomson. OpenMP vs OpenACC vs Serial</p>
    <img src="Images/10.1.NT1.png" alt="Electrolito alrededor de un macroion que tiene su carga en el centro" width="500" style="margin-right: 10px;">
    <p>Electrolito alrededor de un macroion que tiene su carga en el centro</p>
    <img src="Images/10.NT972.png" alt="Thomson problem N = 972" width="500" style="margin-right: 10px;">
    <p>Thomson problem N = 972</p>
    <img src="Images/10.NT80.png" alt="Thomson problem N = 80" width="500" style="margin-right: 10px;">
    <p>Thomson problem N = 80</p>
    <img src="Images/9.GPU.png" alt="GPU distribution" width="500" style="margin-right: 10px;">
    <p>GPU distribution</p>
    <img src="Images/8.CPU.png" alt="CPU distribution" width="500" style="margin-right: 10px;">
    <p>CPU distribution</p>
    <img src="Images/7.ElectrolitoMT.png" alt="Electrolito alrededor de un macroion con carga discreta" width="500" style="margin-right: 10px;">
    <p>Electrolito alrededor de un macroion con carga discreta</p>
    <img src="Images/6.NT1000.png" alt="Thomson problem N = 1000" width="500" style="margin-right: 10px;">
    <p>Thomson problem N = 1000</p>
    <img src="Images/6.NT120.png" alt="Thomson problem N = 120" width="500" style="margin-right: 10px;">
    <p>Thomson problem N = 120</p>
    <img src="Images/5.rho1.png" alt="spherical charge density per unit volume" width="500" style="margin-right: 10px;">
    <p>spherical charge density per unit volume</p>
    <img src="Images/3.Potencial_rc.png" alt="modified Lennard-Jones potential" width="500" style="margin-right: 10px;">
    <p>modified Lennard-Jones potential</p>
    <img src="Images/1.Mov_Browniano.png" alt="Brownian motion / Movimiento browniano" width="500" style="margin-right: 10px;">
    <p>Brownian motion / Movimiento browniano</p>
</div>

