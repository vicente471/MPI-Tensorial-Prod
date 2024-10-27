Calculadora de producto tensorial con algoritmo paralelo en memoria privada
-------------------------
Compilacion: mpicc -o t2.exe t2.c
-------------------------
Ejecucion: mpirun -np k ./t2.exe -P -O data.txt
donde k = Numero de procesos
      P = {V (Particion vertical), H (Particion horizontal)}
      O = {S (Silent), V (Verbose)}
Se incluyen data.txt como datos de prueba.
------------------
Maquina de prueba:

CPU: Ryzen 5 2600x
RAM: 16 gb ram
OS: Ubuntu 22.04.3

And

CPU: Ryzen 5 5600x
RAM: 16 gb ram
OS: Ubuntu 22.04.3

COMPILADORES: gcc 11.4.0, MPICH 4.0

programadores: Vicente Santos y Cristobal Gallardo

