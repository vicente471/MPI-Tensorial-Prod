/***************************************************************************************
 * 
 * t2.c: producto tensorial calculado con mpi  
 *
 * Programmer: Cristobal Gallardo Cubillos & Vicente Santos Varas
 *
 * Santiago de Chile, 7/11/2023
 *
 **************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "mpi.h"

#define MASTER 0
MPI_Status status;
    
float **Matrix1, **Matrix2, **MatrixFinal, **Resultado;

void Usage(char *mess) {

    printf("\nUsage: mpirun -np k %s -P -O data.txt\n",mess);
    printf("K = numero de procesos\n");
    printf("P = {V: particion vertical, H: particion horizontal}\n");
    printf("O = {S: modo silencioso, V: modo vervoso}\n\n");
}


/*
 *
 */ //llena la matrix
float **FillMatrix(unsigned int r, unsigned int c, FILE *archivo) {
    unsigned int i, j;
    float **mat;

    mat = calloc(r, sizeof(float *));
    for (i = 0; i < r; i = i + 1) {
        mat[i] = calloc(c, sizeof(float));
    }
    for (i = 0; i < r; i = i + 1) {
        for (j = 0; j < c; j = j + 1) {
            fscanf(archivo, "%f", &mat[i][j]);
        }
    }
    return mat;
}


/*
 *
 */
void PrintMatrix(unsigned int r, unsigned int c, float **mat) {
    
   unsigned int i, j;
   
   for (i = 0; i < r; i = i + 1) {
      for (j = 0; j < c; j = j + 1)
         printf(" %.2f ",mat[i][j]);
      printf("\n");
   }
}

/*
 *
 */ //proceso de cálculo
float **Process(float **x, float **y, unsigned int r, unsigned c,unsigned rb, unsigned int cb) {

   unsigned int i, j, startRow, startCol, k, l;
   float **res;

   res = calloc(r * rb,sizeof(float* ));
   for (i = 0; i < r * rb; i = i + 1)
      res[i] = calloc(c * cb,sizeof(float*));
   for(i = 0; i < r; i = i + 1){
        for(j = 0; j < c; j = j + 1){
            startRow = i * rb;
            startCol = j * cb;
            for(k = 0; k < rb; k = k + 1){
                for(l = 0; l < cb; l = l + 1){
                    res[startRow + k][startCol + l] = x[i][j] * y[k][l];
                }
            }
        }
    }
    return res;
}
void printMatrix1D(float *matrix, int m, int k){
    printf("Matriz aplanada:\n");
    int i;
    for ( i = 0; i < m * k; i = i + 1){
        printf("%.2f ", matrix[i]);
    }
    printf("\n");
}
void FreeMatrix(unsigned int r, float **mat) {

   unsigned int i;
   
   for (i = 0; i < r; i = i + 1)
      free(mat[i]);
   free(mat);
}
float **GetSubMatrix(float **y, unsigned int start_row, unsigned int start_col, unsigned int end_row, unsigned int end_col) {
    float** x;
    unsigned int i, j;
    x = (float **)malloc((end_row - start_row ) * sizeof(float *));
    for (i = 0; i < (end_row - start_row); i = i + 1) {
        x[i] = (float *)malloc((end_col - start_col ) * sizeof(float));
    }
    for (i = start_row; i < end_row; i = i + 1) {
        for (j = start_col; j < end_col; j = j + 1) {
            x[i - start_row][j - start_col] = y[i][j];
        }
    }
    return x;
}
float* Matrix1D(float **matrix, int filas, int columnas){
    float *matriz_aplanada = (float*) malloc(filas * columnas * sizeof(float));
    int i, j, indice = 0;

    for (i = 0; i < filas; i = i + 1){    
        for ( j = 0; j < columnas; j = j + 1){
            matriz_aplanada[indice] =matrix[i][j];
            indice = indice + 1;
        }
    }
    return matriz_aplanada;
}
float** Matrix2D(float* matriz_aplanada , int filas, int columnas){
    float** matrix = (float**) malloc(filas * sizeof(float*));
    for (int i = 0; i < filas; i = i + 1) {
    matrix[i] = (float*) malloc(columnas * sizeof(float));
    }

    int indice=0, i, j;
    for (i = 0; i < filas; i++){
        for ( j =0; j < columnas; j++){
            matrix[i][j]=matriz_aplanada[indice++];
        }
    }
    return matrix;
}
int main(int argc, char **argv) {
    int m, i, j, f, k, n, mkkn, start_row, end_row, indice,start_col,end_col,p = 0,rows_per_process,colum_per_process, extra_rows, extra_colum, rank, size, me, datos[4];
    float* matrixA_1D, *matrixB_1D, *matrixC = NULL , **arreglo, **MatrixA_2D;
    char  processor_name[MPI_MAX_PROCESSOR_NAME], nombre[20]; 
    float E_cpu;
    long E_wall;  
    FILE *archivo;
    
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name,&me);

    if(rank == MASTER){                //Proceso principal (MASTER)
        time_t  ts, te;
        clock_t cs, ce;
        ts = time(NULL);
        cs = clock();

        if (strcmp(argv[1],"-V") == 0) {       
            strcpy(nombre, "Vertical");
        }else if (strcmp(argv[1],"-H") == 0){
            strcpy(nombre, "Horizontal");
        }
        if (argc == 4 && size > 1) {
            // Leer las dimensiones de la matriz
            archivo = fopen(argv[3],"r");
            if (archivo == NULL) {
                printf("No se pudo abrir el archivo.\n");
                return 0;
            }
            fscanf(archivo, "%d", &m); // Filas A
            fscanf(archivo, "%d", &k); // Columnas A y filas B
            fscanf(archivo, "%d", &n); // Columnas B
            datos[1] = k;
            datos[2] = n;
            MatrixFinal = (float **)calloc(m * k, sizeof(float *));
            for (i = 0; i < m * k; i = i + 1) {
                MatrixFinal[i] = (float *)calloc(k * n, sizeof(float));
            }
            printf("Numero de Procesos: %d, Modo particion: %s, Dimension Matriz A[%d,%d] y Matriz B[%d,%d]\n", size, nombre, m, k, k, n);  
            fflush(stdout);
            // Leer e inicializa las matrices
            Matrix1 = FillMatrix(m, k, archivo);
            Matrix2 = FillMatrix(k, n, archivo);
            fclose(archivo);
            
            Resultado = (float **)calloc(m * k, sizeof(float *));
            for (i = 0; i < m * k; i = i + 1) {
                Resultado[i] = (float *)calloc(k * n, sizeof(float));
            }
            // Muestra las matrices si se desean ver
            if (strcmp(argv[2], "-V") == 0) {
                printf(" Matriz A(%d,%d):\n\n", m, k);
                PrintMatrix(m,k,Matrix1);
                printf("\n");
                printf(" Matriz B(%d,%d):\n\n", k, n);
                PrintMatrix(k, n, Matrix2);
            }
            //Pasar la matrix de 2D a 1D
            matrixA_1D = Matrix1D(Matrix1, m, k);
            matrixB_1D = Matrix1D(Matrix2, k, n);

            if (strcmp(argv[1], "-H") == 0) {       //Modo horizontal
                rows_per_process = m / (size - 1);
                mkkn = rows_per_process * k * k * n;
                extra_rows = m % (size - 1);
                datos[0] = rows_per_process;
                datos[3] = k;
                
                for (i = 1; i < size; i = i + 1){ 
                    start_row = p;
                    end_row = p + rows_per_process;
                    if(extra_rows > 0 && i == size - 1){             //Caso que no quepan todas las filas en los nodos asignados las que sobren se llevaran al ultimo nodo
                        end_row = end_row +extra_rows;
                        rows_per_process = rows_per_process + extra_rows;
                        datos[0] = rows_per_process;
                        mkkn = rows_per_process * k * k * n;
                    }
                    arreglo = GetSubMatrix(Matrix1, start_row, 0, end_row, k);  
                    if(strcmp(argv[2], "-V") == 0){
                        printf("SubMatrix en el proceso %d\n", i);
                        PrintMatrix(rows_per_process, k, arreglo); //Matrix sent to proccess
                        fflush(stdout);
                    }
                    matrixA_1D = Matrix1D(arreglo, rows_per_process, k);
                    // Mandar matrices a los nodos
                    
                    MPI_Send(datos , 4, MPI_INT, i, i, MPI_COMM_WORLD);
                    MPI_Send(matrixA_1D, rows_per_process*k, MPI_FLOAT, i, i, MPI_COMM_WORLD);   //Se envian los datos y las matrices
                    MPI_Send(matrixB_1D, k * n, MPI_FLOAT, i, i, MPI_COMM_WORLD);
                    
                    matrixC = (float*) malloc(mkkn * sizeof(float));
                    MPI_Recv(matrixC , mkkn , MPI_FLOAT,i, i, MPI_COMM_WORLD, &status);  //Se recibe el resultado
                    
                    Resultado = Matrix2D(matrixC, rows_per_process * k, k * n);
                    for (f = 0 ;f < rows_per_process * k; f = f + 1){
                        for( j = 0; j < k * n; j = j + 1 ){
                            MatrixFinal[p * k + f][j]=Resultado[f][j];
                            
                        }
                    }
                    
                    p = p + rows_per_process;
                } 
            }else{              //Modo vertical
                colum_per_process = k / (size - 1);
                mkkn = m * colum_per_process * k * n;
                extra_colum = k % (size - 1);
                datos[0] = m;
                datos[3] = colum_per_process;

                for (i = 1; i < size; i = i + 1){ 
                    start_col = p;
                    end_col = p + colum_per_process;
                    if(extra_colum > 0 && i == size - 1){
                        end_col = end_col + extra_colum;
                        colum_per_process = colum_per_process + extra_colum;
                        datos[3] = colum_per_process;
                        mkkn = m * colum_per_process * k * n;
                    }
                    arreglo = GetSubMatrix(Matrix1, 0, start_col, m, end_col);  
                    if(strcmp(argv[2], "-V") == 0){
                        printf("SubMatrix en el proceso %d\n", i);
                        PrintMatrix(m, colum_per_process, arreglo); //Matrix sent to proccess
                        
                    }
                    matrixA_1D = Matrix1D(arreglo, m, colum_per_process);
                    // Mandar matrices a los nodos
                    
                    MPI_Send(datos, 4, MPI_INT, i, i, MPI_COMM_WORLD);
                    MPI_Send(matrixA_1D, m * colum_per_process, MPI_FLOAT, i, i, MPI_COMM_WORLD);
                    MPI_Send(matrixB_1D, k * n, MPI_FLOAT, i, i, MPI_COMM_WORLD);
                    
                    matrixC = (float*) malloc(mkkn * sizeof(float));
                    MPI_Recv(matrixC, mkkn, MPI_FLOAT,i, i, MPI_COMM_WORLD, &status);
                    
                    Resultado = Matrix2D(matrixC, m * k, colum_per_process * n);
                    for (f = 0 ;f < m * k; f = f + 1){
                        for( j = 0; j < colum_per_process * n; j = j + 1){
                            MatrixFinal[f][p * n + j] = Resultado[f][j];
                            
                        }
                        
                    }
                     
                    p = p + colum_per_process;
                }
                
            }
            // Muestra la matriz resultante si se desea
            if (strcmp(argv[2], "-V") == 0) {
                printf("\n Matriz Resultado:\n\n");
                PrintMatrix(m * k, k * n, MatrixFinal);
                fflush(stdout);
            }
            
            // Limpia la memoria utilizada para las matrices
            free(Matrix1);
            free(Matrix2);
            free(Resultado);
            free(MatrixFinal);

            ce = clock();
            te = time(NULL);
            E_wall = (long) (te - ts);
            E_cpu = (float)(ce - cs) / CLOCKS_PER_SEC;
            MPI_Barrier( MPI_COMM_WORLD);
            printf("From MASTER - Elapsed CPU Time %f Wall Time %ld \n\t Sent Mess: %d, Recv Mess: %d \n",E_cpu,E_wall,3,1); 
        } else {
            Usage(argv[0]);
        }
    }else{           //Procesos 
        //inicia el contador de tiempo=0;
        time_t  ts, te;
        clock_t cs, ce;
        ts = time(NULL);
        cs = clock();
        
        MPI_Recv(datos ,4 , MPI_INT, MASTER , rank , MPI_COMM_WORLD, &status); //Se reciben los datos
        rows_per_process = datos[0];
        k = datos[1];
        n = datos[2];
        colum_per_process = datos[3];
        mkkn = rows_per_process * colum_per_process * k * n;
        matrixA_1D = (float*) malloc(rows_per_process * colum_per_process * sizeof(float));
        matrixB_1D = (float*) malloc(k * n * sizeof(float));

        MPI_Recv(matrixA_1D, rows_per_process*colum_per_process, MPI_FLOAT, MASTER, rank, MPI_COMM_WORLD, &status);  //Se reciben las matrices
        MPI_Recv(matrixB_1D , k * n , MPI_FLOAT, MASTER, rank , MPI_COMM_WORLD, &status);
        Matrix1 = Matrix2D(matrixA_1D, rows_per_process, colum_per_process);
        Matrix2 = Matrix2D(matrixB_1D, k, n);
        Resultado = (float **)calloc(rows_per_process * k, sizeof(float *));
        
        for (i = 0; i < rows_per_process * k; i = i + 1) {
            Resultado[i] = (float *)calloc(colum_per_process * n, sizeof(float));
        }
        
        Resultado = Process(Matrix1, Matrix2, rows_per_process, colum_per_process, k, n);
        matrixC = Matrix1D(Resultado, rows_per_process * k, colum_per_process * n);
        MPI_Send(matrixC  , mkkn, MPI_FLOAT, MASTER,rank, MPI_COMM_WORLD);        // Se envia el resultado al MASTER
            
        ce = clock();
        te = time(NULL);
        E_wall = (long) (te - ts);
        E_cpu = (float)(ce - cs) / CLOCKS_PER_SEC;
        //Termina el contador y se muestra tiempo en pantalla
        MPI_Barrier( MPI_COMM_WORLD);                                       //Barrera para que no se impriman los tiempos en medio de los resultados
        printf("\nFrom %d - Elapsed CPU Time %f Wall Time %ld \n\t Sent Mess: %d, Recv Mess: %d \n",rank, E_cpu, E_wall, 1, 3); 
        
    }
    MPI_Finalize();
    return 0;
}