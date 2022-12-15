/*
    AC - OpenMP -- SERIE
    fun_s.c
     rutinas que se utilizan en el modulo gengrupos_s.c
****************************************************************************/
#include <math.h>
#include <float.h> // DBL_MAX
#include <stdlib.h>
#include <stdio.h>
#include "defineg.h"           // definiciones

/**************************************************************************************
   1 - Funcion para calcular la distancia genetica entre dos elementos (distancia euclidea)
       Entrada:  2 elementos con NCAR caracteristicas (por referencia)
       Salida:  distancia (double)
**************************************************************************************/
double genecticdist(float *elem1, float *elem2)
{
    double total = 0.0f;

    for(int i = 0; i < NCAR; i++)
        total += pow((double)(elem2[i] - elem1[i]), 2);

	return sqrt(total);
}

/****************************************************************************************
   2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)
   Entrada:  nelem  numero de elementos, int
             elem   elementos, una matriz de tamanno MAXE x NCAR, por referencia
             cent   centroides, una matriz de tamanno NGRUPOS x NCAR, por referencia
   Salida:   popul  grupo mas cercano a cada elemento, vector de tamanno MAXE, por referencia
*****************************************************************************************/
void nearest_cluster(int nelem, float elem[][NCAR], float cent[][NCAR], int *popul)
{
	/*
	 * Cada elemento de elem[MAXE][NCAR] se tiene que comparar con los k centroides
	 * de cent[NGRUPOS][NCAR] (vectores). Ejemplo: genecticdist(elem[0], cent[0]).
	 * popul[MAXE] -> posicion 0 de este vector, estará linkado con elem[MAXE].
	 * Ejemplo: popul[0] será el cluster más cercano al elemento elem[0]
	 * Donde popul[0] puede adoptar los valores de cent
	 */
	double dist, mindist;

	for(int i = 0; i < nelem; i++){
        mindist = __DBL_MAX__;
		for(int k = 0; k < nclusters; k++) {
			dist = genecticdist(elem[i], cent[k]);
            //printf("\n Comparando elemento: %i con centroide del cluster k=: %i \n", i, k);
			if(dist < mindist){
				mindist = dist;
                popul[i] = k;
                //printf("\n popul[%i] = %i (mindist=%f)\n", i, k, mindist);
			}
		}
        //printf("\nNearest cluster for muestra %i is %i (mindist=%f)\n", i, popul[i], mindist);
	}
}

/****************************************************************************************
   3 - Funcion para calcular la calidad de la particion de clusteres.
       Ratio entre a y b. El termino a corresponde a la distancia intra-cluster.
       El termino b corresponde a la distancia inter-cluster.
   Entrada:  samples     elementos, una matriz de tamanno MAXE x NCAR, por referencia
             cluster_data   vector de NCLUSTERS structs (informacion de grupos generados), por ref.
             centroids     centroides, una matriz de tamanno NCLUSTERS x NCAR, por referencia
   Salida:   valor del CVI (double): calidad/ bondad de la particion de clusters
*****************************************************************************************/
double silhouette_simple(float samples[][NCAR], struct lista_grupos *cluster_data, float centroids[][NCAR], float a[]){
    float b[nclusters];
    for(int k = 0; k < nclusters; k++) b[k] = 0.0f;

    for(int k = 0; k < nclusters; k++){
        //printf("\nCluster %i) nelements=%i\n", k, cluster_data[k].nelems);
        for(int i = 0; i < cluster_data[k].nelems; i++){
            for(int j = i + 1; j < cluster_data[k].nelems; j++){
                a[k] += (float) genecticdist(samples[cluster_data[k].elem_index[i]],
                                             samples[cluster_data[k].elem_index[j]]);
                //printf("\n elem_index[%i] = %f y con desplazo: %f", i , samples[cluster_data[k].elem_index[i]][0], samples[cluster_data[k].elem_index[j]][0]);
            }
        }

        // Utilizo la Suma Gaussiana para obtener el número de distancias
        // medidas para los n elementos de cada clúster. n(n-1)/2 (Teoría de Grafos, cálculo de aristas totales).
        float ngauss = ((float)(cluster_data[k].nelems * (cluster_data[k].nelems - 1)) / 2);
        a[k] = cluster_data[k].nelems <= 1 ? 0 : a[k] / ngauss;
    }

    // aproximar b[i] de cada cluster
    // Para el cluster b[0] <- tendrá la media de las distancias genéticas
    // con el resto de centroides (nclusters - 1, él mismo)
    for(int k = 0; k < nclusters; k++){
        for(int j = 0; j < nclusters; j++){
            //printf("\nDISTANCE FROM BIG B: %f\n", (float)genecticdist(centroids[k], centroids[j]));
            b[k] += (float) genecticdist(centroids[k], centroids[j]);
            //printf("\nafter: %f\n", b[k]);
        }

        b[k] /=  (float)(nclusters - 1);
    }

    float max, s = 0;
    for(int k = 0; k < nclusters; k++){
        max = MAX(a[k], b[k]);
        if(max != 0){
            //printf("\na[%i]=%f, b[%i]=%f\n", k, a[k], k, b[k]);
            s += (b[k] - a[k]) / max;
        }
    }

    printf("SILHOUETEE: s=%f, nclusters=%i Calidad de cluster: %f", s, nclusters, (double)s / nclusters);
    return (double)s / nclusters;
}

/********************************************************************************************
   4 - Funcion para relizar el analisis de enfermedades
   Entrada:  cluster_data   vector de NGRUPOS structs (informacion de grupos generados), por ref.
             enf      enfermedades, una matriz de tamaño MAXE x TENF, por referencia
   Salida:   prob_enf vector de TENF structs (informacion del análisis realizado), por ref.
*****************************************************************************************/
void analisis_enfermedades(struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf)
{
	// PARA COMPLETAR
	// Realizar el análisis de enfermedades en los grupos:
	//		mediana máxima y el grupo en el que se da este máximo (para cada enfermedad)
	//		mediana mínima y su grupo en el que se da este mínimo (para cada enfermedad)
}



/***************************************************************************************************
   OTRAS FUNCIONES DE LA APLICACION
****************************************************************************************************/

void inicializar_centroides(float cent[][NCAR]){
    //printf("\n CENTROID INITS \n");
	int i, j;
	srand (147);
	for (i=0; i < nclusters; i++)
		for (j=0; j<NCAR/2; j++){
			cent[i][j] = (rand() % 10000) / 100.0;
			cent[i][j+(NCAR/2)] = cent[i][j];
		}
}

int nuevos_centroides(float elem[][NCAR], float cent[][NCAR], int popul[], int nelem){
    //printf("\n NEW CENCTROIDSS \n");
	int i, j, fin;
	double discent;
	double additions[nclusters][NCAR + 1];
	float newcent[nclusters][NCAR];

	for (i=0; i < nclusters; i++)
		for (j=0; j<NCAR+1; j++)
			additions[i][j] = 0.0;

	// acumular los valores de cada caracteristica (100); numero de elementos al final
	for (i=0; i<nelem; i++){
		for (j=0; j<NCAR; j++) additions[popul[i]][j] += elem[i][j];
		additions[popul[i]][NCAR]++;
	}

	// calcular los nuevos centroides y decidir si el proceso ha finalizado o no (en funcion de DELTA)
	fin = 1;
	for (i=0; i < nclusters; i++){
		if (additions[i][NCAR] > 0) { // ese grupo (cluster) no esta vacio
			// media de cada caracteristica
			for (j=0; j<NCAR; j++)
				newcent[i][j] = (float)(additions[i][j] / additions[i][NCAR]);

			// decidir si el proceso ha finalizado
			discent = genecticdist(&newcent[i][0], &cent[i][0]);
			if (discent > DELTA1)
				fin = 0;  // en alguna centroide hay cambios; continuar

			// copiar los nuevos centroides
			for (j=0; j<NCAR; j++)
				cent[i][j] = newcent[i][j];
		}
	}
	return fin;
}

