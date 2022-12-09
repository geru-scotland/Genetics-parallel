/*
    defineg.h
    definiciones utilizadas en los modulos de la aplicacion
***********************************************************/

#define MAXE     211640   // numero de elementos (muestras)
#define MAX_GRUPOS  100      // numero de clusters
#define NCAR     40       // dimensiones de cada muestra
#define TENF     18       // tipos de enfermedad

#define DELTA1    0.01     // convergencia: cambio minimo en un centroide
#define DELTA2    0.01     // convergencia: Best First, Silhouette
#define MAXIT    10000    // convergencia: numero de iteraciones maximo

extern int ngrupos;

// estructuras de datos
//
struct lista_grupos  // informacion de los clusters
{
 int elemg[MAXE];   // indices de los elementos
 int nelemg;        // numero de elementos
};

struct analisis     // resultados del analisis de enfermedades
{
 float mmax, mmin;   // maximo y minimo de las medianas de las probabilidad de cada enfermedad
 int gmax, gmin;    // grupos con los valores maximo y minimoi de las medianas
};

