/*
    fun.h
    cabeceras de las funciones utilizadas en el modulo gengrupos
****************************************************************/

extern void inicializar_centroides(float cent[][NCAR]);
extern int nuevos_centroides(float elem[][NCAR], float cent[][NCAR], int popul[], int nelem);

extern double gendist (float *elem1, float *elem2);
extern void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul);
extern double silhouette_simple(float elem[][NCAR], struct lista_grupos *listag, float cent[][NCAR], float a[]);
extern void analisis_enfermedades (struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf);
