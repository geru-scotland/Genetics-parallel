/*
    fun.h
    cabeceras de las funciones utilizadas en el modulo gengrupos
****************************************************************/

extern void inicializar_centroides(float cent[][NCAR]);
extern int nuevos_centroides(float elem[][NCAR], float cent[][NCAR], int samples[], int nelem);

extern double geneticdist (float *elem1, float *elem2);
extern void nearest_cluster (int nelem, float elem[][NCAR], float cent[][NCAR], int *samples);
extern double silhouette_simple(float samples[][NCAR], struct lista_grupos *cluster_data, float centroids[][NCAR], float a[]);
extern void analisis_enfermedades (struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf);
