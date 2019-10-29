#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

#define numOfThreads 4
pthread_attr_t attr;

//struct to hold data for calculate dist
typedef struct distData distData;
struct distData{
   int dim;
   int size;
   double** array;
   int* idx;
   double* dist1;
   int ref;
};


//calculates the distance between elements of the given part of the array assuming the last element is the vp and using the
//index array as reference
//validated
void *calculateDist(void* data){

      distData* p=(distData*)data;
      double (*ptr)[p->dim]=p->array;

      int i;
      double distSqrd;
      //using i to dereference index table and then access the original holder table of nxd
      for( i=0;i<(p->size);++i){
        distSqrd=0.0;

        int y;

        //for loop to calculate distance for all dimensions

        for(y=0;y<(p->dim);++y){
            distSqrd+=pow((ptr)[(p->idx)[i]][y]-(ptr)[p->ref][y],2);
        }
        (p->dist1)[i]=sqrt(distSqrd);

      }

      return ((void*)0);
}

//function to calculate dist using threads with as minimal thread creation as possible
//Validated
double* calculateDistThreaded(int dim, int size,double arrayList[size][dim],int index[size]){
      int chunk=(size-1)/numOfThreads;
      double* dist=(double*)malloc((size-1)*sizeof(double));

      //variables to be used with threading
      pthread_t callThd[numOfThreads-1];
      int i;
      int rc;
      void* status;
      distData* values=(distData*)malloc(numOfThreads*sizeof(distData));

      //create threads
      for(i=0;i<(numOfThreads-1);++i){
          values[i].dim=dim;
          values[i].size=chunk;
          values[i].array=arrayList;
          values[i].idx=index+chunk*i;
          values[i].dist1=dist+chunk*i;
          values[i].ref=index[size-1];

          if(pthread_create(&callThd[i], &attr,calculateDist, (void *)&values[i])){
            printf("Thread failed to create in calculateDist");
            exit(-1);
          }

      }

      //use main thread to calculate part of dist to avoid creating non-necessary threads
      values[numOfThreads-1].dim=dim;
      values[numOfThreads-1].size=((size-1)-(numOfThreads-1)*chunk);
      values[numOfThreads-1].array=arrayList;
      values[numOfThreads-1].idx=(index+chunk*(numOfThreads-1));
      values[numOfThreads-1].dist1=(dist+chunk*(numOfThreads-1));
      values[numOfThreads-1].ref=index[size-1];


      calculateDist((void*)&values[numOfThreads-1]);

      //join threads
      for(i=0; i<(numOfThreads-1); i++) {
       rc = pthread_join(callThd[i], &status);
       if (rc) {
          printf("ERROR joining thread");
          exit(-1);
          }
      }

      //return dist to the main caller
      return dist;
}

//function for testing
int main(void){
    double holder[13][3]={{1,2,3},{4,5,6},{7,8,9},{10,11,12},{13,14,15},{16,17,18},{19,20,21},{22,23,24},{25,26,27},{28,29,30},{31,32,33},{34,35,36},{37,38,39}};
    int index[13]={0,1,2,3,4,5,6,7,8,9,10,11,12};

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    double* dist=calculateDistThreaded(3,13,holder,index);
    int i;

    for(i=0;i<12;++i)
       {printf("%lf ",dist[i]);}
    return 0;
}

