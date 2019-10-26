#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define SWAP(x, y) { double temp = *x; *x = *y; *y = temp; }

/*definition and initialization of variables such as attr used in threading*/
#define NUM_THREADS 33
int counter;
pthread_attr_t attr;
pthread_mutex_t mutexcounter;

//structure to hold tree element
typedef struct vptree vptree;

struct vptree{
    int VpId;
    vptree* left;
    vptree* right;
    double *coord;
    int dim;
    double mu;
    int index;
};

//structure to pass arguments to createVpTree
typedef struct args args;

struct args{
    int dim;
    int size;
    int* index;
    double** list;
    vptree** position;
};


//calculates the distance between elements of the given part of the array assuming the last element is the vp and using the
//index array as reference
//validated
double * calculateDist(int dim, int size,double arrayList[size][dim],int index[size]){


      double* dist=(double*) malloc((size-1)*sizeof(double));
      int i;

      //using i to dereference index table and then access the original holder table of nxd
      for( i=0;i<size-1;++i){
        double distSqrd=0.0;
        int y;

        //for loop to calculate distance for all dimensions
        for(y=0;y<dim;++y){
            distSqrd+=pow(arrayList[index[i]][y]-arrayList[index[size-1]][y],2);
        }

        dist[i]=sqrt(distSqrd);

      }

       return dist;
}

// Partition using Lomuto partition scheme and parallel update of the initial index table
//validated
double* partition(int size,int* a, double* left, double* right, double* pivotIndex)
{
	// Pick pivotIndex as pivot from the array
	double pivot = *pivotIndex;
	int* pivotA = a+(pivotIndex-left);

	// Move pivot to end
	SWAP(pivotIndex, right);
	SWAP(pivotA,(a+size-1));

	// elements less than pivot will be pushed to the left of pIndex
	// elements more than pivot will be pushed to the right of pIndex
	// equal elements can go either way
	double* pIndex = left;
	int* pIndexA=a;
	int i;

	// each time we finds an element less than or equal to pivot, pIndex
	// is incremented and that element would be placed before the pivot and index table is also updated.
	for (i = 0; i < size-1; i++)
	{
		if (left[i] <= pivot)
		{
			SWAP((left+i), pIndex);
			SWAP((a+i),pIndexA);
			pIndex++;
			pIndexA++;
		}
	}

	// Move pivot to its final place
	SWAP(pIndex, right);
	SWAP(pIndexA,(a+size-1));

	// return pIndex (index of pivot element)
	return pIndex;
}

// Returns the k-th smallest element of list within left..right
// (i.e. left <= k <= right) while updating initial index table using lomuto partition
//validated
double* quickselect(int size,int* a, double* left, double* right, int k)
{
	// If the array contains only one element, return that element
	if (left == right)
		return left;

	// select a pivotIndex between left and right
	double*  pivotIndex = left + (rand() % (right - left + 1));

	pivotIndex = partition(right-left+1,a, left, right, pivotIndex);

	// The pivot is in its final sorted position
	if ((left+k-1) == pivotIndex)
		return pivotIndex;

	// if k is less than the pivot index
	else if ( (left+k-1)< pivotIndex)
		return quickselect((pivotIndex-left), a, left, pivotIndex - 1, k);

	// if k is more than the pivot index
	else
		return quickselect(right-pivotIndex,a+(pivotIndex-left)+1, pivotIndex + 1, right, k-(pivotIndex-left)-1);
}

//finds the median value - median calculated as the first of the middle elements in case of even
//# of elements - and returns that value while having rearranged the index table properly for further use in recursion
//validated
double findMedian(int dim, int size,double arrayList[size][dim],int index[size]){

    //checks whether table has only one element
    if(size==1)
        return 0.0;

    //calculates median and updates index table
    double* dist=calculateDist(dim,size,arrayList,index);
    double mu= *quickselect(size-1,index,dist,dist+(size-2),(size-1)/2);

    return mu;
}

//creates Vptree assuming vp is the last element in index and calls recursively until
//there are no points left
//validated
vptree* createVpTree(int dim,int size,int index[size],double list[size][dim]){

    if(size==0)
       return NULL;

    vptree* node=(vptree*)malloc(sizeof(vptree));
    node->VpId=index[size-1];
    node->coord=list[index[size-1]];
    node->dim=dim;
    node-> mu= findMedian(dim,size,list,index);
    node->index=index[size-1];

    //calls recursively taking into consideration whether size is
    //odd or even number
    node->left=createVpTree(dim,(size-1)/2,index,list);

    if(size%2!=0)
         node->right=createVpTree(dim,(size-1)/2,index+(size-1)/2,list);
    else
         node->right=createVpTree(dim,(size-1)/2+1,index+(size-1)/2,list);

    return node;
}


/*function to create vpTree in threaded mode*/
//nonValidated
void *createVpTreeThreaded(void* data){
     args* argstruct=(args*) data;

     if((argstruct->size)==0){
        *(argstruct->position)=NULL;
        pthread_exit(NULL);
     }

     int s=(argstruct->size);
     int* idx=(argstruct->index);
     double (*list_ptr)[argstruct->dim] = (argstruct->list);

     //create current node info structure
     vptree* node=(vptree*)malloc(sizeof(vptree));
     node->VpId=idx[s-1];
     node->coord=list_ptr[idx[s-1]];
     node->dim=(argstruct->dim);
     node-> mu= findMedian(node->dim,s,list_ptr,idx);
     node->index=idx[s-1];
     *(argstruct->position)=node;


     //lock mutex and check for number of threads
     pthread_mutex_lock(&mutexcounter);

     //procedure when thread overhead is available
     if(counter<NUM_THREADS){
        //update counter and immediately unlock mutex
        counter+=2;
        pthread_mutex_unlock(&mutexcounter);
        pthread_t callThd[2];

        //create argument struct for left and right subtree
        args argsLeft={node->dim,(s-1)/2,idx,argstruct->list,&(node->left)};
        args argsRight={node->dim,(s%2)?(s-1)/2:((s-1)/2+1),idx+(s-1)/2,list_ptr,&(node->right)};

        //create new threads and wait for them to join
        int rc1=pthread_create(&callThd[0], &attr, createVpTreeThreaded, (void *)&argsLeft);
        int rc2 =pthread_create(&callThd[1], &attr, createVpTreeThreaded, (void *)&argsRight);

        if(rc1||rc2){
            printf("error creating thread");
            exit(-1);
        }

        void* status1;
        void* status2;
        rc1=pthread_join(callThd[0], &status1);
        rc2=pthread_join(callThd[1], &status2);

        if(rc1||rc2){
            printf("error joining thread");
            exit(-1);
        }


     }
     else{
        //unlock mutex and proceed sequentially using simple vpTreeCreate
        pthread_mutex_unlock(&mutexcounter);

        node->left=createVpTree(node->dim,(s-1)/2,idx,list_ptr);

        if(s%2!=0)
         node->right=createVpTree(node->dim,(s-1)/2,idx+(s-1)/2,list_ptr);
        else
         node->right=createVpTree(node->dim,(s-1)/2+1,idx+(s-1)/2,list_ptr);
     }

     pthread_exit(NULL);
}

//prints VpTree with reference to its nodes
//Validated
void printTree(vptree* node,int* counter){

     if(node==NULL){
        printf("\n");
        return;
     }
     (*counter)++;
     int temp=*counter;
     printf("printing left->%d \n",temp);
     printTree(node->left,counter);

     printf("%lf %lf %lf ->%d\n",(node->coord)[0],(node->coord)[1],(node->coord)[2],temp);

     printf("printing right->%d\n",temp);
     printTree(node->right,counter);
}

vptree* buildvp(double* X,int n,int d){

     int* indx=(int*)malloc(sizeof(int)*n);
     int i;
     for(i=0;i<n;++i){
        indx[i]=i;
     }

     //initializing values to start executing createVpTreeThreaded
     vptree* root=NULL;

     args* values=(args*)malloc(sizeof(args));
     values->dim=d;
     values->size=n;
     values->index=indx;
     values->list=(double**)X;
     values->position=&root;


     //threading variables
     pthread_attr_init(&attr);
     pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
     pthread_mutex_init(&mutexcounter,NULL);


     counter=0;
     pthread_t starting;
     void* status;

     if(pthread_create(&starting, &attr, createVpTreeThreaded, (void *)values))
        printf("Failed to create starting thread");

     if(pthread_join(starting, &status))
        printf("starting failed to join");

     pthread_attr_destroy(&attr);
     pthread_mutex_destroy(&mutexcounter);

     printf("threads created %d\n",counter);


     return root;
}

//main function meant for testing
int main()
{
    double holder[12][3]={{1,2,3},{4,5,6},{7,8,9},{10,11,12},{13,14,15},{16,17,18},{19,20,21},{22,23,24},{25,26,27},{28,29,30},{31,32,33},{34,35,36}};

    vptree* root=NULL;
    root=buildvp((double*)holder,12,3);
    counter=0;
    printTree(root,&counter);

    return 0;
}

