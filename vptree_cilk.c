#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "vptree.h"

#define SWAP(x, y) { double temp = *x; *x = *y; *y = temp; }


void swapPoints(double* array1,double* array2,int dim){
      double tmp;
      int i;

      for(i=0; i<dim; i++){
         tmp=array1[i];
         array1[i] = array2[i];
         array2[i]=tmp;
      }
}

//calculates the distance between elements of the given part of the array assuming the last element is the vp and using the
//index array as reference
//validated
void calculateDist(int dim, int size,double arrayList[size][dim],double* dist){

     int i;
     double distSqrd;
     int y;
    //using i to dereference index table and then access the original holder table of nxd
    cilk_for( i=0;i<size-1;++i){
        distSqrd=0.0;

        //for loop to calculate distance for all dimensions
        for(y=0;y<dim;++y)
        {
           distSqrd+=pow(arrayList[i][y]-arrayList[size-1][y],2);
        }
        dist[i]=sqrt(distSqrd);
    }
}


// Partition using Lomuto partition scheme and parallel update of the initial index table
//validated
double* partition(int size,int dim,int* a,double list[size][dim],double* left, double* right, double* pivotIndex)
{
	// Pick pivotIndex as pivot from the array
	double pivot = *pivotIndex;
	int* pivotA = a+(pivotIndex-left);
	double (*pivotArray)[dim]=list+(pivotIndex-left);

	// Move pivot to end
	SWAP(pivotIndex, right);
	SWAP(pivotA,(a+size-1));
	swapPoints(*pivotArray,*(list+size-1),dim);

	// elements less than pivot will be pushed to the left of pIndex
	// elements more than pivot will be pushed to the right of pIndex
	// equal elements can go either way
	pivotIndex = left;
	pivotA=a;
	pivotArray=list;

	int i;

	// each time we finds an element less than or equal to pivot, pIndex
	// is incremented and that element would be placed before the pivot and index table is also updated.
	for (i = 0; i < size-1; i++)
	{
		if (left[i] <= pivot)
		{
			SWAP((left+i), pivotIndex);
			SWAP((a+i),pivotA);
			swapPoints(*(list+i),*pivotArray,dim);
			pivotIndex++;
			pivotA++;
			pivotArray++;
		}
	}

	// Move pivot to its final place
	SWAP(pivotIndex, right);
	SWAP(pivotA,(a+size-1));
	swapPoints(*(pivotArray),*(list+size-1),dim);

	// return pIndex (index of pivot element)
	return pivotIndex;
}

// Returns the k-th smallest element of list within left..right
// (i.e. left <= k <= right) while updating initial index table using lomuto partition
//validated
double* quickselect(int size,int dim,int* a,double list[size][dim],double* left, double* right, int k)
{
	// If the array contains only one element, return that element
	if (left == right)
		return left;

	// select a pivotIndex between left and right
	double*  pivotIndex = left + (rand() % (right - left + 1));

	pivotIndex = partition(right-left+1,dim,a,list, left, right, pivotIndex);

	// The pivot is in its final sorted position
	if ((left+k-1) == pivotIndex)
		return pivotIndex;

	// if k is less than the pivot index
	else if ( (left+k-1)< pivotIndex)
		return quickselect((pivotIndex-left),dim, a,list, left, pivotIndex - 1, k);

	// if k is more than the pivot index
	else
		return quickselect(right-pivotIndex,dim,a+(pivotIndex-left)+1,list+(pivotIndex-left)+1, pivotIndex + 1, right, k-(pivotIndex-left)-1);
}


//finds the median value - median calculated as the first of the middle elements in case of even
//# of elements - and returns that value while having rearranged the index table properly for further use in recursion
//validated
double findMedian(int dim, int size,double arrayList[size][dim],int index[size],double* dist){

    //checks whether table has only one element
    if(size==1)
        return 0.0;
    //calculates median and updates index table
    calculateDist(dim,size,arrayList,dist);
    double mu= *quickselect(size-1,dim,index,arrayList, dist,dist+(size-2),(size-1)/2);

    return mu;
}

vptree * createVpTree(int dim,int size,int index[size],double list[size][dim],double* dist){

    if(size==0)
       return NULL;

    vptree* node=(vptree*)malloc(sizeof(vptree));
    node->VpId=index[size-1];
    node->coord=list[size-1];
    node-> mu= findMedian(dim,size,list,index,dist);


    //calls recursively taking into consideration whether size is
    //odd or even number
    if(size%2!=0){
         node->left= (cilk_spawn createVpTree(dim,(size-1)/2,index,list,dist));
         node->right=createVpTree(dim,(size-1)/2,index+(size-1)/2,list+(size-1)/2,dist+(size-1)/2);
         cilk_sync;
    }
    else{
         node->left=(cilk_spawn createVpTree(dim,(size-1)/2+1,index,list,dist));
         node->right=createVpTree(dim,(size-1)/2,index+(size-1)/2+1,list+(size-1)/2+1,dist+(size-1)/2+1);
         cilk_sync;
    }
    return node;
}

//function written only to meet the header criteria of the project
//thus simply making use of createVpTree plus creating the index table
//validated
vptree* buildvp(double* X,int n,int d){

     double* list=(double*)malloc(sizeof(double)*n*d);
     memcpy(list,X,sizeof(double)*n*d);

     double* dist=(double*)malloc(n*sizeof(double));

     int* index=(int*)malloc(sizeof(int)*n);
     int i;
     for(i=0;i<n;++i){
        index[i]=i;
     }

     //creating vp
     vptree* root=createVpTree(d,n,index,(double**)list,dist);

     free(dist);
     free(list);
     free(index);

     return root;
}

vptree * getInner(vptree * T)
{
	return T->left;
}
vptree * getOuter(vptree * T)
{
	return T->right;
}
double getMD(vptree * T)
{
	return T->mu;
}
double * getVP(vptree * T)
{
	return T->coord;
}
int getIDX(vptree * T)
{
	return T->VpId;
}