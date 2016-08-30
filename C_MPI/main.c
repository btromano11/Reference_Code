//This program creates a central difference scheme stencil and solves the resulting Ax=b problem using the Jacobi Method in parallel. Message Passing Interface and pthreads is used for parallel implementation.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MyHeader.h"

int main(int argc, char *argv[])
{
    int *local_val;
    int *local_col;
    int *local_row_ptr;
    
    int my_rank, comm_sz;
    MPI_Comm comm;
    
    MPI_Init(&argc, &argv);
    
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);
    
    int n  = atoi(argv[1]); //converts n-size input into an integer
    int local_n = n/comm_sz;
    
    int p = atoi(argv[2]);
    
    int nnz = 3*n; //calcuate number of non-zero entries for spring
    int local_nnz;
    if(my_rank == 0 || my_rank == comm_sz-1)
    {
        local_nnz = nnz/comm_sz - 1;
    }
    else
    {
        local_nnz = nnz/comm_sz;
    }
    
    //allocate memory for CRS variables
    local_val = (int*) malloc(sizeof(int)*local_nnz);
    local_col = (int*) malloc(sizeof(int)*local_nnz);
    local_row_ptr = (int*) malloc(sizeof(int)*local_n + 1);
    
    val_intialize(local_val, local_nnz, my_rank, comm_sz, comm); //fills values vector
    col_intialize(local_col, n, local_n, local_nnz, my_rank, comm_sz, comm); //fills columns vector
    row_ptr_intialize(local_row_ptr, local_n, my_rank, comm_sz, comm); //fills row_ptr vector
    calc_matrix(local_val, local_col, local_row_ptr, n, p, local_n, my_rank, comm_sz, comm); //completes Jacobi method  
    
 free(local_val);
 free(local_col);
 free(local_row_ptr);   
    
MPI_Finalize();
    
return 0;    
}



