#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MyHeader.h"

//adds the non-zero values of the A matrix
void val_intialize(int *local_val, int local_nnz, int my_rank, int comm_sz, MPI_Comm comm)
{
    int bshift = 0;
    int eshift = 0;
    
    //fills in beginning values for process 0
    if(my_rank == 0)
    {
    (local_val)[0] = 2;
    (local_val)[1] = -1;
    bshift = 2;
    }
    
    //fills end values for last process
    if(my_rank == comm_sz-1)
    {
    (local_val)[local_nnz - 2] = -1;
    (local_val)[local_nnz - 1] = 2;
    eshift = 2;
    }

    
    //fills middle(repeated) values
    int i;
    for(i=bshift; i< local_nnz-eshift; i+=3)
    {
        if(i+2 < local_nnz)
        {
            (local_val)[i] = -1;
            (local_val)[i+1] = 2;
            (local_val)[i+2] = -1;
        }
        else
        {
            break;
        }
    }
}

//Fills the column number for each coresponding value
void col_intialize(int *local_col, int n, int local_n, int local_nnz, int my_rank, int comm_sz, MPI_Comm comm)
{
    int bshift = 0;
    int eshift = 0;
    int jshift = 1;
    
    //columns for beginning values in process 0
    if(my_rank == 0)
    {
    (local_col)[0] = 0;
    (local_col)[1] = 1;
    bshift = 2;
    jshift = 0;
    }
    
    //fills in end values for last process
    if(my_rank == comm_sz-1)
    {
    (local_col)[local_nnz - 2] = n - 2;
    (local_col)[local_nnz - 1] = n - 1;
    eshift = 2;
    }
    
    int i;
    int j = local_n*my_rank - jshift;
    
    //columns for middle values
    for(i=bshift; i< local_nnz-eshift; i+=3)
    {
        if(i+2 < local_nnz)
        {
            (local_col)[i] = j;
            (local_col)[i+1] = j + 1;
            (local_col)[i+2] = j + 2;
            
            j++;
        }
        else
        {
            break;
        }
    }
}

//fills the row pointers for each row
void row_ptr_intialize(int *local_row_ptr, int local_n, int my_rank, int comm_sz, MPI_Comm comm)
{
    int bshift = 1;
    int eshift = 0;
    
    //row pointers for first row
    (local_row_ptr)[0] = 0;
    if(my_rank == 0)
    {
    (local_row_ptr)[1] = 2;
    bshift++;
    }
    
    if(my_rank == comm_sz-1)
    { 
    eshift++;
    }
    
    int i;
    for(i=bshift; i<local_n+1-eshift; i++)
    {
        (local_row_ptr)[i] = (local_row_ptr)[i-1] + 3;
    }
    
    //row pointer for last row
    if(my_rank == comm_sz-1)
    {
    (local_row_ptr)[local_n] = (local_row_ptr)[i-1] + 2;
    }
}

//implements the Jacobi method
void calc_matrix(int *local_val, int *local_col, int *local_row_ptr, int n, int p, int local_n, int my_rank, int comm_sz, MPI_Comm comm)
{   
    double D_inv = .5; //D inverse matrix
    int b = 1; //b matrix
    double *local_x = (double*) malloc(sizeof(double)*local_n); //holds x values calculated by process
    double *calc_x = (double*) malloc(sizeof(double)*(local_n + 2)); //holds x values used for multiplication
    double *local_resids = (double*) malloc(sizeof(double)*local_n); //holds r values calculated by process
    double *local_y = (double*) malloc(sizeof(double)*local_n); //temporary holder for x matrix data
    double *local_Ax = (double*) malloc(sizeof(double)*local_n); //Ax matrix
    
    int i; //loop counter
    int k; //loop counter
    int M=8000*n; //max iterations
    double local_r = 0; //b-Ax squares counter
    double global_r;
    double residual=1; //current residual value
    double tau = .0001; //tolerance
    double residual_min = tau*sqrt(n); //convergence criterion
    
    int shift=0;
    int hold = 1;
    
    double *global_x;
   
    //timer variables
    double itimet = 0;
    double itimee;
    double itimes;
    
    double mtimet = 0;
    double mtimee;
    double mtimes;
    
    double utimet = 0;
    double utimee;
    double utimes;
    
    double itimew;
    double *itime = (double*) malloc(sizeof(double)*comm_sz);
    double mtimew;
    double *mtime = (double*) malloc(sizeof(double)*comm_sz);
    double utimew;
    double *utime = (double*) malloc(sizeof(double)*comm_sz);
    
    //initializes arrays
    for(i=0; i<local_n; i++)
    {
        local_Ax[i] = 0;
        local_x[i] = 0;
    }
    
    //calculates iterations of x
    itimes = get_time();
    int z;
    for(z=1; z<M; z++)
    {
        //matrix multiplication 
        //temp multiplication variable gets local values of x
        for(i=1; i<local_n+1; i++)
        {
            calc_x[i] = local_x[i-1];
        }
        
        //passes messages only if more than one process is running
        if(comm_sz>1)
        {
            //temp multiplication variable gets one value of x from above process
            if(my_rank==0)
            {
                MPI_Send(&local_x[local_n-1], 1, MPI_DOUBLE, my_rank+1, 0, comm);
            }
            else if(my_rank == comm_sz-1)
            {
                MPI_Recv(&calc_x[0], 1, MPI_DOUBLE, my_rank-1, 0, comm, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Recv(&calc_x[0], 1, MPI_DOUBLE, my_rank-1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Send(&local_x[local_n-1], 1, MPI_DOUBLE, my_rank+1, 0, comm);
            }

            //temp multiplication variable gets one value of x from below process
            if(my_rank==0)
            {
                MPI_Recv(&calc_x[local_n+1], 1, MPI_DOUBLE, my_rank+1, 0, comm, MPI_STATUS_IGNORE);
            }
            else if(my_rank == comm_sz-1)
            {
                MPI_Send(&local_x[0], 1, MPI_DOUBLE, my_rank-1, 0, comm);
            }
            else
            {
                MPI_Recv(&calc_x[local_n+1], 1, MPI_DOUBLE, my_rank+1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Send(&local_x[0], 1, MPI_DOUBLE, my_rank-1, 0, comm);
            }
        }
        
        mtimes = get_time();
        #pragma omp parallel for num_threads(p) \
         private(k, shift)
        for(i=0; i<local_n; i++)
        {
            int k;
            for(k=local_row_ptr[i]; k<local_row_ptr[i+1]; k++)
            {
                if(i==0 && my_rank==0)
                {
                    local_Ax[i] = local_Ax[i] + (local_val)[k]*calc_x[i+shift+hold];
                    shift++;
                }
                else
                {
                    local_Ax[i] = local_Ax[i] + (local_val)[k]*calc_x[i+shift];
                    shift++; 
                }
            }
            shift=0;   
        }

        //calculates residual
        for(i=0; i<local_n; i++)
        {
            local_resids[i] = (b - local_Ax[i]);
            local_r += pow(local_resids[i], 2);
        }
        mtimee = get_time();
        mtimet += mtimee - mtimes;
        
        MPI_Allreduce(&local_r, &global_r, 1, MPI_DOUBLE, MPI_SUM, comm);
        
        if(sqrt(global_r)<residual_min)
        {
            break;
        }
        
        local_r = 0;
        global_r = 0;
        
        //updates x
        utimes = get_time();
        for(i=0; i<local_n; i++)
        {
            local_y[i] = local_x[i] + D_inv*local_resids[i]; //calculates new x value, temporarily places it in y
        }
        
        //fills x values, resets temp placeholder and Ax
        for(i=0; i<local_n; i++)
        {
            local_x[i] = local_y[i];
            local_Ax[i] = 0; 
        }
        utimee = get_time();
        utimet += utimee - utimes;
        
        //writes files for first 2 iterations
        if(z==1 || z==2)
        {   
            if(my_rank==0)
            {
            global_x = (double*) malloc(sizeof(double)*n); 
            MPI_Gather(local_x, local_n, MPI_DOUBLE, global_x, local_n, MPI_DOUBLE, 0, comm);
            write_x_file(global_x, n, z);
            free(global_x);
            }
            else
            {
                MPI_Gather(local_x, local_n, MPI_DOUBLE, global_x, local_n, MPI_DOUBLE, 0, comm);
            }
        }
    }
    itimee = get_time();
    itimet += itimee - itimes;
    
    //writes final solution
    if(my_rank==0)
    {
    global_x = (double*) malloc(sizeof(double)*n); 
    MPI_Gather(local_x, local_n, MPI_DOUBLE, global_x, local_n, MPI_DOUBLE, 0, comm);
    write_x_file(global_x, n, z);
    free(global_x);
    }
    else
    {
        MPI_Gather(local_x, local_n, MPI_DOUBLE, global_x, local_n, MPI_DOUBLE, 0, comm);
    }

    free(local_x);
    free(local_y);
    free(local_Ax); 
    free(calc_x);
    free(local_resids);
    
   //gathers times from all processes
    MPI_Gather(&itimet, 1, MPI_DOUBLE, itime, 1, MPI_DOUBLE, 0, comm);
    MPI_Gather(&mtimet, 1, MPI_DOUBLE, mtime, 1, MPI_DOUBLE, 0, comm);
    MPI_Gather(&utimet, 1, MPI_DOUBLE, utime, 1, MPI_DOUBLE, 0, comm);
    
    //gets max time for each timer, writes to file
    if(my_rank == 0)
    {   
        itimew = itime[0];
        mtimew = mtime[0];
        utimew = utime[0];
        
        for(i=1; i<comm_sz; i++)
        {
            if(itime[i]>itimew)
            {
                itimew = itime[i];
            }
            if(mtime[i]>mtimew)
            {
                mtimew = mtime[i];
            }
            if(utime[i]>utimew)
            {
                utimew = utime[i];
            }
        }
        
        write_time_file(itimew, mtimew, utimew);
    }
    
    free(itime);
    free(mtime);
    free(utime);
}

//writes files for the x1, x2 and x* matrices
void write_x_file(double *x, int n, int z)
{
    FILE *fout;
    if(z==1)
    {
        fout = fopen("x1.txt", "w");
    }
    else if(z==2)
    {
        fout = fopen("x2.txt", "w");
    }
    else
    {
        fout = fopen("xfinal.txt", "w");
    }
    
    int i;
    for(i=0; i<n; i++)
    {
        fprintf(fout, "%le\n", (x)[i]);
    }
    fclose(fout);
}

double get_time(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec/1000000.0;
}

void write_time_file(double i, double m, double u)
{
    FILE *fout;
    fout = fopen("time.txt", "w");
    fprintf(fout, "%e %e %e\n", i, m, u);
    fclose(fout);
}



