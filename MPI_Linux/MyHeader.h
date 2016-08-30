void val_intialize(int *local_val, int local_nnz, int my_rank, int comm_sz, MPI_Comm comm);
void col_intialize(int *local_col, int n, int local_n, int local_nnz, int my_rank, int comm_sz, MPI_Comm comm);
void row_ptr_intialize(int *local_row_ptr, int local_n, int my_rank, int comm_sz, MPI_Comm comm);
void calc_matrix(int *local_val, int *local_col, int *local_row_ptr, int n, int p, int local_n, int my_rank, int comm_sz, MPI_Comm comm);
void write_x_file(double *x, int n, int z);
double get_time(void);
void write_time_file(double i, double m, double u);