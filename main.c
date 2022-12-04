// вариант 7 u(x, y, z, t) = f1(x, y, z) * f2(t)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

const int N = 129;
const int K = 1000;

double f1(double l_x, double l_y, double l_z, double x, double y, double z) {
    return sin(2.0 * M_PI * x / l_x + M_PI) * sin(2.0 * M_PI * y / l_y) * sin(M_PI * z / l_z);
}

double f2(double alpha, double t) {
    return cos(alpha * t + M_PI);
}

double laplas(double *arr, int i, int j, int k, double h) {
    return (arr[i-1 + j * N + k * N * N] + arr[i+1 + j * N + k * N * N] +
           arr[i + (j-1) * N + k * N * N] + arr[i + (j+1) * N + k * N * N] +
           arr[i + j * N + (k-1) * N * N] + arr[i + j * N + (k+1) * N * N] - 6.0 * arr[i + j * ydim + k * zdim * zdim]) / (h * h);
}

int main(int argc, char* argv[]) {
    char *str_end;
    double l_x = strtod(argv[1], &str_end);
    double l_y = strtod(argv[2], &str_end);
    double l_z = strtod(argv[3], &str_end);
    double h = l_x / (N - 1);
    double tau = 1.0 / K;
    double alpha = M_PI * sqrt(4.0 / (l_x * l_x) + 4.0 / (l_y * l_y) + 1.0 / (l_z * l_z));


    int rank, size;
    MPI_Comm comm;
    int dim[3], period[3], reorder;
    int coord[3], id;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    period[0] = 0; period[1] = 0; period[2] = 0;
    dim[0] = 0; dim[1] = 0; dim[2] = 0;
    reorder = 0;
    MPI_Dims_create(size, 3, dim);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, reorder, &comm);
    MPI_Cart_coords(comm, rank, 3, coord);

    int xdim = N / dim[0];
    int ydim = N / dim[1];
    int zdim = N / dim[2];

    double *analytical1;
    double *grid_1;
    double *grid_2;
    analytical1 = malloc(sizeof(l_x) * xdim * ydim * zdim);
    grid_1 = malloc(sizeof(l_x) * xdim * ydim * zdim);
    grid_2 = malloc(sizeof(l_x) * xdim * ydim * zdim);
    // u_0 and analytical
    for (int i = 0; i < xdim; i++)
        for (int j = 0; j < ydim; j++)
            for (int k = 0; k < zdim; k++) {
                double val = f1(l_x, l_y, l_z, (i + xdim * coord[0]) * h, (j + ydim * coord[1]) * h, (k + xdim * coord[2]) * h);
                analytical1[i + j * ydim + k * zdim * zdim] = val;
                grid_2[i + j * ydim + k * zdim * zdim] = -val;
                grid_1[i + j * ydim + k * zdim * zdim] = -val;
            }

    double *top;
    top = malloc(sizeof(l_x) * xdim * zdim);
    double *bottom;
    bottom = malloc(sizeof(l_x) * xdim * zdim);
    double *front;
    front = malloc(sizeof(l_x) * xdim * ydim);
    double *back;
    back = malloc(sizeof(l_x) * xdim * ydim);
    double *left;
    left = malloc(sizeof(l_x) * zdim * ydim);
    double *right;
    right = malloc(sizeof(l_x) * zdim * ydim);

    int xsrc, xdst, ysrc, ydst, zsrc, zdst;

    MPI_Cart_shift(comm, 0, 1, &xsrc, &xdst);
    MPI_Cart_shift(comm, 1, 1, &ysrc, &ydst);
    MPI_Cart_shift(comm, 2, 1, &zsrc, &zdst);

    double t = tau;
    for (int m = 1; m <= K; t += tau, m++) {
        double max_delta = 0.0;
        double *tmp;
        tmp = malloc(sizeof(l_x) * xdim * ydim * zdim);

        MPI_Status status;
        if (dim[0] == 2) {
            if (coord[0] == 0) {
                MPI_Send(right, ydim * zdim, MPI_DOUBLE, xdst, 0, comm);
                MPI_Recv(right, ydim * zdim, MPI_DOUBLE, xdst, 0, comm, &status);
            } else {
                MPI_Send(left, ydim * zdim, MPI_DOUBLE, xsrc, 0, comm);
                MPI_Recv(left, ydim * zdim, MPI_DOUBLE, xsrc, 0, comm, &status);
            }
        }
        if (dim[1] == 2) {
            if (coord[1] == 0) {
                MPI_Send(bottom, xdim * zdim, MPI_DOUBLE, ydst, 0, comm);
                MPI_Recv(bottom, xdim * zdim, MPI_DOUBLE, ydst, 0, comm, &status);
            } else {
                MPI_Send(top, xdim * zdim, MPI_DOUBLE, ysrc, 0, comm);
                MPI_Recv(top, xdim * zdim, MPI_DOUBLE, ysrc, 0, comm, &status);
            }
        }
        if (dim[2] == 2) {
            if (coord[2] == 0) {
                MPI_Send(front, xdim * ydim, MPI_DOUBLE, zdst, 0, comm);
                MPI_Recv(front, xdim * ydim, MPI_DOUBLE, zdst, 0, comm, &status);
            } else {
                MPI_Send(back, xdim * ydim, MPI_DOUBLE, zsrc, 0, comm);
                MPI_Recv(back, xdim * ydim, MPI_DOUBLE, zsrc, 0, comm, &status);
            }
        }
        if (xsrc != MPI_PROC_NULL && xdst != MPI_PROC_NULL) {
            MPI_Sendrecv_replace(top, xdim * zdim, MPI_DOUBLE, xsrc, 0, xdst, 0, comm, &status);
            MPI_Sendrecv_replace(bottom, xdim * zdim, MPI_DOUBLE, xdst, 0, xsrc, 0, comm, &status);
            double *arr;
            arr = top; top = bottom; bottom = arr;
        }
        if (ysrc != MPI_PROC_NULL && ydst != MPI_PROC_NULL) {
            MPI_Sendrecv_replace(front, xdim * ydim, MPI_DOUBLE, ysrc, 0, ydst, 0, comm, &status);
            MPI_Sendrecv_replace(back, xdim * ydim, MPI_DOUBLE, ydst, 0, ysrc, 0, comm, &status);
            double *arr;
            arr = front; front = back; back = arr;
        }
        if (zsrc != MPI_PROC_NULL && zdst != MPI_PROC_NULL) {
            MPI_Sendrecv_replace(left, ydim * zdim, MPI_DOUBLE, zsrc, 0, zdst, 0, comm, &status);
            MPI_Sendrecv_replace(right, ydim * zdim, MPI_DOUBLE, zdst, 0, zsrc, 0, comm, &status);
            double *arr;
            arr = left; left = right; right = arr;
        }

        for (int i = 0; i < xdim; i++)
            for (int j = 0; j < ydim; j++)
                for (int k = 0; k < zdim; k++) {
                    // edge of cube
                    if (i == 0 && coord[0] == 0)
                        tmp[i + j * ydim + k * zdim * zdim] = analytical1[i + j * ydim + k * zdim * zdim] * f2(alpha, t);
                    else if (i == N-1 && coord[0] == dim[0] - 1)
                        tmp[i + j * ydim + k * zdim * zdim] = analytical1[i + j * ydim + k * zdim * zdim] * f2(alpha, t);
                    else if (j == 0 && coord[1] == 0)
                        tmp[i + j * ydim + k * zdim * zdim] = analytical1[i + j * ydim + k * zdim * zdim] * f2(alpha, t);
                    else if (j == N-1 && coord[1] == dim[1] - 1)
                        tmp[i + j * ydim + k * zdim * zdim] = analytical1[i + j * ydim + k * zdim * zdim] * f2(alpha, t);
                    else if (k == 0 && coord[2] == 0)
                        tmp[i + j * ydim + k * zdim * zdim] = 0.0;
                    else if (k == N-1 && coord[2] == dim[2] - 1)
                        tmp[i + j * ydim + k * zdim * zdim] = 0.0;
                    else {
                        double res = 2.0 * grid_1[i + j * ydim + k * zdim * zdim] - grid_2[i + j * ydim + k * zdim * zdim] + tau * tau * laplas(grid_1, i, j, k, h);
                        tmp[i + j * ydim + k * zdim * zdim] = res;
                        double delta = fabs(res - analytical1[i + j * ydim + k * zdim * zdim] * f2(alpha, t));
                        if (delta > max_delta)
                            max_delta = delta;
                    }
                }
        printf("%f\n", max_delta);
        free(grid_2);
        grid_2 = grid_1;
        grid_1 = tmp;
    }
    free(grid_2);
    free(grid_1);
    free(top);
    free(bottom);
    free(front);
    free(back);
    free(left);
    free(right);
    MPI_Comm_free( &comm );

    MPI_Finalize();
    return 0;
}