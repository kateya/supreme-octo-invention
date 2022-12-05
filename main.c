// вариант 7 u(x, y, z, t) = f1(x, y, z) * f2(t)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

const int N = 32;
const int K = 1000;

double f1(double l_x, double l_y, double l_z, double x, double y, double z) {
    return sin(2.0 * M_PI * x / l_x + M_PI) * sin(2.0 * M_PI * y / l_y) * sin(M_PI * z / l_z);
}

double f2(double alpha, double t) {
    return cos(alpha * t + M_PI);
}

double laplas(double val, double top, double bottom, double left, double right, double front, double back, double h) {
    return (top + bottom + left + right + front + back - 6.0 * val) / (h * h);
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
    int coord[3];

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

    //printf("%d, %d, %d\n", xdim, ydim, zdim);

    double *analytical1;
    double *grid_1;
    double *grid_2;
    analytical1 = malloc(sizeof(l_x) * xdim * ydim * zdim);
    grid_1 = malloc(sizeof(l_x) * xdim * ydim * zdim);
    grid_2 = malloc(sizeof(l_x) * xdim * ydim * zdim);

    double *top;
    top = malloc(sizeof(l_x) * ydim * zdim);
    double *bottom;
    bottom = malloc(sizeof(l_x) * ydim * zdim);
    double *front;
    front = malloc(sizeof(l_x) * xdim * ydim);
    double *back;
    back = malloc(sizeof(l_x) * xdim * ydim);
    double *left;
    left = malloc(sizeof(l_x) * zdim * xdim);
    double *right;
    right = malloc(sizeof(l_x) * zdim * xdim);

    int xsrc, xdst, ysrc, ydst, zsrc, zdst;

    MPI_Cart_shift(comm, 0, 1, &xsrc, &xdst);
    MPI_Cart_shift(comm, 1, 1, &ysrc, &ydst);
    MPI_Cart_shift(comm, 2, 1, &zsrc, &zdst);
    // u_0 and analytical
    for (int i = 0; i < xdim; i++)
        for (int j = 0; j < ydim; j++)
            for (int k = 0; k < zdim; k++) {
                double val = f1(l_x, l_y, l_z, (i + xdim * coord[0]) * h, (j + ydim * coord[1]) * h, (k + zdim * coord[2]) * h);
                analytical1[i + j * xdim + k * xdim * ydim] = val;
                grid_2[i + j * xdim + k * xdim * ydim] = -val;
                grid_1[i + j * xdim + k * xdim * ydim] = -val;
                if (i == 0 && xsrc != MPI_PROC_NULL) {
                    top[j * zdim + k] = -val;
                }
                else if (i == xdim-1 && xdst != MPI_PROC_NULL) {
                    bottom[j * zdim + k] = -val;
                }
                if (j == 0 && ysrc != MPI_PROC_NULL) {
                    left[i * zdim + k] = -val;
                }
                else if (j == ydim-1 && ydst != MPI_PROC_NULL) {
                    right[i * zdim + k] = -val;
                }
                if (k == 0 && zsrc != MPI_PROC_NULL) {
                    front[i * ydim + j] = -val;
                }
                else if (k == zdim-1 && zdst != MPI_PROC_NULL) {
                    back[i * ydim + j] = -val;
                }
            }

    double t = tau;
    for (int m = 1; m <= 20; t += tau, m++) {
        double root_delta, max_delta = 0.0;
        double *tmp;
        tmp = malloc(sizeof(l_x) * xdim * ydim * zdim);

        MPI_Status status;
        MPI_Request request;
        if (xsrc == MPI_PROC_NULL && xdst != MPI_PROC_NULL) {
            MPI_Isend(bottom, ydim * zdim, MPI_DOUBLE, xdst, 0, comm, &request);
            MPI_Recv(bottom, ydim * zdim, MPI_DOUBLE, xdst, 0, comm, &status);
        } else if (xsrc != MPI_PROC_NULL && xdst == MPI_PROC_NULL) {
            MPI_Isend(top, ydim * zdim, MPI_DOUBLE, xsrc, 0, comm, &request);
            MPI_Recv(top, ydim * zdim, MPI_DOUBLE, xsrc, 0, comm, &status);
        }
        if (ysrc == MPI_PROC_NULL && ydst != MPI_PROC_NULL) {
            MPI_Isend(right, xdim * zdim, MPI_DOUBLE, ydst, 0, comm, &request);
            MPI_Recv(right, xdim * zdim, MPI_DOUBLE, ydst, 0, comm, &status);
        } else if (ysrc != MPI_PROC_NULL && ydst == MPI_PROC_NULL) {
            MPI_Isend(left, xdim * zdim, MPI_DOUBLE, ysrc, 0, comm, &request);
            MPI_Recv(left, xdim * zdim, MPI_DOUBLE, ysrc, 0, comm, &status);
        }
        if (zsrc == MPI_PROC_NULL && zdst != MPI_PROC_NULL) {
            MPI_Isend(front, xdim * ydim, MPI_DOUBLE, zdst, 0, comm, &request);
            MPI_Recv(front, xdim * ydim, MPI_DOUBLE, zdst, 0, comm, &status);
        } else if ((zsrc != MPI_PROC_NULL && zdst == MPI_PROC_NULL)) {
            MPI_Isend(back, xdim * ydim, MPI_DOUBLE, zsrc, 0, comm, &request);
            MPI_Recv(back, xdim * ydim, MPI_DOUBLE, zsrc, 0, comm, &status);
        }

        if (xsrc != MPI_PROC_NULL && xdst != MPI_PROC_NULL) {
            MPI_Sendrecv_replace(top, ydim * zdim, MPI_DOUBLE, xsrc, 0, xdst, 0, comm, &status);
            MPI_Sendrecv_replace(bottom, ydim * zdim, MPI_DOUBLE, xdst, 0, xsrc, 0, comm, &status);
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
            MPI_Sendrecv_replace(left, xdim * zdim, MPI_DOUBLE, zsrc, 0, zdst, 0, comm, &status);
            MPI_Sendrecv_replace(right, xdim * zdim, MPI_DOUBLE, zdst, 0, zsrc, 0, comm, &status);
            double *arr;
            arr = left; left = right; right = arr;
        }

        for (int i = 0; i < xdim; i++)
            for (int j = 0; j < ydim; j++)
                for (int k = 0; k < zdim; k++) {
                    // edge of whole cube
                    if (i == 0 && xsrc == MPI_PROC_NULL)
                        tmp[i + j * xdim + k * xdim * ydim] = analytical1[i + j * xdim + k * xdim * ydim] * f2(alpha, t);
                    else if (i == xdim-1 && xdst == MPI_PROC_NULL)
                        tmp[i + j * xdim + k * xdim * ydim] = analytical1[i + j * xdim + k * xdim * ydim] * f2(alpha, t);
                    else if (j == 0 && ysrc == MPI_PROC_NULL)
                        tmp[i + j * xdim + k * xdim * ydim] = analytical1[i + j * xdim + k * xdim * ydim] * f2(alpha, t);
                    else if (j == ydim-1 && ydst == MPI_PROC_NULL)
                        tmp[i + j * xdim + k * xdim * ydim] = analytical1[i + j * xdim + k * xdim * ydim] * f2(alpha, t);
                    else if (k == 0 && zsrc == MPI_PROC_NULL)
                        tmp[i + j * xdim + k * xdim * ydim] = 0.0;
                    else if (k == zdim-1 && zdst == MPI_PROC_NULL)
                        tmp[i + j * xdim + k * xdim * ydim] = 0.0;
                    else {
                        double top_val, bottom_val, right_val, left_val, front_val, back_val;
                        if (i == 0) {
                            top_val = top[j * zdim + k];
                            bottom_val = grid_1[i+1 + j * xdim + k * xdim * ydim];
                        }
                        else if (i == xdim-1) {
                            bottom_val = bottom[j * zdim + k];
                            top_val = grid_1[i-1 + j * xdim + k * xdim * ydim];
                        } else {
                            top_val = grid_1[i-1 + j * xdim + k * xdim * ydim];
                            bottom_val = grid_1[i+1 + j * xdim + k * xdim * ydim];
                        }
                        if (j == 0) {
                            left_val = left[i * zdim + k];
                            right_val = grid_1[i + (j+1) * xdim + k * xdim * ydim];
                        }
                        else if (j == ydim-1) {
                            right_val = right[i * zdim + k];
                            left_val = grid_1[i + (j-1) * xdim + k * xdim * ydim];
                        } else {
                            right_val = grid_1[i + (j+1) * xdim + k * xdim * ydim];
                            left_val = grid_1[i + (j-1) * xdim + k * xdim * ydim];
                        }
                        if (k == 0) {
                            front_val = front[i * ydim + j];
                            back_val = grid_1[i + j * xdim + (k+1) * xdim * ydim];
                        }
                        else if (k == zdim-1) {
                            back_val = back[i * ydim + j];
                            front_val = grid_1[i + j * xdim + (k-1) * xdim * ydim];
                        } else {
                            back_val = grid_1[i + j * xdim + (k+1) * xdim * ydim];
                            front_val = grid_1[i + j * xdim + (k-1) * xdim * ydim];
                        }
                        double res = 2.0 * grid_1[i + j * xdim + k * xdim * ydim] - grid_2[i + j * xdim + k * xdim * ydim] + \
                                     tau * tau * laplas(grid_1[i + j * xdim + k * xdim * ydim], top_val, bottom_val, left_val, right_val, front_val, back_val, h);
                        tmp[i + j * xdim + k * xdim * ydim] = res;
                        double delta = fabs(res - analytical1[i + j * xdim + k * xdim * ydim] * f2(alpha, t));
                        if (delta > max_delta)
                            max_delta = delta;

                        if (i == 0) {
                            top[j * zdim + k] = res;
                        }
                        else if (i == xdim-1) {
                            bottom[j * zdim + k] = res;
                        }
                        if (j == 0) {
                            left[i * zdim + k] = res;
                        }
                        else if (j == ydim-1) {
                            right[i * zdim + k] = res;
                        }
                        if (k == 0) {
                            front[i * ydim + j] = res;
                        }
                        else if (k == zdim-1) {
                            back[i * ydim + j] = res;
                        }
                    }
                }
        MPI_Reduce(&max_delta, &root_delta, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        if (rank == 0)
            printf("t = %d, eps = %f\n", m, root_delta);
        free(grid_2);
        grid_2 = grid_1;
        grid_1 = tmp;
    }
    free(grid_2);
    free(grid_1);
    free(analytical1);
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