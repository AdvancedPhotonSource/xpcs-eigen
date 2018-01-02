/*
* @Author: Faisal Khan
* @Date:   2017-11-17 10:55:13
* @Last Modified by:   Faisal Khan
* @Last Modified time: 2017-11-17 13:34:24
*/
#include <stdio.h>
#include <omp.h>
#include <chrono>

static long num_steps = 1000000;
double step;

#define NUM_THREADS 10
int main()
{
    int i, nthreads; double pi, sum = 0.0;
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);

    auto t0 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
        int i, id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        if (id == 0) nthreads = nthrds;

        double t_sum = 0;
        for (i = id; i < num_steps; i=i+nthrds)
        {
            x = (i + 0.5) * step;
            t_sum += 4.0/(1.0 + x * x);
        }

        #pragma omp critical
            sum += t_sum;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;

    pi = sum * step;

    printf("Compute %f in %f secs\n", pi, (elapsed.count()));
}