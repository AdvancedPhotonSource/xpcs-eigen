/*
* @Author: Faisal Khan
* @Date:   2017-11-17 10:55:13
* @Last Modified by:   Faisal Khan
* @Last Modified time: 2017-11-17 13:33:19
*/
#include <stdio.h>
#include <omp.h>

static long num_steps = 1000000;
double step;

int main()
{
    double pi;
    step = 1.0/(double) num_steps;
    double *t_sum = 0;
    int steps_per_threads = 1;
    
    double t0 = omp_get_wtime();
    int threads = 0;

    #pragma omp parallel
    #pragma omp master
    {
        threads = omp_get_num_threads();
        steps_per_threads = num_steps / threads;
        t_sum = new double[threads];
        printf("Threads = %d, Steps per Threads = %d\n", threads, steps_per_threads);
    }

    #pragma omp parallel
    {
        
        int i = 0;
        double x = 0;
        int id = omp_get_thread_num();
        t_sum[id] = 0.0;

        for (i = id * steps_per_threads; i < (id*steps_per_threads + steps_per_threads); i++) 
        {
            x = (i + 0.5) * step;
            t_sum[id] = t_sum[id] + 4.0/(1.0+x*x);
        }
    }

    double t1 = omp_get_wtime();

    double sum = 0.0;
    int i;
    for (i = 0 ; i < threads; i++)
        sum += t_sum[i];

    pi = step * sum;

    printf("Compute %f in %f time\n", pi, (t1-t0));
}