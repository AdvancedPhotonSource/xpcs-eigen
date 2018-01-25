/*
* @Author: Faisal Khan
* @Date:   2018-01-18 14:30:19
* @Last Modified by:   Faisal Khan
* @Last Modified time: 2018-01-23 11:46:57
*/

#include <omp.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#define NUM_THREADS 1

struct Task {
  int id;
};

void do_work(struct Task* t)
{
  printf("Working on task # %d\n", t->id);
  usleep(500 * 1000);
  printf("Done with task  # %d\n", t->id);
}


int main(int argc, char** argv)
{
  #pragma omp parallel
  {
    printf("total number of thread %d\n", omp_get_num_threads());
    #pragma omp single
    {
      int total_tasks = 10;
      struct Task** tasks = new Task*[total_tasks];

      for (int i = 0; i < total_tasks; i++)
      {
        tasks[i] = new Task;
        tasks[i]->id = (i+1);

        struct Task* t = tasks[i];
        #pragma omp task firstprivate(t)
          do_work(t);
      }
    }

  }
}

