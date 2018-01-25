/*
* @Author: Faisal Khan
* @Date:   2018-01-18 14:30:19
* @Last Modified by:   Faisal Khan
* @Last Modified time: 2018-01-22 18:21:51
*/

#include <omp.h>
#include <stdio.h>
#include <assert.h>


#define NUM_THREADS 2

struct Node {
  int value;
  struct Node* next = NULL;
};

struct Queue {
  struct Node* head = NULL;
  struct Node* tail = NULL;
  int enqueued = 0;
  int dequeued = 0;
};

void Enqueue(struct Queue* q, int value)
{
  struct Node* n = new Node;
  n->value = value;
  
  if (q->head == NULL)
    q->head = n;
  else
    q->tail->next = n;

  q->tail = n;
  q->enqueued++;
}

struct Node* Dequeue(struct Queue *q)
{
  assert(q->enqueued != 0);

  struct Node* res = q->head;
  struct Node* next = q->head->next;
  q->head = next;
  q->dequeued++;

  return res;
}

int main(int argc, char** argv)
{
  omp_set_num_threads(NUM_THREADS);

  struct Queue q;

  bool isDone = false;

  int queue_size = 0;

  #pragma omp parallel
  {
    #pragma omp sections nowait
    {
      #pragma omp section
      {
        for (int i = 0; i < 10000000; i++)
          Enqueue(&q, i+1);
        isDone = true;
      }

      #pragma omp section
      while (1)
      {
        queue_size = q.enqueued - q.dequeued;

        if (queue_size == 0 && isDone) 
          break;
        
        if (queue_size == 0)
          continue;

        struct Node* n = NULL;
        if (queue_size == 1)
        {
          #pragma omp critical
          n = Dequeue(&q);
        }
        else
          n = Dequeue(&q);

        // printf("%d\n", n->value);
        delete(n);
        
      }

    } 
  }

}

