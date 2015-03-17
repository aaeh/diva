/* rand example: guess the number */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

int main (int argc,char *argv[])
{
  float r1,r2,r3,r4, r5, min , max;
  int rounds; 

  /* initialize random seed: */
  srand (time(NULL));

 min = atof(argv[1]);
 max = atof(argv[2]);
 rounds =  atof(argv[3]);

 printf("%f \n",  min);
 printf("%f \n",   max);
 printf(" %d \n", rounds);



  /* generate secret number between 1 and 10: */
 for(int i=0; i< rounds; i++){
 r1 = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
 r2 = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
 r3 = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
 r4 = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
 r5 = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;

  
  
  printf("%d , %f \n", i,  r1+r2+r3+r4+r5);
}
  return 0;
}
