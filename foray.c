#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

struct tablo {
  int * tab;
  int size;
};

int max(int a, int b) {
  if (a >= b) {
    return a;
  }
  return b;
}

void printArray(struct tablo * tmp) {
  printf("---- Array of size %i ---- \n", tmp->size);
  int size = tmp->size;
  int i;
  for (i = 0; i < size; ++i) {
    printf("%i ", tmp->tab[i]);
  }
  printf("\n");
}

struct tablo * allocateTablo(int size) {
  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->size = size;
  tmp->tab = malloc(size*sizeof(int));
  tmp->tab[0]=0;
  return tmp;
}

void sum_prefix_montee(struct tablo * source, struct tablo * destination) {
  //TODO : fill the destination array of size 2*n by copying the 
  // source array at the end
  // You can assume the malloc for the destination array has been performed.
  int a_size = source->size;
  int b_size = destination->size;
  for (; a_size > 0; ) {
    destination->tab[b_size-1] = source->tab[a_size-1];
    a_size--;
    b_size--;
  }

  // TODO: Implement the up phase of the algorithm
  int m = log2(source->size);
  for (int l = m-1; l >= 0; l--) {
    int begin = pow(2.0,l);
    int end = pow(2.0,l+1)-1;
    #pragma omp parallel for
    for (int j = begin; j <= end; j++) {
      destination->tab[j] = destination->tab[2*j] + destination->tab[2*j+1];
    }
  }
}

void sum_prefix_descente(struct tablo * a, struct tablo * b) {
  // TODO : implement the down phase of the algorithm
  b->tab[1] = 0;
  int m = log2((a->size+1)/2);
  for (int l = 1; l <= m; l++) {
    int begin = pow(2.0, l);
    int end = pow(2.0,l+1) - 1;
    #pragma omp parallel for
    for (int j = begin; j <= end; j++) {
      if (j%2 == 0) {
        b->tab[j] = b->tab[j/2];
      } else {
        b->tab[j] = b->tab[j/2] + a->tab[j-1];
      }
    }
  }
}

void sum_prefix_final(struct tablo * a, struct tablo *b) {
  // TODO : implement the final phase
  int m = log2((a->size+1)/2);
  int begin = pow(2.0,m);
  int end = pow(2.0,m+1)-1;
  #pragma omp parallel for
  for (int j = begin; j <= end; j++) {
    b->tab[j] = b->tab[j] + a->tab[j];
  }
}

void max_prefix_montee(struct tablo * source, struct tablo * destination) {
  //TODO : fill the destination array of size 2*n by copying the 
  // source array at the end
  // You can assume the malloc for the destination array has been performed.
  int a_size = source->size;
  int b_size = destination->size;
  for (; a_size > 0; ) {
    destination->tab[b_size-1] = source->tab[a_size-1];
    a_size--;
    b_size--;
  }

  // TODO: Implement the up phase of the algorithm
  int m = log2(source->size);
  for (int l = m-1; l >= 0; l--) {
    int begin = pow(2.0,l);
    int end = pow(2.0,l+1)-1;
    #pragma omp parallel for
    for (int j = begin; j <= end; j++) {
      destination->tab[j] = max(destination->tab[2*j], destination->tab[2*j+1]);
    }
  }
}

void max_prefix_descente(struct tablo * a, struct tablo * b) {
  // TODO : implement the down phase of the algorithm
  b->tab[1] = -100000;
  int m = log2((a->size+1)/2);
  for (int l = 1; l <= m; l++) {
    int begin = pow(2.0, l);
    int end = pow(2.0,l+1) - 1;
    #pragma omp parallel for
    for (int j = begin; j <= end; j++) {
      if (j%2 == 0) {
        b->tab[j] = b->tab[j/2];
      } else {
        b->tab[j] = max(b->tab[j/2], a->tab[j-1]);
      }
    }
  }
}

void max_prefix_final(struct tablo * a, struct tablo *b) {
  // TODO : implement the final phase
  int m = log2((a->size+1)/2);
  int begin = pow(2.0,m);
  int end = pow(2.0,m+1)-1;
  #pragma omp parallel for
  for (int j = begin; j <= end; j++) {
    b->tab[j] = max(b->tab[j], a->tab[j]);
  }
}

void sum_suffix_descente(struct tablo * a, struct tablo * b) {
  // Cette phase de descente pour le suffix marche est inspirée par la descente de prefix
  // Racine élément neutre
  // La différence est donc : 
  //  si le noeud est le fils droit de son père, il prend la valeur de son père
  //  si le noeud est le fils gauche de son père, il prend la valeur de son père (*) valeur de son frère dans l'arbre de la montée
  //    -> C'est pourquoi la condition est "inversée" et l'indice pour a->tab devient j+1 au lieu de j-1 car puisque c'est le miroir, 
  //      le frère est de l'autre coté par rapport au père

  b->tab[1] = 0;
  int m = log2((a->size+1)/2);
  for (int l = 1; l <= m; l++) {
    int begin = pow(2.0, l);
    int end = pow(2.0,l+1) - 1;
    #pragma omp parallel for
    for (int j = begin; j <= end; j++) {
      if (j%2 == 0) {
        b->tab[j] = b->tab[j/2] + a->tab[j+1];
      } else {
        b->tab[j] = b->tab[j/2];
      }
    }
  }
}

void max_suffix_descente(struct tablo * a, struct tablo * b) {
  // Cette phase de descente pour le suffix marche est inspirée par la descente de prefix
  // Racine élément neutre
  // La différence est donc : 
  //  si le noeud est le fils droit de son père, il prend la valeur de son père
  //  si le noeud est le fils gauche de son père, il prend la valeur de son père (*) valeur de son frère dans l'arbre de la montée
  //    -> C'est pourquoi la condition est "inversée" et l'indice pour a->tab devient j+1 au lieu de j-1 car puisque c'est le miroir, 
  //      le frère est de l'autre coté par rapport au père

  b->tab[1] = -100000;
  int m = log2((a->size+1)/2);
  for (int l = 1; l <= m; l++) {
    int begin = pow(2.0, l);
    int end = pow(2.0,l+1) - 1;
    #pragma omp parallel for
    for (int j = begin; j <= end; j++) {
      if (j%2 == 0) {
        b->tab[j] = max(b->tab[j/2], a->tab[j+1]);
      } else {
        b->tab[j] = b->tab[j/2];
      }
    }
  }
}

void generateArray(struct tablo * s) {
  s->size=16;
  s->tab=malloc(s->size*sizeof(int));
  s->tab[0]=3;
  s->tab[1]=2;
  s->tab[2]=-7;
  s->tab[3]=11;
  s->tab[4]=10;
  s->tab[5]=-6;
  s->tab[6]=4;
  s->tab[7]=9;
  s->tab[8]=-6;
  s->tab[9]=1;
  s->tab[10]=-2;
  s->tab[11]=-3;
  s->tab[12]=4;
  s->tab[13]=-3;
  s->tab[14]=0;
  s->tab[15]=2;
}

void getGoodArray(struct tablo * final, struct tablo * apresFinal) {
  int j = apresFinal->size/2;
  for (int i = 0; i < final->size; i++) {
    final->tab[i] = apresFinal->tab[j+i];
  }
}

void sum_prefix(struct tablo source, struct tablo * goodValues) {
  //printf("Montee : \n");
  struct tablo * a = allocateTablo(source.size*2);
  sum_prefix_montee(&source, a);
  /*printArray(a);
  printf("__________________________________\n");*/

  //printf("Descente : \n");
  struct tablo * b = allocateTablo(source.size*2);
  sum_prefix_descente(a, b);
  /*printArray(b);
  printf("__________________________________\n");*/
   
  //printf("Tableau final total: \n");
  sum_prefix_final(a,b);
  //printArray(b);

  getGoodArray(goodValues, b);
}
void sum_suffix(struct tablo source, struct tablo * goodValues) {
  //printf("Montee : \n");
  struct tablo * a = allocateTablo(source.size*2);
  sum_prefix_montee(&source, a);
  /*printArray(a);
  printf("__________________________________\n");*/

  //printf("Descente : \n");
  struct tablo * b = allocateTablo(source.size*2);
  sum_suffix_descente(a, b);
  /*printArray(b);
  printf("__________________________________\n");*/
   
  //printf("Tableau final total: \n");
  sum_prefix_final(a,b);
  //printArray(b);

  getGoodArray(goodValues, b);
}

void max_prefix(struct tablo source, struct tablo * goodValues) {
  //printf("Montee : \n");
  struct tablo * a = allocateTablo(source.size*2);
  max_prefix_montee(&source, a);
  /*printArray(a);
  printf("__________________________________\n");*/

 //printf("Descente : \n");
  struct tablo * b = allocateTablo(source.size*2);
  max_prefix_descente(a, b);
  /*printArray(b);
  printf("__________________________________\n");*/
   
  //printf("Tableau final total: \n");
  max_prefix_final(a,b);
  //printArray(b);

  getGoodArray(goodValues, b);
}

void max_suffix(struct tablo source, struct tablo * goodValues) {
  //printf("Montee : \n");
  struct tablo * a = allocateTablo(source.size*2);
  max_prefix_montee(&source, a);
  /*printArray(a);
  printf("__________________________________\n");*/

  //printf("Descente : \n");
  struct tablo * b = allocateTablo(source.size*2);
  max_suffix_descente(a, b);
  /*printArray(b);
  printf("__________________________________\n");*/
   
  //printf("Tableau final total: \n");
  max_prefix_final(a,b);
  //printArray(b);

  getGoodArray(goodValues, b);
  /*printf("Tableau final max_suffix: \n");
  printArray(goodValues);*/
}

int main(int argc, char **argv) {
  struct tablo source;

  generateArray(&source);
  /*printf("Tableau source : \n");
  printArray(&source);*/

  struct tablo * psum = allocateTablo(source.size);
  struct tablo * ssum = allocateTablo(source.size);
  struct tablo * smax = allocateTablo(source.size);
  struct tablo * pmax = allocateTablo(source.size);

  /*printf("############################sum_prefix#######################\n");
  sum_prefix(source, psum);
  printf("Tableau final sum_prefix: \n");
  printArray(psum);

  printf("\n############################sum_suffix#######################\n"); 
  sum_suffix(source, ssum);
  printf("Tableau final sum_suffix: \n");
  printArray(ssum);

  printf("\n############################max_suffix#######################\n"); 
  max_suffix(*psum, smax);
  printf("Tableau final max_suffix: \n");
  printArray(smax);

  printf("\n############################max_prefix#######################\n"); 
  max_prefix(*ssum, pmax); 
  printf("Tableau final max_prefix\n");
  printArray(pmax);*/

  sum_prefix(source, psum);
  sum_suffix(source, ssum);
  max_suffix(*psum, smax);
  max_prefix(*ssum, pmax);

  printf("\n############################    M    #######################\n"); 
  struct tablo * M = allocateTablo(source.size);
  for (int i = 0; i < source.size; i++) {
    int Ms = pmax->tab[i] - ssum->tab[i] + source.tab[i];
    int Mp = smax->tab[i] - psum->tab[i] + source.tab[i];
    M->tab[i] = Ms + Mp - source.tab[i];
  }
  printArray(M);
  
}