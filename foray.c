#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>

//Partie réalisée pour le TP, partie fournie
struct tablo {
  long * tab;
  long size;
};

long max(long a, long b) {
  if (a >= b) {
    return a;
  }
  return b;
}

void printArray(struct tablo * tmp) {
  printf("---- Array of size %ld ---- \n", tmp->size);
  long size = tmp->size;
  int i;
  for (i = 0; i < size; ++i) {
    printf("%ld ", tmp->tab[i]);
  }
  printf("\n");
}

struct tablo * allocateTablo(long size) {
  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->size = size;
  tmp->tab = malloc(size*sizeof(long));
  tmp->tab[0]=0;
  return tmp;
}

void sum_prefix_montee(struct tablo * source, struct tablo * destination) {
  // TODO : fill the destination array of size 2*n by copying the 
  // source array at the end
  // You can assume the malloc for the destination array has been performed.
  long a_size = source->size;
  long b_size = destination->size;
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

//Fonction non-utilisée pour le projet puisqu'on lit le tableau dans un fichier
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

//Partie réalisée pour le projet :

//Les trois fonctions suivantes servent à créer le tableau max en prefix
void max_prefix_montee(struct tablo * source, struct tablo * destination) {
  long a_size = source->size;
  long b_size = destination->size;
  for (; a_size > 0; ) {
    destination->tab[b_size-1] = source->tab[a_size-1];
    a_size--;
    b_size--;
  }

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
  b->tab[1] = INT_MIN;
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
  int m = log2((a->size+1)/2);
  int begin = pow(2.0,m);
  int end = pow(2.0,m+1)-1;
  #pragma omp parallel for
  for (int j = begin; j <= end; j++) {
    b->tab[j] = max(b->tab[j], a->tab[j]);
  }
}

//Les deux fonctions suivantes font "suffix descente" avec des operations differentes (sum et max). Elles suivent ce procédé:
// Racine élément neutre
  // La différence est donc : 
  //  si le noeud est le fils droit de son père, il prend la valeur de son père
  //  si le noeud est le fils gauche de son père, il prend la valeur de son père (*) valeur de son frère dans l'arbre de la montée
  //    -> C'est pourquoi la condition est "inversée" et l'indice pour a->tab devient j+1 au lieu de j-1 car puisque c'est le miroir, 
  //      le frère est de l'autre coté par rapport au père

//Cette phase de descente pour le suffix marche est inspirée par la descente de prefix
void sum_suffix_descente(struct tablo * a, struct tablo * b) {
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

//Cette phase de descente pour le suffix marche est inspirée par la descente de prefix
void max_suffix_descente(struct tablo * a, struct tablo * b) {
  b->tab[1] = INT_MIN;
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

//Fonction qui récupère la partie du tableau qui contient le resultat des fonctions précédentes
void getGoodArray(struct tablo * final, struct tablo * apresFinal) {
  int j = apresFinal->size/2;
  for (int i = 0; i < final->size; i++) {
    final->tab[i] = apresFinal->tab[j+i];
  }
}

//Les 4 fonctions suivantes appliques les algorithmes de prefix et suffix avec les operations de sum et max
void sum_prefix(struct tablo source, struct tablo * goodValues) {
  struct tablo * a = allocateTablo(source.size*2);
  sum_prefix_montee(&source, a);

  struct tablo * b = allocateTablo(source.size*2);
  sum_prefix_descente(a, b);
   
  sum_prefix_final(a,b);

  getGoodArray(goodValues, b);
}
void sum_suffix(struct tablo source, struct tablo * goodValues) {
  struct tablo * a = allocateTablo(source.size*2);
  sum_prefix_montee(&source, a);

  struct tablo * b = allocateTablo(source.size*2);
  sum_suffix_descente(a, b);

  sum_prefix_final(a,b);

  getGoodArray(goodValues, b);
}

void max_prefix(struct tablo source, struct tablo * goodValues) {
  struct tablo * a = allocateTablo(source.size*2);
  max_prefix_montee(&source, a);

  struct tablo * b = allocateTablo(source.size*2);
  max_prefix_descente(a, b);
   
  max_prefix_final(a,b);

  getGoodArray(goodValues, b);
}
void max_suffix(struct tablo source, struct tablo * goodValues) {
  struct tablo * a = allocateTablo(source.size*2);
  max_prefix_montee(&source, a);

  struct tablo * b = allocateTablo(source.size*2);
  max_suffix_descente(a, b);
   
  max_prefix_final(a,b);

  getGoodArray(goodValues, b);
}

//Fonction du cours (équivalente à Montée) pour récupérer le max de façon parallele
void getMaxParallel(struct tablo * source, struct tablo * tabMax) {
	//On copie source à la fin de tabMax
	long a_size = source->size;
	long b_size = tabMax->size;
	for (; a_size > 0; ) {
		tabMax->tab[b_size-1] = source->tab[a_size-1];
		a_size--;
		b_size--;
	}

	int m = log2((source->size+1));
	for (int l = m-1; l >= 0; l--) {
    int borne = (pow(2.0, l+1)-1);
    #pragma omp parallel for
		for (int j = pow(2.0, l); j <= borne; j++) {
			tabMax->tab[j] = max(tabMax->tab[2*j], tabMax->tab[2*j+1]);
		}
	}
}

//Fonction qui récupère les indices du tableau initial où se trouvent les valeurs qui forment le max subarray
void getMaxSubArrayIndices(struct tablo * source, struct tablo * tabIndices, long maximum) {
	int j = 0;

	for(int i = 0; i < source->size; i++) {
		//Puisque les sum max sont tous regroupés dans M (voir main), pas besoin vérifier autre chose
		if (source->tab[i] == maximum) {
			tabIndices->tab[j] = i;
			j++;
		}
	}

}

//Fonction qui print les valeurs correspondantes au max subarray
void printMaxSubArray(struct tablo source, struct tablo * tabIndices) {
	for (int i = 0; i < tabIndices->size; i++) {
		if (tabIndices->tab[i] >= 0) {
			printf("%ld ", source.tab[tabIndices->tab[i]]);
		}
	}
	printf("\n");
}

void negativeInit(struct tablo * tabIndices) {
	for (int i = 0; i < tabIndices->size; i++) {
		tabIndices->tab[i] = -1;
	}
}

//Fonction qui crée le tableau initial à partir d'un fichier passé en argument
void generateArrayFromFile(char* filename, FILE* fichier, struct tablo * source) {
  fichier = fopen(filename, "r");
  int i = 0;
  long c = 0;
  int potentialEOF;
  if (fichier != NULL) {
    source->tab = malloc(100*sizeof(long));
    potentialEOF = fscanf(fichier, "%ld", &c);
    while (potentialEOF != EOF) {
      source->tab[i] = c;
      i++;
      potentialEOF = fscanf(fichier, "%ld", &c);
      if (i%100 == 0) {
        source->tab = realloc(source->tab, (i+100)*sizeof(long));
      }
    }
    source->size = i;
    fclose(fichier);
  } else {
    printf("Impossible de lire le fichier\n");
  }
}

int main(int argc, char **argv) {
  struct tablo source;
  FILE* fichier = NULL;

  //generateArray(&source);
  
  generateArrayFromFile(argv[1], fichier, &source);
  
  struct tablo * psum = allocateTablo(source.size);
  struct tablo * ssum = allocateTablo(source.size);
  struct tablo * smax = allocateTablo(source.size);
  struct tablo * pmax = allocateTablo(source.size);

  sum_prefix(source, psum);
  sum_suffix(source, ssum);
  max_suffix(*psum, smax);
  max_prefix(*ssum, pmax);

  struct tablo * M = allocateTablo(source.size);
  #pragma omp parallel for
  for (int i = 0; i < source.size; i++) {
    long Ms = pmax->tab[i] - ssum->tab[i] + source.tab[i];
    long Mp = smax->tab[i] - psum->tab[i] + source.tab[i];
    M->tab[i] = Ms + Mp - source.tab[i];
  }
  
  struct tablo * tabMax = allocateTablo(M->size*2);
  getMaxParallel(M, tabMax);

  struct tablo * tabIndices = allocateTablo(M->size);
  negativeInit(tabIndices);
  getMaxSubArrayIndices(M, tabIndices, tabMax->tab[1]);
  printf("%ld ", tabMax->tab[1]);
  printMaxSubArray(source, tabIndices);
}