#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

int main(void) {
	FILE * fichier;

	fichier = fopen("montest", "w+");
	if (fichier == NULL) {
		printf("oops\n");
		return 1;
	}
	int rand_number = rand()%(214748364 - (-214748364 + 1)) + -214748364;
	//int rand_number = rand()%(1000 - (-1000 + 1)) + -1000;
	for (int i = 0; i < 4096; i++) {
		printf("%d\n", rand_number);
		fprintf(fichier, "%d", rand_number);
		fprintf(fichier, "%s", " ");
		rand_number = rand()%(214748364 - (-214748364 + 1)) + -21410748364;
		//rand_number = rand()%(1000 - (-1000 + 1)) + -1000;
	}
	fclose(fichier);
}