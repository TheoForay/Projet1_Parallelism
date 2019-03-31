long max(long a, long b)																		: non parrelele

void printArray(struct tablo * tmp)																: non parrelele

struct tablo * allocateTablo(long size)															: non parrelele

void sum_prefix_montee(struct tablo * source, struct tablo * destination)						: parrelele

void sum_prefix_descente(struct tablo * a, struct tablo * b)									: parrelele

void sum_prefix_final(struct tablo * a, struct tablo *b)										: parrelele

void generateArray(struct tablo * s)															: non parrelele

void max_prefix_montee(struct tablo * source, struct tablo * destination)						: parrelele

void max_prefix_descente(struct tablo * a, struct tablo * b)									: parrelele

void max_prefix_final(struct tablo * a, struct tablo *b)										: parrelele

void sum_suffix_descente(struct tablo * a, struct tablo * b)									: parrelele

void max_suffix_descente(struct tablo * a, struct tablo * b)									: parrelele

void getGoodArray(struct tablo * final, struct tablo * apresFinal)								: non parrelele

void sum_prefix(struct tablo source, struct tablo * goodValues)									: non parrelele

void sum_suffix(struct tablo source, struct tablo * goodValues)									: non parrelele

void max_prefix(struct tablo source, struct tablo * goodValues)									: non parrelele

void max_suffix(struct tablo source, struct tablo * goodValues)									: non parrelele

void getMaxParallel(struct tablo * source, struct tablo * tabMax)								: parrelele

void getMaxSubArrayIndices(struct tablo * source, struct tablo * tabIndices, long maximum)		: non parrelele

void printMaxSubArray(struct tablo source, struct tablo * tabIndices) 							: non parrelele

void negativeInit(struct tablo * tabIndices)													: non parrelele

void generateArrayFromFile(char* filename, FILE* fichier, struct tablo * source)				: non parrelele

int main(int argc, char **argv)																	: parrelele