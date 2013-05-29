#include "utils.h"

char *remove_ext(const char* mystr)
{
	char *retstr;
	char *lastdot;
	if (mystr == NULL)
		return NULL;
	if ((retstr = (char*)malloc (strlen (mystr) + 1)) == NULL)
		return NULL;
	strcpy (retstr, mystr);
	lastdot = strrchr (retstr, '.');
	if (lastdot != NULL)
		*lastdot = '\0';
	return retstr;
}

char* get_first_match(const char *str_regex, const char *str)
{
	regex_t preg;
	regmatch_t *pmatch = NULL;

	//Compilation de l'expression réguliére :
	if( regcomp(&preg, str_regex, REG_EXTENDED) != 0 )
	{
		fprintf(stderr, "%s::%s::%d :: Bug dans la recherche de l'indice !\n", __FILE__, __func__, __LINE__);
		regfree(&preg);
		return NULL;
	}

	//Allocation du tableau contenant les résultats de l'expression :
	if( (pmatch = malloc(sizeof(*pmatch) * preg.re_nsub)) == NULL )
	{
		fprintf(stderr, "%s::%s::%d :: Bug dans la recherche de l'indice !\n", __FILE__, __func__, __LINE__);
		regfree(&preg);
		return NULL;
	}

	//On applique l'expression :
	int res = regexec (&preg, str, preg.re_nsub, pmatch, 0);

	//Si on ne trouve aucun motif correspondant :
	if( res == REG_NOMATCH )
	{
		//fprintf(stderr, "%s don't match %s.\n", str_regex, filename);
		regfree(&preg);
		free(pmatch);
		return NULL;
	}
	//Ou si l'erreur est plus profonde :
	else if( res != 0 )
	{
		char *text; size_t size;
		size = regerror(res, &preg, NULL, 0);
		if( (text = malloc (sizeof (*text) * size)) == NULL )
		{
			perror("Memoire insuffisante\n");
			regfree(&preg);
			free(pmatch);
			exit(EXIT_FAILURE);
		}
		regerror(res, &preg, text, size);
		fprintf(stderr, "%s\n", text);
		free(text);
		return NULL;
	}

	char *site = NULL;
	int start = pmatch[0].rm_so;
	int end = pmatch[0].rm_eo;
	size_t size = end - start;

	if( (site = malloc (sizeof (*site) * (size + 1))) == NULL )
	{
		perror("Mémoire insuffisante\n");
		regfree(&preg);
		exit(EXIT_FAILURE);
	}

	strncpy (site, &str[start], size);
	site[size] = '\0';

	free(pmatch);
	regfree(&preg);

	return site;
}

int get_id(const char *filename)
{
	int start = 0;
	char str_regex[] = "([[:digit:]]+)$", *site = NULL;

	site = get_first_match(str_regex, filename);
	if( site == NULL )
		return -1;

	while(site[start] != '0')
		start++;
	int ind = atoi(&site[start]);

	free (site);

	return ind;
}

void lissage(double *tab, const int N)
{
	if( tab[0] == 0.0 )
	{
		tab[0] = ( tab[0] + tab[1] )/2.0;
#ifdef PREV_LISSAGE_P
		fprintf(stderr, "%s::[%d]::Lissage\n", __func__, 0);
#endif
	}

	for (int i = 2; i < N-1; i++)
	{
		if( tab[i] == 0.0 )
		{
			tab[i] = (tab[i-1] + tab[i] + tab[i+1]) / 3.0;
#ifdef PREV_LISSAGE_P
			fprintf(stderr, "%s::[%d]::Lissage\n", __func__, i);
#endif
		}
	}

	if( tab[N-1] == 0.0 )
	{
		tab[N-1] = ( tab[N-1] + tab[N-2] ) / 2.0;
#ifdef PREV_LISSAGE_P
		fprintf(stderr, "%s::[%d]::Lissage\n", __func__, N-1);
#endif
	}
}

unsigned int* unsignedint1d(const unsigned int n)
{
	unsigned int * ptf=NULL;
	if (n >0) {
		ptf = (unsigned int *) calloc(n, sizeof(unsigned int)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation unsigned int1d\n");
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void unsignedint1d_libere(unsigned int *ptf)
{
	/* liberation d'un tableau 1D de unsigned ints */
	free(ptf);
	return ;
}

int* int1d(const int n)
{
	int * ptf=NULL;
	if (n >0) {
		ptf = (int *) calloc(n, sizeof(int)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation int1d\n");
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void int1d_libere(int *ptf)
{
	/* liberation d'un tableau 1D de ints */
	free(ptf);
	return ;
}

float* float1d(const int n)
{
	float * ptf=NULL;
	if (n >0) {
		ptf = (float *) calloc(n, sizeof(float)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation float1d\n");
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void float1d_libere(float *ptf)
{
	/* liberation d'un tableau 1D de floats */
	free(ptf);
	return ;
}

double* double1d(const int n)
{
	double * ptf=NULL;
	if (n >0) {
		ptf = (double *) calloc(n, sizeof(double)) ;
	}
	if (ptf==NULL) {
		perror("Erreur allocation double1d ");
		//fprintf(stderr, "Erreur allocation double1d : %p\n", ptf);
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void double1d_libere(double *ptf)
{
	/* liberation d'un tableau 1D de doubles */
	free(ptf);
	return ;
}

int maxlocdouble1d(const double * x, const int n)
{
	double tmp;
	int imax;
	int i;
	tmp = x[0];
	imax = 0 ;
	for (i=0; i<n; i++){
		if( x[i] > tmp) {
			tmp = x[i];
			imax = i ;
		}
	}
	return imax;
}

double maxdouble2d(const double **tab, const int NbPart, const int col)
{
	double max = tab[0][col];
	for(int i=1; i<NbPart; i++)
		if( tab[i][col] > max )
			max = tab[i][col];
	return max;
}

int minlocdouble1d(const double * x, const int n)
{
	double tmp;
	int imax;
	int i;
	tmp = x[0];
	imax = 0 ;
	for (i=1; i<n; i++){
		if( x[i] < tmp) {
			tmp = x[i];
			imax = i ;
		}
	}
	return imax;
}

void double2d_libere(double ** mat)
{
	free(mat[0]); /* liberation du pointeur de la matrice */
	free(mat) ; /* liberation du vecteur des pointeurs de debut de ligne */
}

double ** double2d(const int nlignes, const int mcolonnes)
{
	/* allocation d'un tableau 2D de nlignes x mcolonnes doubles */
	double * ptf=NULL, ** pt=NULL;
	if (nlignes >0 && mcolonnes >0) { /* Allocation globale de la matrice */
		ptf = (double *) calloc(nlignes*mcolonnes, sizeof(double)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation double2d \n");
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	/* allocation des n pointeurs de debut de ligne */
	pt = (double **) calloc(nlignes, sizeof(double*));
	if (pt==NULL) {
		fprintf(stderr, "erreur allocation double2d \n");
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	/* affectation des pointeurs de debut de ligne */
	for(int k=0; k < nlignes; k++) {
		pt[k] = &(ptf[k*mcolonnes]);
	}
	return pt;
}

