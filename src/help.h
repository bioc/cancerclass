

/* functions */
double testx(double tmp);
int vergleich(int *i1, int *i2);
 
int compare_doubles (const void *x, const void *y);

double median(double *values, int from, int to);
double madian(double *values, int from, int to);
double sum(double *values, int from, int to);
double quantile(double *values, int from, int to, double r);
double quantile7(double *values, int from, int to, double r);
double mean(double *values, int from, int to);
double var(double *array, int from, int to);
double IQR(double *values, int from, int to);
double student_t(double *values1, double *values2, int ng1, int ng2);
double welch_t(double *values1, double *values2, int ng1, int ng2);
double welch_df(double *values1, double *values2, int ng1, int ng2);
double p_value(double stat, double df);
double p_wilcoxon(double stat, double m, double n);
double getmax(double *values, int form, int to);

void showarray(double *array, int from, int to);
void showvalue(double x);

void sort(double *values, int nsamples);
void sorttmp(int);
void bsortdesc(double *array, int from, int to);
void bsort(double *array, int from, int to);

void writedata(double *array, int size, char *filename);
void writeheader(char *filename);

int mini(double a, double b);
//void R_qsort(double *v, int i, int j);
