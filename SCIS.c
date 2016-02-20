#include <libgimp/gimp.h>
#include <stdbool.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>


/**
 * To know program state and printing error message
 */
static int STATE;
static bool interactive_run_mode;
static char* mallocError;
void errror_handler(){

    char* msg;
    switch(STATE){
    case 1:
        /*No seed found*/
        msg="SCIS doesn't find any seeds. Please check your seeds layer.";
        break;
    case 2:
        /*Too much seeds*/
        msg="SCIS find more than 100 different colors for seeds : there is probably something wrong with your seed layers.";

        break;
    case 3:
        msg="Give at least 2 layers  for image and seeds.";
        break;
    case 4:
        msg="Segmentation fail.";
        break;
    case 5:
        msg="Something wrong during SVM training.";
        break;
    default:
        msg="Unknow error";
    }


    if(interactive_run_mode){
        g_message("%s\n",msg);
    }else{
        printf("%s\n",msg);
    }

    if(STATE==6){

        if(interactive_run_mode){
            g_message("%s\n",mallocError);
        }else{
            printf("%s\n",mallocError);
        }
    }

}


/**
 * @class SuperpixelFeatures
 * @brief Superpixel descriptor
 *
 * A superpixel is describe by :
 *  - its average color
 *  - its center of mass coordinates
 *  - its number of pixels
 *  - its class
 */
typedef struct {
    double colorAtt1;/**< first average color channel*/
    double colorAtt2;/**< second average color channel*/
    double colorAtt3;/**< third average color channel*/
    double r;/**< center of mass ordinate */
    double c;/**< center of mass abscissa */
    int nbPixels;/**< number of pixels belonging to the superpixel */
    int classLabel;/**< class label of the superpixel*/
}SuperpixelFeatures;

/************************************************************************************************************************
 *
 * Superpixel segmentation : first part
 *
 * Adapted from C++ implementation of Pedro Felzenszwalb
 * Copyright (c) Pedro Felzenszwalb
 *
 * http://cs.brown.edu/~pff/segment/
 *
 ************************************************************************************************************************/


/* threshold function for initialise graph edges*/
#define THRESHOLD(size, c) (c/size)

typedef struct {
    int rank;/**< set rank */
    int p;/**< set parent*/
    int size;/**< cardinality of subset */
} DisjointSet_elt;

typedef struct{
    DisjointSet_elt *elts; /**< subsets */
    int num;/**< number of subsets */
} DisjointSet;

/**
 * @brief create_disjoint_set create data structure disjoint set
 * @param nbElements number of set elements
 * @return
 */
DisjointSet* create_disjoint_set(int nbElements) {
    DisjointSet * res=(DisjointSet *)malloc(sizeof(DisjointSet));

    res->elts = (DisjointSet_elt *) malloc(sizeof(DisjointSet_elt)*nbElements);
    res->num = nbElements;
    {
        int i;
        for (i = 0; i < nbElements; i++) {
            res->elts[i].rank = 0;
            res->elts[i].size = 1;
            res->elts[i].p = i;
        }
    }
    return res;
}

/**
 * @brief delete_disjoint_set clean up memory
 * @param dSet
 */
void delete_disjoint_set(DisjointSet** dSet){
    free((*dSet)->elts);
    free(*dSet);
}

/**
 * @brief disjointSet_find determine which subset a particular element is in
 * @param dSet set
 * @param x element
 * @return
 */
int disjointSet_find(DisjointSet* dSet,int x) {

    int y = x;
    /*convert recursion of algorithm to iteration*/
    while (y != dSet->elts[y].p)
        y = dSet->elts[y].p;
    /*update parent of the subset of index x */
    dSet->elts[x].p = y;
    /*return parent*/
    return y;
}

/**
 * @brief disjointSet_join join two subsets into a single subset.
 * @param dSet set
 * @param subSet1 subset
 * @param subSet2 subset
 */
void disjointSet_join(DisjointSet * dSet,int subSet1, int subSet2) {
    /* for efficiency reasons :
       always attach the smaller tree
       to the root of the larger tree */
    if (dSet->elts[subSet1].rank > dSet->elts[subSet2].rank) {
        /* add subset y to the subset x */
        dSet->elts[subSet2].p = subSet1;
        dSet->elts[subSet1].size += dSet->elts[subSet2].size;
    } else {
        /* add subset x to the subset y */
        dSet->elts[subSet1].p = subSet2;
        dSet->elts[subSet2].size += dSet->elts[subSet1].size;
        /* choose y subset to have the higher rank */
        if (dSet->elts[subSet1].rank == dSet->elts[subSet2].rank)
            dSet->elts[subSet2].rank++;
    }
    dSet->num--;
}

/**
 * @brief get_subset_size : return subset cardinality
 * @param set data set
 * @param subSet
 * @return
 */
int get_subset_size(DisjointSet* dSet,int subSet)  {
    return dSet->elts[subSet].size;
}

/**
 * @brief get_nb_subsets return number of subset
 * @param dSet
 * @return
 */
int get_nb_subsets(DisjointSet* dSet)  {
    return dSet->num;
}


/**
 * Disjoint graph edge between two subsets
 */
typedef struct {
    float w;
    int subSet1, subSet2;
} GraphEdge;


/* compare to edge) */
static int compare_edges (void const *subSet1, void const *subSet2)
{
    /* definir des pointeurs type's et initialise's
      avec les parametres */
    GraphEdge const *p1 = subSet1;
    GraphEdge const *p2 = subSet2;

    /* evaluer et retourner l'etat de l'evaluation (tri croissant) */
    return p1->w - p2->w;
}

/**
 * @brief felzenszwalb_oversegmentation return a disjoint-set forest representing the segmentation
 * @param num_vertices[i] number of vertices in graph
 * @param num_edges[i] number of edges in graph
 * @param edges[i] array of edges
 * @param c[i] constant for treshold function
 * @return a disjoint-set forest representing the segmentation
 */
DisjointSet * felzenszwalb_oversegmentation(int num_vertices, int num_edges, GraphEdge *edges,
                                            float c) {
    int i;
    int j;
    DisjointSet* dSet;
    float * threshold;
    /* sort edges by weight */
    qsort(edges,num_edges, sizeof(GraphEdge), compare_edges);

    /* make a disjoint-set forest */
    dSet = create_disjoint_set(num_vertices);

    /* init thresholds */
    threshold = (float*) malloc(sizeof(float)*num_vertices);

    for ( i = 0; i < num_vertices; i++){
        threshold[i] = THRESHOLD(1,c);
    }
    /* for each edge, in ascending weight order : */
    for ( i = 0; i < num_edges; i++) {
        GraphEdge *pedge = &edges[i];

        /* components conected by this edge */
        int subSet1 = disjointSet_find(dSet,pedge->subSet1);
        int subSet2 = disjointSet_find(dSet,pedge->subSet2);
        if (subSet1 != subSet2) {
            if ((pedge->w <= threshold[subSet1]) &&
                    (pedge->w <= threshold[subSet2])) {
                disjointSet_join(dSet,subSet1, subSet2);
                subSet1 = disjointSet_find(dSet,subSet1);
                threshold[subSet1] = pedge->w + THRESHOLD(get_subset_size(dSet,subSet1), c);
            }
        }
    }

    /* clean up memory */
    free(threshold);
    return dSet;
}



/**
 *
 * SVM of type C_SVC with RBF kernel implementation
 *
 * Adapted of LibSVM c++ implementation
 *
 * Copyright (c) Chih-Chung Chang and Chih-Jen Lin
 * https://www.csie.ntu.edu.tw/~cjlin/libsvm/
 *
 **/

typedef struct
{
    int index;
    double value;
}svm_node;

typedef struct
{
    int l;
    double *y;
    svm_node **x;
}svm_problem;


typedef struct
{
    int degree;	/**< for poly */
    double gamma;	/**< for poly/rbf/sigmoid */

    /* these are for training only */
    double cache_size; /**< in MB */
    double eps;	/**< stopping criteria */
    double C;	/**< for C_SVC, EPSILON_SVR and NU_SVR */
    int nr_weight;		/**< for C_SVC */
    int *weight_label;	/**< for C_SVC */
    double* weight;		/**< for C_SVC */
    int shrinking;	/**< use the shrinking heuristics */
} svm_parameter;


typedef struct
{
    svm_parameter param;	/**< parameter */
    int nr_class;		/**< number of classes, = 2 in regression/one class svm */
    int l;			/**< total #SV */
    svm_node **SV;		/**< SVs (SV[l]) */
    double **sv_coef;	/**< coefficients for SVs in decision functions (sv_coef[k-1][l]) */
    double *rho;		/**< constants in decision functions (rho[k*(k-1)/2]) */
    double *probA;		/**< pariwise probability information */
    double *probB;
    int *sv_indices;        /**< sv_indices[0,...,nSV-1] are values in [1,...,num_traning_data] to indicate SVs in the training set */

    /* for classification only */

    int *label;		/**< label of each class (label[k]) */
    int *nSV;		/**< number of SVs for each class (nSV[k])
                        nSV[0] + nSV[1] + ... + nSV[k-1] = l */
    /* XXX */
    int free_sv;		/**< 1 if svm_model is created by svm_load_model
                             0 if svm_model is created by svm_train */
} svm_model;


svm_model * svm_train(svm_problem *prob, svm_parameter *param);


int svm_get_svm_type( svm_model *model);
int svm_get_nr_class( svm_model *model);
void svm_get_labels( svm_model *model, int *label);
void svm_get_sv_indices( svm_model *model, int *sv_indices);


void svm_free_model_content( svm_model *model_ptr);
void svm_free_and_destroy_model( svm_model **model_ptr_ptr);
void svm_destroy_param( svm_parameter *param);


void svm_set_print_string_function(void (*print_func)(const char *));

typedef float Qfloat;
typedef signed char schar;


double svm_predict_values( svm_model *model,  svm_node *x, double* dec_values);
double svm_predict( svm_model *model,  svm_node *x);
double svm_predict_probability( svm_model *model,  svm_node *x, double* prob_estimates);

double min(double x,double y){
    return (x<y)?x:y;
}


double max(double x,double y){
    return (x>y)?x:y;
}

void swap_double(double* x,double* y){
    double t=*x;
    *x=*y;
    *y=t;
}

void swap_int(int* x,int* y){
    int t=*x;
    *x=*y;
    *y=t;
}

void swap_schar(schar* x,schar* y){
    int t=*x;
    *x=*y;
    *y=t;
}

void swap_char(char* x,char* y){
    int t=*x;
    *x=*y;
    *y=t;
}

void swap_svm_node( svm_node **x, svm_node **y){
    svm_node *t=*x;
    *x=*y;
    *y=t;
}
typedef enum { LOWER_BOUND, UPPER_BOUND, FREE }SOLVER_STATUS;

void swap_SOLVER_STATUS( SOLVER_STATUS *x, SOLVER_STATUS*y){
    SOLVER_STATUS t=*x;
    *x=*y;
    *y=t;
}


double powi(double base, int times)
{
    double tmp = base, ret = 1.0;

    {
        int t;
        for(t=times; t>0; t/=2)
        {
            if(t%2==1) ret*=tmp;
            tmp = tmp * tmp;
        }
    }
    return ret;
}
#define INF HUGE_VAL
#define TAU 1e-12

void print_string_stdout(const char *s)
{
    fputs(s,stdout);
    fflush(stdout);
}

void (*svm_print_string) (const char *) = &print_string_stdout;


/* kernel Cache
 * l is the number of total data items
 * size is the cache size limit in bytes
 */


typedef struct head_t
{
    struct head_t *prev, *next;	/* a circular list */
    Qfloat *data;
    int len;
}head_t;

typedef struct {
    int l;
    long int size;
    head_t *head;
    head_t lru_head;

}Cache;

Cache* create_cache(int l,long int size){
    Cache* res=( Cache*)malloc(sizeof( Cache));
    res->head = (head_t *)calloc(l,sizeof(head_t));
    res->size /= sizeof(Qfloat);
    res->size -= l * sizeof(head_t) / sizeof(Qfloat);
    res->size = max(size, 2 * (long int) l);	/* cache must be large enough for two columns */
    res->lru_head.next = res->lru_head.prev = &(res->lru_head);
    return res;
}
void delete_cache( Cache** cache){
    head_t *h ;
    for(h = (*cache)->lru_head.next; h != &((*cache)->lru_head); h=h->next)
        free(h->data);
    free((*cache)->head);
    free(*cache);
}

/* delete from current location */
void cache_lru_delete(head_t *h)
{
    h->prev->next = h->next;
    h->next->prev = h->prev;
}

/* insert to last position */
void cache_lru_insert( Cache* cache,head_t *h)
{
    h->next = &(cache->lru_head);
    h->prev = cache->lru_head.prev;
    h->prev->next = h;
    h->next->prev = h;
}


/* request data [0,len)
   return some position p where [p,len) need to be filled
   (p >= len if nothing needs to be filled) */
int cache_get_data(Cache* cache,const int index, Qfloat **data, int len){
    head_t *h = &(cache->head[index]);
    int more;
    if(h->len) cache_lru_delete(h);
    more = len - h->len;

    if(more > 0)
    {
        /* free old space */
        while(cache->size < more)
        {
            head_t *old = cache->lru_head.next;
            cache_lru_delete(old);
            free(old->data);
            cache->size += old->len;
            old->data = 0;
            old->len = 0;
        }

        /* allocate new space */
        h->data = (Qfloat *)realloc(h->data,sizeof(Qfloat)*len);
        cache->size -= more;
        swap_int(&(h->len),&(len));
    }

    cache_lru_insert(cache,h);
    *data = h->data;
    return len;
}


void cache_swap_index( Cache* cache,int i, int j)
{
    Qfloat* saveData;
    head_t *h;
    if(i==j) return;

    if(cache->head[i].len) cache_lru_delete(&(cache->head[i]));
    if(cache->head[j].len) cache_lru_delete(&(cache->head[j]));
    saveData=cache->head[i].data;
    cache->head[i].data=cache->head[j].data;
    cache->head[j].data=saveData;
    swap_int(&(cache->head[i].len),&(cache->head[j].len));
    if(cache->head[i].len) cache_lru_insert(cache,&(cache->head[i]));
    if(cache->head[j].len) cache_lru_insert(cache, &(cache->head[j]));

    if(i>j) swap_int(&i,&j);
    for(h = cache->lru_head.next; h!=&(cache->lru_head); h=h->next)
    {
        if(h->len > i)
        {
            if(h->len > j){

                saveData=cache->head[i].data;
                cache->head[i].data=cache->head[j].data;
                cache->head[j].data=saveData;
            }
            else
            {
                /* give up */
                cache_lru_delete(h);
                free(h->data);
                cache->size += h->len;
                h->data = 0;
                h->len = 0;
            }
        }
    }
}



typedef struct {
    svm_node **x;
    double *x_square;

    /* svm_parameter */
    int degree;
    double gamma;
    double coef0;
    int l;
    schar *y;
    Cache *cache;
    double *QD;
}Kernel;


double kernel_k_function( svm_node *x,  svm_node *y,
                          svm_parameter* param){

    double sum = 0;

    while(x->index != -1 && y->index !=-1)
    {


        if(x->index == y->index)
        {
            double d = x->value - y->value;
            sum += d*d;
            ++x;
            ++y;

        }
        else
        {
            if(x->index > y->index)
            {
                sum += y->value * y->value;
                ++y;
            }
            else
            {
                sum += x->value * x->value;
                ++x;
            }
        }
    }

    while(x->index != -1)
    {
        sum += x->value * x->value;
        ++x;
    }

    while(y->index != -1)
    {
        sum += y->value * y->value;
        ++y;
    }

    return exp(-param->gamma*sum);
}

double kernel_dot(const svm_node *px, const svm_node *py){
    double sum = 0;
    while(px->index != -1 && py->index != -1)
    {
        if(px->index == py->index)
        {
            sum += px->value * py->value;
            ++px;
            ++py;
        }
        else
        {
            if(px->index > py->index)
                ++py;
            else
                ++px;
        }
    }
    return sum;
}

/**
 * @brief kernel_k_function only compute RBF kernel
 * @param i
 * @param j
 * @return
 */
double kernel_k_function2( Kernel * kernel,int i,int j){
    return exp(-kernel->gamma*(kernel->x_square[i]+kernel->x_square[j]-2*kernel_dot(kernel->x[i],kernel->x[j])));

}

Qfloat * kernel_get_Q( Kernel* kernel,int i, int len)
{
    Qfloat *data;
    int start, j;

    if((start = cache_get_data(kernel->cache,i,&data,len)) < len)
    {
        for(j=start;j<len;j++)
            data[j] = (Qfloat)(kernel->y[i]*kernel->y[j]*kernel_k_function2(kernel,i,j));
    }
    return data;
}

double * kernel_get_QD( Kernel* kernel)
{
    return kernel->QD;
}


void kernel_swap_index( Kernel* kernel,int i,int j){
    cache_swap_index(kernel->cache,i,j);

    swap_svm_node(&(kernel->x[i]),&(kernel->x[j]));
    if(kernel->x_square) swap_double(&(kernel->x_square[i]),&(kernel->x_square[j]));


    swap_schar(&(kernel->y[i]),&(kernel->y[j]));
    swap_double(&(kernel->QD[i]),&(kernel->QD[j]));

}

void clone(svm_node** dst, svm_node* src, int n)
{
    *dst = (svm_node*) malloc(sizeof(svm_node)*n);
    memcpy((void *)*dst,(void *)src,sizeof(svm_node)*n);
}


void clone_schar(schar** dst, schar* src, int n)
{
    *dst = (schar*) malloc(sizeof(schar)*n);
    memcpy((void *)dst,(void *)src,sizeof(schar)*n);
}



Kernel * create_kernel(  svm_problem* prob,  svm_parameter* param, const schar *y_){
    Kernel* res=( Kernel*)malloc(sizeof( Kernel));

    int i;

    res->degree=param->degree;
    res->gamma=param->gamma;
    res->x=(svm_node**)malloc(sizeof(svm_node*)*prob->l);
    for( i=0;i<prob->l;i++){
        res->x[i]=prob->x[i];
    }


    res->x_square = (double*) malloc(sizeof(double)*prob->l);
    for(i=0;i<prob->l;i++){
        res->x_square[i] = kernel_dot(prob->x[i],prob->x[i]);
    }

    res->l=prob->l;

    res->y=(schar*)malloc(sizeof(schar)*prob->l);
    {
        int i;
        for( i=0;i<prob->l;i++){
            res->y[i]=y_[i];
        }
    }

    res->cache = create_cache(prob->l,(long int)(param->cache_size*(1<<20)));
    res->QD = (double*) malloc(sizeof(double)*prob->l);
    {
        int i;
        for(i=0;i<prob->l;i++){
            res->QD[i] = kernel_k_function2(res,i,i);
        }
    }
    return res;
}

void delete_kernel( Kernel** kernel){


    free((*kernel)->x_square);
    free((*kernel)->x);
    free((*kernel)->y);
    free((*kernel)->QD);
    delete_cache(&((*kernel)->cache));
    free(*kernel);

}


/*An SMO algorithm in Fan et al., JMLR 6(2005), p. 1889--1918
 Solves:

   min 0.5(\alpha^T Q \alpha) + p^T \alpha

       y^T \alpha = \delta
       y_i = +1 or -1
       0 <= alpha_i <= Cp for y_i = 1
       0 <= alpha_i <= Cn for y_i = -1

Given:

   Q, p, y, Cp, Cn, and an initial feasible point \alpha
   l is the size of vectors and matrices
   eps is the stopping tolerance

solution will be put in \alpha, objective value will be put in obj*/

typedef struct  {
    double obj;
    double rho;
    double upper_bound_p;
    double upper_bound_n;
} SolutionInfo;



typedef struct{
    int active_size;
    schar *y;
    double *G;		/* gradient of objective function */
    SOLVER_STATUS *alpha_status;	/* LOWER_BOUND, UPPER_BOUND, FREE */
    double *alpha;
    Kernel *Q;
    const double *QD;
    double eps;
    double Cp,Cn;
    double *p;
    int *active_set;
    double *G_bar;		/* gradient, if we treat free variables as 0 */
    int l;
    bool unshrink;

} Solver;

Solver * create_solver(){
    Solver* res=( Solver*)malloc(sizeof( Solver));

    return res;
}



void delete_solver( Solver** solver){
    free(*solver);
}


double solver_get_C( Solver* solver,int i)
{
    return (solver->y[i] > 0)? solver->Cp : solver->Cn;
}
void solver_update_alpha_status( Solver* solver,int i)
{
    if(solver->alpha[i] >= solver_get_C(solver,i))
        solver->alpha_status[i] = UPPER_BOUND;
    else if(solver->alpha[i] <= 0)
        solver->alpha_status[i] = LOWER_BOUND;
    else solver->alpha_status[i] = FREE;
}
bool solver_is_upper_bound( Solver* solver,int i) {
    return solver->alpha_status[i] == UPPER_BOUND;
}
bool solver_is_lower_bound( Solver* solver,int i) {
    return solver->alpha_status[i] == LOWER_BOUND;
}
bool solver_is_free( Solver* solver, int i) {
    return solver->alpha_status[i] == FREE;
}
void solver_swap_index( Solver* solver,int i, int j)
{
    kernel_swap_index(solver->Q,i,j);
    swap_schar(&(solver->y[i]),&(solver->y[j]));
    swap_double(&(solver->G[i]),&(solver->G[j]));
    swap_SOLVER_STATUS(&(solver->alpha_status[i]),&(solver->alpha_status[j]));
    swap_double(&(solver->alpha[i]),&(solver->alpha[j]));
    swap_double(&(solver->p[i]),&(solver->p[j]));
    swap_int(&(solver->active_set[i]),&(solver->active_set[j]));
    swap_double(&(solver->G_bar[i]),&(solver->G_bar[j]));
}

void solver_reconstruct_gradient( Solver* solver){
    int i,j,nrFree;
    /* reconstruct inactive elements of G from G_bar and free variables */

    if(solver->active_size ==solver->l) return;

    nrFree = 0;

    for(j=solver->active_size;j<solver->l;j++)
        solver->G[j] = solver->G_bar[j] + solver->p[j];

    for(j=0;j<solver->active_size;j++)
        if(solver_is_free(solver,j))
            nrFree++;

    if(2*nrFree < solver->active_size)
        printf("\nWARNING: using -h 0 may be faster\n");

    if (nrFree*solver->l > 2*solver->active_size*(solver->l-solver->active_size))
    {
        for(i=solver->active_size;i<solver->l;i++)
        {
            const Qfloat *Q_i = kernel_get_Q(solver->Q,i,solver->active_size);
            for(j=0;j<solver->active_size;j++)
                if(solver_is_free(solver,j))
                    solver->G[i] += solver->alpha[j] * Q_i[j];
        }
    }
    else
    {
        for(i=0;i<solver->active_size;i++)
            if(solver_is_free(solver,i))
            {
                const Qfloat *Q_i = kernel_get_Q(solver->Q,i,solver->l);
                double alpha_i = solver->alpha[i];
                for(j=solver->active_size;j<solver->l;j++)
                    solver->G[j] += alpha_i * Q_i[j];
            }
    }
}



/**
 * @brief solver_select_working_set return 1 if already optimal, return 0 otherwise
 *
 * return i,j such that
 * i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
 * j: minimizes the decrease of obj value
 * (if quadratic coefficeint <= 0, replace it with tau)
 * -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
 *
 * @param solver
 * @param out_i
 * @param out_j
 * @return
 */
int solver_select_working_set( Solver* solver,int * out_i, int *out_j)
{

    double Gmax = -INF;
    double Gmax2 = -INF;
    int Gmax_idx = -1;
    int Gmin_idx = -1;
    double obj_diff_min = INF;
    int t;

    int i;
    const Qfloat *Q_i = NULL;
    int j;
    for( t=0;t<solver->active_size;t++)
        if(solver->y[t]==+1)
        {
            if(!solver_is_upper_bound(solver,t))
                if(-solver->G[t] >= Gmax)
                {
                    Gmax = -solver->G[t];
                    Gmax_idx = t;
                }
        }
        else
        {
            if(!solver_is_lower_bound(solver,t))
                if(solver->G[t] >= Gmax)
                {
                    Gmax = solver->G[t];
                    Gmax_idx = t;
                }
        }

    i = Gmax_idx;
    if(i != -1) /* NULL Q_i not accessed: Gmax=-INF if i=-1 */
        Q_i = kernel_get_Q(solver->Q,i,solver->active_size);

    for( j=0;j<solver->active_size;j++)
    {
        if(solver->y[j]==+1)
        {
            if (!solver_is_lower_bound(solver,j))
            {
                double grad_diff=Gmax+solver->G[j];
                if (solver->G[j] >= Gmax2)
                    Gmax2 = solver->G[j];
                if (grad_diff > 0)
                {
                    double obj_diff;
                    double quad_coef = solver->QD[i]+solver->QD[j]-2.0*solver->y[i]*Q_i[j];
                    if (quad_coef > 0)
                        obj_diff = -(grad_diff*grad_diff)/quad_coef;
                    else
                        obj_diff = -(grad_diff*grad_diff)/TAU;

                    if (obj_diff <= obj_diff_min)
                    {
                        Gmin_idx=j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
        else
        {
            if (!solver_is_upper_bound(solver,j))
            {
                double grad_diff= Gmax-solver->G[j];
                if (-solver->G[j] >= Gmax2)
                    Gmax2 = -solver->G[j];
                if (grad_diff > 0)
                {
                    double obj_diff;
                    double quad_coef = solver->QD[i]+solver->QD[j]+2.0*solver->y[i]*Q_i[j];
                    if (quad_coef > 0)
                        obj_diff = -(grad_diff*grad_diff)/quad_coef;
                    else
                        obj_diff = -(grad_diff*grad_diff)/TAU;

                    if (obj_diff <= obj_diff_min)
                    {
                        Gmin_idx=j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
    }

    if(Gmax+Gmax2 < solver->eps){
        return 1;
    }

    *out_i = Gmax_idx;
    *out_j = Gmin_idx;
    return 0;
}

bool solver_be_shrunk( Solver* solver,int i, double Gmax1, double Gmax2)
{
    if(solver_is_upper_bound(solver,i))
    {
        if(solver->y[i]==+1)
            return(-solver->G[i] > Gmax1);
        else
            return(-solver->G[i] > Gmax2);
    }
    else if(solver_is_lower_bound(solver,i))
    {
        if(solver->y[i]==+1)
            return(solver->G[i] > Gmax2);
        else
            return(solver->G[i] > Gmax1);
    }
    else
        return(false);
}

void solver_do_shrinking( Solver* solver)
{
    int i;
    double Gmax1 = -INF;
    double Gmax2 = -INF;

    /* find maximal violating pair first */
    for(i=0;i<solver->active_size;i++)
    {
        if(solver->y[i]==+1)
        {
            if(!solver_is_upper_bound(solver,i))
            {
                if(-solver->G[i] >= Gmax1)
                    Gmax1 = -solver->G[i];
            }
            if(!solver_is_lower_bound(solver,i))
            {
                if(solver->G[i] >= Gmax2)
                    Gmax2 = solver->G[i];
            }
        }
        else
        {
            if(!solver_is_upper_bound(solver,i))
            {
                if(-solver->G[i] >= Gmax2)
                    Gmax2 = -solver->G[i];
            }
            if(!solver_is_lower_bound(solver,i))
            {
                if(solver->G[i] >= Gmax1)
                    Gmax1 = solver->G[i];
            }
        }
    }

    if(solver->unshrink == false && Gmax1 + Gmax2 <= solver->eps*10)
    {
        solver->unshrink = true;
        solver_reconstruct_gradient(solver);
        solver->active_size = solver->l;
    }

    for(i=0;i<solver->active_size;i++)
        if (solver_be_shrunk(solver,i, Gmax1, Gmax2))
        {
            solver->active_size--;
            while (solver->active_size > i)
            {
                if (!solver_be_shrunk(solver,solver->active_size, Gmax1, Gmax2))
                {
                    solver_swap_index(solver,i,solver->active_size);
                    break;
                }
                solver->active_size--;
            }
        }
}

double solver_calculate_rho( Solver* solver)
{
    double r;
    int nr_free = 0;
    double ub = INF, lb = -INF, sum_free = 0;
    int i;
    for(i=0;i<solver->active_size;i++){
        double yG = solver->y[i]*solver->G[i];

        if(solver_is_upper_bound(solver,i))
        {
            if(solver->y[i]==-1)
                ub = min(ub,yG);
            else
                lb = max(lb,yG);
        }
        else if(solver_is_lower_bound(solver,i))
        {
            if(solver->y[i]==+1)
                ub = min(ub,yG);
            else
                lb = max(lb,yG);
        }
        else
        {
            ++nr_free;
            sum_free += yG;
        }
    }

    if(nr_free>0)
        r = sum_free/nr_free;
    else
        r = (ub+lb)/2;

    return r;
}

void solver_solve( Solver * solver,int l,  svm_problem* prob,  svm_parameter* param, const double *p_, const schar *y_,
                   double *alpha_, double Cp, double Cn, double eps,
                   SolutionInfo* si, int shrinking)
{

    int iter ;
    int max_iter ;
    int counter ;

    solver->l = l;
    solver->Q = create_kernel(prob,param,y_);
    solver->QD=kernel_get_QD(solver->Q);
    solver->p=(double*)malloc(sizeof(double)*l);
    solver->y=(schar*)malloc(sizeof(schar)*l);
    solver->alpha=(double*)malloc(sizeof(double)*l);
    {
        int i;
        for(i=0;i<l;i++){
            solver->p[i]=p_[i];
            solver->y[i]=y_[i];
            solver->alpha[i]=alpha_[i];
        }
    }
    solver->Cp = Cp;
    solver->Cn = Cn;
    solver->eps = eps;
    solver->unshrink = false;

    /* initialize alpha_status */
    {
        int i;
        solver->alpha_status = (SOLVER_STATUS*) malloc(sizeof(SOLVER_STATUS)*l);
        for(i=0;i<l;i++)
            solver_update_alpha_status(solver,i);
    }

    /* initialize active set (for shrinking) */
    {
        int i;
        solver->active_set = (int *) malloc(sizeof(int)*l);
        for( i=0;i<l;i++) solver->active_set[i] = i;
        solver->active_size = l;
    }

    /* initialize gradient */
    {
        int i;
        solver->G = (double*) malloc(sizeof(double)*l);
        solver->G_bar = (double*) malloc(sizeof(double)*l);

        for(i=0;i<l;i++)
        {
            solver->G[i] = solver->p[i];
            solver->G_bar[i] = 0;
        }
        for(i=0;i<l;i++)
            if(!solver_is_lower_bound(solver,i))
            {
                const Qfloat *Q_i = kernel_get_Q(solver->Q,i,l);
                double alpha_i = solver->alpha[i];
                int j;
                for(j=0;j<l;j++)
                    solver->G[j] += alpha_i*Q_i[j];
                if(solver_is_upper_bound(solver,i))
                    for(j=0;j<l;j++)
                        solver->G_bar[j] += solver_get_C(solver,i) * Q_i[j];
            }
    }

    /* optimization step */

    iter = 0;
    max_iter = max(10000000, l>INT_MAX/100 ? INT_MAX : 100*l);
    counter = min(l,1000)+1;

    while(iter < max_iter)
    {
        int i,j;
        const Qfloat* Q_i, *Q_j;
        double C_i, C_j, old_alpha_i, old_alpha_j,delta_alpha_i,delta_alpha_j;
        /* show progress and do shrinking */

        if(--counter == 0)
        {
            counter = min(l,1000);
            if(shrinking) solver_do_shrinking(solver);
        }

        if(solver_select_working_set(solver,&i,&j)!=0)
        {
            /* reconstruct the whole gradient */
            solver_reconstruct_gradient(solver);
            /* reset active set size and check */
            solver->active_size = l;
            if(solver_select_working_set(solver,&i,&j)!=0)
                break;
            else
                counter = 1;	/* do shrinking next iteration */
        }

        ++iter;

        /* update alpha[i] and alpha[j], handle bounds carefully  */

        Q_i = kernel_get_Q(solver->Q,i,solver->active_size);
        Q_j = kernel_get_Q(solver->Q,j,solver->active_size);

        C_i = solver_get_C(solver,i);
        C_j = solver_get_C(solver,j);

        old_alpha_i = solver->alpha[i];
        old_alpha_j = solver->alpha[j];

        if(solver->y[i]!=solver->y[j])
        {
            double delta, diff;
            double quad_coef = solver->QD[i]+solver->QD[j]+2*Q_i[j];

            if (quad_coef <= 0)
                quad_coef = TAU;
            delta = (-solver->G[i]-solver->G[j])/quad_coef;
            diff = solver->alpha[i] - solver->alpha[j];
            solver->alpha[i] += delta;
            solver->alpha[j] += delta;

            if(diff > 0)
            {
                if(solver->alpha[j] < 0)
                {
                    solver->alpha[j] = 0;
                    solver->alpha[i] = diff;
                }
            }
            else
            {
                if(solver->alpha[i] < 0)
                {
                    solver->alpha[i] = 0;
                    solver->alpha[j] = -diff;
                }
            }
            if(diff > C_i - C_j)
            {
                if(solver->alpha[i] > C_i)
                {
                    solver->alpha[i] = C_i;
                    solver->alpha[j] = C_i - diff;
                }
            }
            else
            {
                if(solver->alpha[j] > C_j)
                {
                    solver->alpha[j] = C_j;
                    solver->alpha[i] = C_j + diff;
                }
            }
        }
        else
        {
            double delta,sum;
            double quad_coef = solver->QD[i]+solver->QD[j]-2*Q_i[j];
            if (quad_coef <= 0)
                quad_coef = TAU;
            delta = (solver->G[i]-solver->G[j])/quad_coef;
            sum = solver->alpha[i] + solver->alpha[j];
            solver->alpha[i] -= delta;
            solver->alpha[j] += delta;

            if(sum > C_i)
            {
                if(solver->alpha[i] > C_i)
                {
                    solver->alpha[i] = C_i;
                    solver->alpha[j] = sum - C_i;
                }
            }
            else
            {
                if(solver->alpha[j] < 0)
                {
                    solver->alpha[j] = 0;
                    solver->alpha[i] = sum;
                }
            }
            if(sum > C_j)
            {
                if(solver->alpha[j] > C_j)
                {
                    solver->alpha[j] = C_j;
                    solver->alpha[i] = sum - C_j;
                }
            }
            else
            {
                if(solver->alpha[i] < 0)
                {
                    solver->alpha[i] = 0;
                    solver->alpha[j] = sum;
                }
            }
        }

        /* update G */

        delta_alpha_i = solver->alpha[i] - old_alpha_i;
        delta_alpha_j = solver->alpha[j] - old_alpha_j;
        {
            int k;
            for( k=0;k<solver->active_size;k++)
            {
                solver->G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
            }
        }

        /* update alpha_status and G_bar */

        {
            int k;
            bool ui = solver_is_upper_bound(solver,i);
            bool uj = solver_is_upper_bound(solver,j);
            solver_update_alpha_status(solver,i);
            solver_update_alpha_status(solver,j);
            if(ui != solver_is_upper_bound(solver,i))
            {
                Q_i = kernel_get_Q(solver->Q,i,l);
                if(ui)
                    for(k=0;k<l;k++)
                        solver->G_bar[k] -= C_i * Q_i[k];
                else
                    for(k=0;k<l;k++)
                        solver->G_bar[k] += C_i * Q_i[k];
            }

            if(uj != solver_is_upper_bound(solver,j))
            {
                Q_j = kernel_get_Q(solver->Q,j,l);
                if(uj)
                    for(k=0;k<l;k++)
                        solver->G_bar[k] -= C_j * Q_j[k];
                else
                    for(k=0;k<l;k++)
                        solver->G_bar[k] += C_j * Q_j[k];
            }
        }
    }

    if(iter >= max_iter)
    {
        if(solver->active_size < l)
        {
            /* reconstruct the whole gradient to calculate objective value */
            solver_reconstruct_gradient(solver);
            solver->active_size = l;
        }
        fprintf(stderr,"\nWARNING: reaching max number of iterations\n");
    }

    /* calculate rho */

    si->rho = solver_calculate_rho(solver);

    /* calculate objective value */
    {
        double v = 0;
        int i;
        for(i=0;i<l;i++)
            v += solver->alpha[i] * (solver->G[i] + solver->p[i]);

        si->obj = v/2;
    }

    /* put back the solution */
    {
        int i;
        for(i=0;i<l;i++)
            alpha_[solver->active_set[i]] = solver->alpha[i];
    }


    si->upper_bound_p = Cp;
    si->upper_bound_n = Cn;

    free(solver->p);
    free(solver->y);
    free(solver->alpha);
    free(solver->alpha_status);
    free(solver->active_set);
    free(solver->G);
    free(solver->G_bar);
    delete_kernel(&(solver->Q));
}




/**
 * @brief solve_c_svc construct and solve various formulations
 * @param prob
 * @param param
 * @param alpha
 * @param si
 * @param Cp
 * @param Cn
 */
void solve_c_svc(
        svm_problem *prob,  svm_parameter* param,
        double *alpha,  SolutionInfo* si, double Cp, double Cn)
{
    int l = prob->l;
    double *minus_ones = (double*) malloc(sizeof(double)*l);
    schar *y = (schar*) malloc(sizeof(schar)*l);

    int i;
    Solver * solver;
    double sum_alpha;

    for(i=0;i<l;i++)
    {
        alpha[i] = 0;
        minus_ones[i] = -1;
        if(prob->y[i] > 0) y[i] = +1; else y[i] = -1;
    }
    solver= create_solver();

    solver_solve(solver,l,prob,param, minus_ones, y,
                 alpha, Cp, Cn, param->eps, si, param->shrinking);


    sum_alpha=0;
    for(i=0;i<l;i++)
        sum_alpha += alpha[i];
    if (Cp==Cn)
        for(i=0;i<l;i++)
            alpha[i] *= y[i];


    free(minus_ones);
    free(y);
    delete_solver(&solver);
}


typedef struct {
    double *alpha;
    double rho;
}decision_function;

decision_function svm_train_one(
        svm_problem *prob,  svm_parameter *param,
        double Cp, double Cn)
{
    double *alpha = (double *) malloc(sizeof(double)*prob->l);
    SolutionInfo si;
    int nSV = 0;
    int nBSV = 0;
    int i;
    decision_function f;
    solve_c_svc(prob,param,alpha,&si,Cp,Cn);


    /* output SVs */

    for( i=0;i<prob->l;i++)
    {
        alpha[i];

        if(fabs(alpha[i]) > 0)
        {
            ++nSV;
            if(prob->y[i] > 0)
            {
                if(fabs(alpha[i]) >= si.upper_bound_p)
                    ++nBSV;
            }
            else
            {
                if(fabs(alpha[i]) >= si.upper_bound_n)
                    ++nBSV;
            }
        }
    }

    f.alpha = alpha;
    f.rho = si.rho;
    return f;
}



/**
 * @brief multiclass_probability Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
 * @param k
 * @param r
 * @param p
 */
void multiclass_probability(int k, double **r, double *p)
{
    int t,j;
    int iter = 0, max_iter=max(100,k);
    double **Q=(double**)malloc(sizeof(double *)*k);
    double *Qp=(double*)malloc(sizeof(double*)*k);

    double pQp, eps=0.005/k;

    for (t=0;t<k;t++){
        p[t]=1.0/k;  /* Valid if k = 1 */
        Q[t]=(double*)malloc(sizeof(double)*k);
        Q[t][t]=0;
        for (j=0;j<t;j++)
        {
            Q[t][t]+=r[j][t]*r[j][t];
            Q[t][j]=Q[j][t];
        }
        for (j=t+1;j<k;j++)
        {
            Q[t][t]+=r[j][t]*r[j][t];
            Q[t][j]=-r[j][t]*r[t][j];
        }
    }
    for (iter=0;iter<max_iter;iter++)
    {
        double max_error;
        /* stopping condition, recalculate QP,pQP for numerical accuracy */
        pQp=0;
        for (t=0;t<k;t++)
        {
            Qp[t]=0;
            for (j=0;j<k;j++)
                Qp[t]+=Q[t][j]*p[j];
            pQp+=p[t]*Qp[t];
        }
        max_error=0;
        for (t=0;t<k;t++)
        {
            double error=fabs(Qp[t]-pQp);
            if (error>max_error)
                max_error=error;
        }
        if (max_error<eps) break;

        for (t=0;t<k;t++)
        {
            double diff=(-Qp[t]+pQp)/Q[t][t];
            p[t]+=diff;
            pQp=(pQp+diff*(diff*Q[t][t]+2*Qp[t]))/(1+diff)/(1+diff);
            for (j=0;j<k;j++)
            {
                Qp[j]=(Qp[j]+diff*Q[t][j])/(1+diff);
                p[j]/=(1+diff);
            }
        }
    }
    if (iter>=max_iter)
        for(t=0;t<k;t++) free(Q[t]);
    free(Q);
    free(Qp);
}


/**
 * @brief svm_group_classes
 * label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
 * perm, length l, must be allocated before calling this subroutine
 * @param prob
 * @param nr_class_ret
 * @param label_ret
 * @param start_ret
 * @param count_ret
 * @param perm
 */
void svm_group_classes(const svm_problem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm)
{
    int l = prob->l;
    int max_nr_class = 16;
    int nr_class = 0;
    int *label = (int*)malloc(sizeof(int)*max_nr_class);
    int *count = (int*)malloc(sizeof(int)*max_nr_class);
    int *data_label = (int*)malloc(sizeof(int)*l);
    int i;
    int * start;

    for(i=0;i<l;i++)
    {
        int this_label = (int)prob->y[i];
        int j;
        for(j=0;j<nr_class;j++)
        {
            if(this_label == label[j])
            {
                ++count[j];
                break;
            }
        }
        data_label[i] = j;
        if(j == nr_class)
        {
            if(nr_class == max_nr_class)
            {
                max_nr_class *= 2;
                label = (int *)realloc(label,max_nr_class*sizeof(int));
                count = (int *)realloc(count,max_nr_class*sizeof(int));
            }
            label[nr_class] = this_label;
            count[nr_class] = 1;
            ++nr_class;
        }
    }


    /* Labels are ordered by their first occurrence in the training set.
     However, for two-class sets with -1/+1 labels and -1 appears first,
     we swap labels to ensure that internally the binary SVM has positive data corresponding to the +1 instances.
    */
    if (nr_class == 2 && label[0] == -1 && label[1] == 1)
    {
        swap_int(&label[0],&label[1]);
        swap_int(&count[0],&count[1]);
        for(i=0;i<l;i++)
        {
            if(data_label[i] == 0)
                data_label[i] = 1;
            else
                data_label[i] = 0;
        }
    }

    assert(nr_class > 0);
    start = (int*)malloc(sizeof(int)*nr_class);
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+count[i-1];
    for(i=0;i<l;i++)
    {
        perm[start[data_label[i]]] = i;
        ++start[data_label[i]];
    }
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+count[i-1];

    *nr_class_ret = nr_class;
    *label_ret = label;
    *start_ret = start;
    *count_ret = count;
    free(data_label);
}

svm_model *svm_train( svm_problem *prob,  svm_parameter *param)
{
    svm_node ** x;
    int i,p;
    double * probA,*probB;
    double *weighted_C ;
    bool *nonzero;
    decision_function *f;
    int total_sv;
    int *nz_count;

    /* classification */
    int l = prob->l;
    int nr_class;
    int *label = NULL;
    int *start = NULL;
    int *count = NULL;
    int *perm = (int*)malloc(sizeof(int)*l);
    int *nz_start;

    svm_model *model = (svm_model*)malloc(sizeof(svm_model)*1);
    model->param = *param;
    model->free_sv = 0;




    /* group training data of the same class */
    svm_group_classes(prob,&nr_class,&label,&start,&count,perm);
    if(nr_class == 1) printf("WARNING: training data in only one class. See README for details.\n");

    x = (svm_node **)malloc(sizeof(svm_node *)*l);
    i;
    for(i=0;i<l;i++){
        x[i] = prob->x[perm[i]];
    }

    /* calculate weighted C */
    weighted_C = (double*)malloc(sizeof(double)* nr_class);
    for(i=0;i<nr_class;i++){
        weighted_C[i] = param->C;
    }
    for(i=0;i<param->nr_weight;i++){
        int j;
        for(j=0;j<nr_class;j++)
            if(param->weight_label[i] == label[j])
                break;
        if(j == nr_class)
            fprintf(stderr,"WARNING: class label %d specified in weight is not found\n", param->weight_label[i]);
        else
            weighted_C[j] *= param->weight[i];
    }

    /* train k*(k-1)/2 models */
    nonzero = (bool*)malloc(sizeof(bool)*l);
    for(i=0;i<l;i++)
        nonzero[i] = false;
    f = (decision_function *)malloc(sizeof(decision_function)*nr_class*(nr_class-1)/2);
    probA=NULL;
    probB=NULL;

    p = 0;
    for(i=0;i<nr_class;i++){
        int j;
        for( j=i+1;j<nr_class;j++){
            int k;
            svm_problem sub_prob;
            int si = start[i], sj = start[j];
            int ci = count[i], cj = count[j];
            sub_prob.l = ci+cj;
            sub_prob.x = (svm_node **)malloc(sizeof(svm_node *)*sub_prob.l);
            sub_prob.y = (double*)malloc(sizeof(double)*sub_prob.l);
            for(k=0;k<ci;k++){
                sub_prob.x[k] = x[si+k];
                sub_prob.y[k] = +1;
            }
            for(k=0;k<cj;k++){
                sub_prob.x[ci+k] = x[sj+k];
                sub_prob.y[ci+k] = -1;
            }
            f[p] = svm_train_one(&sub_prob,param,weighted_C[i],weighted_C[j]);


            for(k=0;k<ci;k++)
                if(!nonzero[si+k] && fabs(f[p].alpha[k]) > 0)
                    nonzero[si+k] = true;
            for(k=0;k<cj;k++)
                if(!nonzero[sj+k] && fabs(f[p].alpha[ci+k]) > 0)
                    nonzero[sj+k] = true;
            free(sub_prob.x);
            free(sub_prob.y);
            ++p;

        }
    }


    /* build output */
    model->nr_class = nr_class;

    model->label = (int*)malloc(sizeof(int)*nr_class);
    for(i=0;i<nr_class;i++)
        model->label[i] = label[i];

    model->rho = (double*)malloc(sizeof(double)*nr_class*(nr_class-1)/2);
    for(i=0;i<nr_class*(nr_class-1)/2;i++)
        model->rho[i] = f[i].rho;


    model->probA=NULL;
    model->probB=NULL;

    total_sv = 0;
    nz_count = (int*)malloc(sizeof(int)*nr_class);
    model->nSV = (int*)malloc(sizeof(int)*nr_class);
    for(i=0;i<nr_class;i++)
    {
        int nSV = 0;
        int j;
        for( j=0;j<count[i];j++)
            if(nonzero[start[i]+j])
            {
                ++nSV;
                ++total_sv;
            }
        model->nSV[i] = nSV;
        nz_count[i] = nSV;
    }

    model->l = total_sv;
    model->SV = (svm_node **)malloc(sizeof(svm_node *)*total_sv);
    model->sv_indices = (int*)malloc(sizeof(int)*total_sv);
    p = 0;
    for(i=0;i<l;i++)
        if(nonzero[i])
        {
            model->SV[p] = x[i];
            model->sv_indices[p++] = perm[i] + 1;
        }

    nz_start = (int*)malloc(sizeof(int)*nr_class);
    nz_start[0] = 0;
    for(i=1;i<nr_class;i++)
        nz_start[i] = nz_start[i-1]+nz_count[i-1];

    model->sv_coef = (double**)malloc(sizeof(double*)*nr_class-1);
    for(i=0;i<nr_class-1;i++)
        model->sv_coef[i] = (double*)malloc(sizeof(double)*total_sv);

    p = 0;
    for(i=0;i<nr_class;i++){
        int j;
        for( j=i+1;j<nr_class;j++)
        {
            /* classifier (i,j): coefficients with
               i are in sv_coef[j-1][nz_start[i]...],
               j are in sv_coef[i][nz_start[j]...] */

            int si = start[i];
            int sj = start[j];
            int ci = count[i];
            int cj = count[j];

            int q = nz_start[i];
            int k;
            for(k=0;k<ci;k++)
                if(nonzero[si+k])
                    model->sv_coef[j-1][q++] = f[p].alpha[k];
            q = nz_start[j];
            for(k=0;k<cj;k++)
                if(nonzero[sj+k])
                    model->sv_coef[i][q++] = f[p].alpha[ci+k];
            ++p;
        }
    }


    free(label);
    free(probA);
    free(probB);
    free(count);
    free(perm);
    free(start);
    free(x);
    free(weighted_C);
    free(nonzero);
    for(i=0;i<nr_class*(nr_class-1)/2;i++)
        free(f[i].alpha);
    free(f);
    free(nz_count);
    free(nz_start);


    return model;
}



int svm_get_nr_class( svm_model *model)
{
    return model->nr_class;
}

void svm_get_labels( svm_model *model, int* label)
{
    int i;
    if (model->label != NULL)
        for( i=0;i<model->nr_class;i++)
            label[i] = model->label[i];
}

void svm_get_sv_indices( svm_model *model, int* indices)
{
    int i;
    if (model->sv_indices != NULL)
        for( i=0;i<model->l;i++)
            indices[i] = model->sv_indices[i];
}

int svm_get_nr_sv(const svm_model *model)
{
    return model->l;
}

double svm_predict_values( svm_model *model,  svm_node *x, double* dec_values)
{
    int i,p,vote_max_idx;

    int nr_class = model->nr_class;
    int l = model->l;
    int * start,*vote;

    double *kvalue = (double*)malloc(sizeof(double)*l);
    for(i=0;i<l;i++)
        kvalue[i] = kernel_k_function(x,model->SV[i],&model->param);

    start = (int*)malloc(sizeof(int)*nr_class);
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+model->nSV[i-1];

    vote = (int*)malloc(sizeof(int)*nr_class);
    for(i=0;i<nr_class;i++)
        vote[i] = 0;

    p=0;
    for(i=0;i<nr_class;i++){
        int j;
        for( j=i+1;j<nr_class;j++)
        {
            double sum = 0;
            int si = start[i];
            int sj = start[j];
            int ci = model->nSV[i];
            int cj = model->nSV[j];

            int k;
            double *coef1 = model->sv_coef[j-1];
            double *coef2 = model->sv_coef[i];
            for(k=0;k<ci;k++)
                sum += coef1[si+k] * kvalue[si+k];
            for(k=0;k<cj;k++)
                sum += coef2[sj+k] * kvalue[sj+k];
            sum -= model->rho[p];
            dec_values[p] = sum;

            if(dec_values[p] > 0)
                ++vote[i];
            else
                ++vote[j];
            p++;
        }
    }

    vote_max_idx = 0;
    for(i=1;i<nr_class;i++)
        if(vote[i] > vote[vote_max_idx])
            vote_max_idx = i;

    free(kvalue);
    free(start);
    free(vote);
    return model->label[vote_max_idx];

}

double svm_predict( svm_model *model,  svm_node *x)
{
    int nr_class = model->nr_class;
    double *dec_values;
    double pred_result;
    dec_values = (double*)malloc(sizeof(double)* nr_class*(nr_class-1)/2);
    pred_result = svm_predict_values(model, x, dec_values);
    free(dec_values);
    return pred_result;
}


void svm_free_model_content(svm_model* model_ptr)
{
    int i;
    if(model_ptr->free_sv && model_ptr->l > 0 && model_ptr->SV != NULL)
        free((void *)(model_ptr->SV[0]));
    if(model_ptr->sv_coef)
    {
        for( i=0;i<model_ptr->nr_class-1;i++)
            free(model_ptr->sv_coef[i]);
    }

    free(model_ptr->SV);
    model_ptr->SV = NULL;

    free(model_ptr->sv_coef);
    model_ptr->sv_coef = NULL;

    free(model_ptr->rho);
    model_ptr->rho = NULL;

    free(model_ptr->label);
    model_ptr->label= NULL;

    free(model_ptr->probA);
    model_ptr->probA = NULL;

    free(model_ptr->probB);
    model_ptr->probB= NULL;

    free(model_ptr->sv_indices);
    model_ptr->sv_indices = NULL;

    free(model_ptr->nSV);
    model_ptr->nSV = NULL;
}

void svm_free_and_destroy_model(svm_model** model_ptr_ptr)
{
    if(model_ptr_ptr != NULL && *model_ptr_ptr != NULL)
    {
        svm_free_model_content(*model_ptr_ptr);
        free(*model_ptr_ptr);
        *model_ptr_ptr = NULL;
    }
}

void svm_destroy_param(svm_parameter* param)
{
    free(param->weight_label);
    free(param->weight);
}

char *svm_check_parameter( svm_problem *prob,  svm_parameter *param)
{

    if(param->gamma < 0)
        return "gamma < 0";

    if(param->degree < 0)
        return "degree of polynomial kernel < 0";


    if(param->cache_size <= 0)
        return "cache_size <= 0";

    if(param->eps <= 0)
        return "eps <= 0";

    if(param->C <= 0)
        return "C <= 0";


    if(param->shrinking != 0 &&
            param->shrinking != 1)
        return "shrinking != 0 and shrinking != 1";



    return NULL;
}



void svm_set_print_string_function(void (*print_func)(const char *))
{
    if(print_func == NULL)
        svm_print_string = &print_string_stdout;
    else
        svm_print_string = print_func;
}



/**
 * @brief SCIS implementation : first part
 * Author : Bérengère MATHIEU
 *
 * Copyright (c) Bérengère MATHIEU
 *
 * Under CeCILL-C licence
 *
 **/


/**
 * @brief to store color in Lab color space
 **/
typedef struct{
    double l;
    double a;
    double b;
}LAB_Color;

/**
 * @brief to store color in RGB  color space
 **/
typedef struct{
    double red;
    double green;
    double blue;
}RGB_Color;

/**
 * @brief simple 2D Matrix of double
 **/
typedef struct{
    double ** data;
    unsigned int height;
    unsigned int width;
} _Matrix;

typedef _Matrix* Matrix;

/**
 * @brief allocate memory for 2D Matrix but not initialize it
 **/
Matrix create_matrix(unsigned int height,unsigned int width){
    int y;
    Matrix res= (Matrix) malloc(sizeof(_Matrix));
    res->height = height;
    res->width = width;

    res->data = (double ** ) malloc(sizeof(double*)*height);

    for(y=0;y<height;y++){
        res->data[y] = (double*) malloc(sizeof(double)*width);
    }
    return res;
}

/**
 * @brief clean up memory
 **/
void delete_matrix(Matrix* matrix){
    int y;
    for(y=0;y<(*matrix)->height;y++){
        free((*matrix)->data[y]);
    }

    free((*matrix)->data);

    free(*matrix);

}

/**
 * To store a RGB image
 **/
typedef struct{
    unsigned char** red;
    unsigned char** green;
    unsigned char** blue;
    unsigned int nbRow;
    unsigned int nbCol;
}_Image;
typedef _Image* Image;

/**
 * @brief allocate memory but not initialize Image data
 **/
Image create_Image(unsigned int nbRow,unsigned int nbCol){
    Image im=(Image)malloc(sizeof(_Image));
    int r;
    im->nbCol=nbCol;
    im->nbRow=nbRow;

    im->red=(unsigned char**)malloc(sizeof(unsigned char*)*nbRow);
    im->green=(unsigned char**)malloc(sizeof(unsigned char*)*nbRow);
    im->blue=(unsigned char**)malloc(sizeof(unsigned char*)*nbRow);

    for(r=0;r<nbRow;r++){
        im->red[r]=(unsigned char*)malloc(sizeof(unsigned char)*nbCol);
        im->green[r]=(unsigned char*)malloc(sizeof(unsigned char)*nbCol);
        im->blue[r]=(unsigned char*)malloc(sizeof(unsigned char)*nbCol);

    }

    return im;
}

/**
 * @brief clean up memory
 **/
void delete_Image(Image* im){
    int r;

    for(r=0;r<(*im)->nbRow;r++){
        free((*im)->red[r]) ;
        free((*im)->green[r]);
        free((*im)->blue[r]);

    }

    free((*im)->red);
    free((*im)->green);
    free((*im)->blue);

    free(*im);

}

/**
 * Copy an image
 **/
Image ImCopy(Image im){
    int r;

    Image res=create_Image(im->nbRow,im->nbCol);
    for(r=0;r<im->nbRow;r++){
        int c;
        for(c=0;c<im->nbCol;c++){
            res->red[r][c]=im->red[r][c];
            res->green[r][c]=im->green[r][c];
            res->blue[r][c]=im->blue[r][c];
        }
    }


    return res;
}

/**
 * @brief convert color from Lab color space to rgb color space
 **/
RGB_Color convert_lab_to_rgb(LAB_Color lab){
    double r,g,b;

    double x,y,z;

    double beta,gamma,epsilon,kappa,Xr,Yr,Zr;

    RGB_Color res;

    /* convert from Lab color space to XYZ color space
       http://en.wikipedia.org/wiki/Lab_color_space */
    double alpha=1./116.;
    y=alpha*(lab.l+16.);
    x=y+(1./500.)*lab.a;
    z=y-(1./200.)*lab.b;

    alpha=6./29.;
    beta=3.*pow(alpha,2);
    gamma=4./29.;

    if(y>alpha) y = pow(y,3);
    else y = beta*(y-gamma);

    if(x>alpha) x = pow(x,3);
    else x = beta*(x-gamma);

    if(z>alpha) z = pow(z,3);
    else z = beta*(z-gamma);

    /* scale x,y,z values according to reference white (D65)*/
    epsilon = 0.008856;

    kappa   = 903.3;

    Xr = 0.9504; /* x value for reference white point D65 */

    Yr = 1.0; /* y value for reference white point D65 */

    Zr = 1.0889; /* z value for reference white point D65 */
    /* scale according to reference white (D65) */
    y=y*Yr;
    x=x*Xr;
    z=z*Zr;

    /* convert XYZ to normilized sRGB linear */
    r=3.240479*x -1.537150*y -0.498535*z;
    g=-0.969256*x + 1.875992*y + 0.041556*z;
    b= 0.055648*x -0.204043*y +  1.057311*z;
    /* get gamma corrected normalized sRGB values */
    if ( r > 0.0031308 ) r = 1.055 * ( pow(r, ( 1 / 2.4 )) ) - 0.055;
    else                     r = 12.92 * r;
    if ( g > 0.0031308 ) g = 1.055 * ( pow(g,( 1 / 2.4 )) ) - 0.055;
    else                     g = 12.92 * g;
    if ( b > 0.0031308 ) b = 1.055 * ( pow(b, ( 1 / 2.4 ) )) - 0.055;
    else                     b = 12.92 * b;

    /* scale normalized sRGB values from 0 to 255 */
    res.red=r*255;
    res.green=g*255;
    res.blue=b*255;

    return res;

}

/**
 * @brief convert color from RGB color space to Lab color space
 **/
LAB_Color convert_rgb_to_lab(unsigned char red,unsigned char green,unsigned blue){
    /* assuming we are in sRGB color space
       with gamma corrected sRGB values
       in range [0,255]
       http://en.wikipedia.org/wiki/SRGB
      use CIE Standard Illuminant D65 */

    double beta,gamma,epsilon,kappa,Xr,Yr,Zr,xr,yr,zr, fx, fy, fz;

    LAB_Color res;

    /*convert from sRGB to XYZ color space */
    double x, y, z;

    /* normalize value of red, green */
    double r = (double)red/255.0;
    double g = (double)green/255.0;
    double b = (double)blue/255.0;

    /* convert gamma corrected sRGB values to linear sRGB values */
    if(r <= 0.04045)	r = r/12.92;
    else				r = pow((r+0.055)/1.055,2.4);
    if(g <= 0.04045)	g = g/12.92;
    else				g = pow((g+0.055)/1.055,2.4);
    if(b <= 0.04045)	b = b/12.92;
    else				b = pow((b+0.055)/1.055,2.4);

    /* convert from linear sRGB values to XYZ value */
    x = r*0.4124564 + g*0.3575761 + b*0.1804375;
    y = r*0.2126729 + g*0.7151522 + b*0.0721750;
    z = r*0.0193339 + g*0.1191920 + b*0.9503041;

    /* scale x,y,z values according to reference white (D65) */
    epsilon = 0.008856;

    kappa   = 903.3;

    Xr = 0.9504; /* x value for reference white point D65 */

    Yr = 1.0; /* y value for reference white point D65 */

    Zr = 1.0889; /* z value for reference white point D65 */
    xr = x/Xr;
    yr = y/Yr;
    zr = z/Zr;

    /* convert from XYZ color space to Lab colorspace
    http://en.wikipedia.org/wiki/Lab_color_space */

    if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
    else				fx = (kappa*xr + 16.0)/116.0;

    if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
    else				fy = (kappa*yr + 16.0)/116.0;

    if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
    else				fz = (kappa*zr + 16.0)/116.0;

    res.l = 116.0*fy-16.0;
    res.a = 500.0*(fx-fy);
    res.b = 200.0*(fy-fz);

    return res;

}


/**
 * @brief extract color corresponding to the various semantic classes defined by user
 *
 * Simply return an Image with one row by color present in the image with the seeds given by the user
 **/
Image extractClassColors(Image seeds){

    int height = seeds->nbRow;
    int width = seeds->nbCol;

    /* allocate memory for 10 colors */
    Image classColors=create_Image(100,1);
    RGB_Color pixelColor;
    bool knownClass;
    int nextColor=0;
    int r;


    /* search all seeds colors */
    for( r=0;r<height;++r){
        int c;
        for( c=0;c<width;++c){
            pixelColor.red=seeds->red[r][c];
            pixelColor.green=seeds->green[r][c];
            pixelColor.blue=seeds->blue[r][c];
            knownClass=false;
            if(pixelColor.red>0 || pixelColor.green>0 || pixelColor.blue>0){ /* not black */
                int i;
                /* check if color is known */
                for( i=0;i<nextColor;++i){
                    if(pixelColor.red==classColors->red[i][0] &&
                            pixelColor.green==classColors->green[i][0] &&
                            pixelColor.blue==classColors->blue[i][0] ){
                        knownClass=true;
                        break;
                    }
                }

                if(!knownClass){
                    /* add new color */
                    classColors->red[nextColor][0]=pixelColor.red;
                    classColors->green[nextColor][0]=pixelColor.green;
                    classColors->blue[nextColor][0]=pixelColor.blue;
                    nextColor++;
                }
                /*too much different colors found */
                if(nextColor>=classColors->nbRow){
                    STATE=2;
                    errror_handler();
                    return NULL;
                }
            }
        }


    }

    /* remove useless rows */
    if(classColors->nbRow!=nextColor){
        int j;
        Image temp=ImCopy(classColors);
        delete_Image(&classColors);
        classColors=create_Image(nextColor,1);

        for( j=0;j<nextColor;j++){
            classColors->red[j][0]=temp->red[j][0];
            classColors->green[j][0]=temp->green[j][0];
            classColors->blue[j][0]=temp->blue[j][0];
        }

        delete_Image(&temp);
    }



    return classColors;

}




/************************************************************************************************************************
 *
 * Superpixel segmentation : second part
 *
 * Adapted from C++ implementation of Pedro Felzenzswalb
 * Copyright (c) Pedro Felzenszwalb
 *
 * http://cs.brown.edu/~pff/segment/
 *
 ************************************************************************************************************************/
/**
 * Compute distance between two pixels according to they RGB colors
 **/
float felzenszwalbAlgorithm_diff(double **r, double **g, double **b, int x1, int y1, int x2, int y2){

    return sqrt(pow(r[y1][x1]- r[y2][x2],2) +
                pow(g[y1][x1]- g[y2][x2],2) +
                pow(b[y1][x1]- b[y2][x2],2) );
}


/**
 * @brief compute image over-segmentation
 **/
Matrix felzenszwalbAlgorithm_segment(float c,int minSize,Image image, int* nbSuperpixels, bool gimpPlugin){
    DisjointSet *u ;
    int i;
    GraphEdge* edges ;
    int num;

    int width = image->nbCol;
    int height = image->nbRow;
    int nbPixels = width*height;


    /* split color channels */
    Matrix red = create_matrix(height,width);
    Matrix green = create_matrix(height,width);
    Matrix blue = create_matrix(height,width);


    /* to update gimp progression bar */
    double totalActions=2.*height;
    double doneActions=0;
    int r;

    Matrix res;

    for( r=0;r<height;++r){
        int c;

        for(c=0;c<width;++c){
            red->data[r][c] = image->red[r][c];
            green->data[r][c] = image->green[r][c];
            blue->data[r][c] = image->blue[r][c];
        }
        doneActions++;
        if(gimpPlugin)  gimp_progress_update(doneActions/totalActions);


    }


    edges = (GraphEdge*) malloc(sizeof(GraphEdge)*nbPixels*4);
    num =0;
    for(r=0;r<height;++r){
        int c;
        for(c=0;c<width;++c){

            if(c<width-1){
                /* pixel has a east neighbor */
                edges[num].subSet1 = r*width+c;
                edges[num].subSet2 = r*width+(c+1);
                edges[num].w = felzenszwalbAlgorithm_diff(red->data, green->data, blue->data, c, r, c+1, r);
                num++;

            }

            if(r<height-1) {
                /* pixel has a south neighbor */
                edges[num].subSet1 = r*width+c;
                edges[num].subSet2 = (r+1)*width+c;
                edges[num].w = felzenszwalbAlgorithm_diff(red->data, green->data, blue->data, c, r, c, r+1);
                num++;

            }

            if((c<width-1) && (r<height-1)){
                /* pixel has a south-east neighbor */
                edges[num].subSet1 = r*width+c;
                edges[num].subSet2 = (r+1)*width+(c+1);
                edges[num].w = felzenszwalbAlgorithm_diff(red->data, green->data, blue->data, c, r, c+1, r+1);
                num++;

            }

            if((c<width-1) && (r>0)){
                /* pixel has a north-east neighbor */
                edges[num].subSet1 = r*width+c;
                edges[num].subSet2 = (r-1)*width+(c+1);
                edges[num].w = felzenszwalbAlgorithm_diff(red->data, green->data, blue->data, c, r, c+1, r-1);
                num++;

            }
        }
        doneActions++;
        if(gimpPlugin)  gimp_progress_update(doneActions/totalActions);
    }

    delete_matrix(&red);
    delete_matrix(&green);
    delete_matrix(&blue);


    u = felzenszwalb_oversegmentation(nbPixels, num, edges, c);

    /* post process small components */
    totalActions=num+2*height;
    doneActions=0;
    i;
    for (i = 0; i < num; i++) {
        int a = disjointSet_find(u,edges[i].subSet1);
        int b = disjointSet_find(u,edges[i].subSet2);
        if ((a != b) && ((get_subset_size(u,a) < minSize) || (get_subset_size(u,b) < minSize)))
            disjointSet_join(u,a, b);

        doneActions++;
        if(gimpPlugin)  gimp_progress_update(doneActions/totalActions);
    }
    free(edges);
    *nbSuperpixels = get_nb_subsets(u);

    /* ensure that superpixels labels are consecutive integers */
    {
        int label=0;
        int* labels = (int* ) malloc(sizeof(int)*(*nbSuperpixels));
        for(r=0;r<height;++r){
            int c;
            for(c=0;c<width;++c){
                int k;
                int comp = disjointSet_find(u,r * width + c);
                bool find=false;
                for(k=0;k<label;k++){
                    if(labels[k]==comp){
                        find=true;
                        break;
                    }
                }
                if(!find){
                    labels[label]=comp;
                    ++label;
                }
            }
            doneActions++;
            if(gimpPlugin)  gimp_progress_update(doneActions/totalActions);
        }

        /* create a matrix containing for each pixel
       label of its superpixel */
        res = create_matrix(height,width);
        for( r=0;r<height;++r){
            int c;
            for( c=0;c<width;++c){
                int comp = disjointSet_find(u,r * width + c);
                int newLabel=0;
                while(comp!=labels[newLabel]) newLabel++;
                res->data[r][c] = newLabel;
            }
            doneActions++;
            if(gimpPlugin)  gimp_progress_update(doneActions/totalActions);

        }

        delete_disjoint_set(&u);
        free(labels);
        return res;
    }
}


/**
 * @brief SCIS implementation : second part
 * Author : Bérengère MATHIEU
 *
 * Copyright (c) Bérengère MATHIEU
 *
 * Under CeCILL-C licence
 *
 **/

/**
 * @class SVMAndSuperpixels
 * @brief Interactive segmentation of an image using SVM classifier.
 *
 * First image is oversegmented into superpixels using Felzenszwalb algorithm
 * Then average color in RGB color space and center of mass is computed for each superpixels
 * Next superpixels with pixels labeled by user are used to train SVM
 * Finally remaining superpixels are classified using SVM
 *
 */
typedef struct   {

    float gamma_; /**< kernel parameter */
    float C_;/**< regularization parameter */

    svm_model * model;/**< trained SVM*/

    svm_problem prob; /**< to store training data information */

    svm_node *x_space; /**< to store training data information*/

    int nbFeatures;/**< number of features used to describe a superpixel */
    int minSize_felzenswalb_;/**< minimum component size*/

    SuperpixelFeatures* trainingData_;/**< descripors of superpixels used for training classifier */
    int nbTrainingData;

    Image classColors;/**< link to each class a color, according to seeds given by user*/
} SVMAndSuperpixels;


/**
 * @brief train SVM classifier
 **/
svm_model * SVMAndSuperpixels_trainSVM(SVMAndSuperpixels* algo){


    svm_model * model ;
    svm_parameter param;
    param.degree = 3;
    param.gamma =algo->gamma_;	/* 1/num_features */
    param.cache_size = 100;
    param.C = algo->C_;
    param.eps = 1e-3;
    param.shrinking = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    /* initialize prob with training data */
    algo->prob.l=algo->nbTrainingData;
    /* to store class labels */
    algo->prob.y = (double*) malloc(sizeof(double)*algo->prob.l);
    /* features access */
    algo->prob.x=( svm_node**) malloc(sizeof( svm_node*)*algo->prob.l);
    /*to store all features and their index
      if a feature is equal to 0, it is not stored
      index of the feature is used to know each value corresponding to the same feature */
    algo->x_space=( svm_node*) malloc(sizeof( svm_node)*algo->prob.l*(algo->nbFeatures+1));

    {
        int j=0;/* index of the svm_node for store the next feature */
        int i=0;
        int index;
        for(index=0;index<algo->nbTrainingData;index++){
            int k=0;/* current feature index */
            /* store class label of superpixel */
            algo->prob.y[i]=algo->trainingData_[i].classLabel;
            /* store features of superpixel */
            algo->prob.x[i]=&(algo->x_space[j]);
            /* if a feature is equal to 0, it is not stored */
            if(algo->trainingData_[index].colorAtt1 !=0 ){ /* first color coordinate */
                algo->x_space[j].index=k+1;
                algo->x_space[j].value=algo->trainingData_[index].colorAtt1;
                ++j;
            }
            k++;
            if(algo->trainingData_[index].colorAtt2 !=0 ){ /* second color coordinate */
                algo->x_space[j].index=k+1;
                algo->x_space[j].value=algo->trainingData_[index].colorAtt2;
                ++j;
            }
            k++;
            if(algo->trainingData_[index].colorAtt3 !=0 ){/* third color coordinate */
                algo->x_space[j].index=k+1;
                algo->x_space[j].value=algo->trainingData_[index].colorAtt3;
                ++j;
            }
            k++;

            if(algo->trainingData_[index].r !=0 ){/* ordinate */
                algo->x_space[j].index=k+1;
                algo->x_space[j].value=algo->trainingData_[index].r;
                ++j;
            }
            k++;
            if(algo->trainingData_[index].c !=0 ){/* abscissa */
                algo->x_space[j].index=k+1;
                algo->x_space[j].value=algo->trainingData_[index].c;
                k++;
                ++j;
            }

            algo->x_space[j].index=-1;
            ++j;
            ++i;

        }
    }

    /* train SVM using features of superpixels */
    model = svm_train(&(algo->prob),&param);

    /* return trained SVM */
    return model;


}



/**
 * @brief compute for each supeprixel its features and sort superpixel between training superpixels (superpixel with seeds
 * for a unique class) and test superpixels.
 **/
void SVMAndSuperpixels_extractTrainingAndTestData(SVMAndSuperpixels* algo,Image image, Image seeds, Matrix superpixels, int nbSuperpixels, SuperpixelFeatures** data,int* nbData, bool gimpPlugin){
    /* get image size */
    int height = image->nbRow;
    int width = image->nbCol;
    int i,r,idxTraining;
    double totalActions, doneActions;

    /* Compute superpixels features
       used features are color and position
       to access to pixels labels */
    *nbData=nbSuperpixels;
    /* create a new vector to store each superpixel features */
    *data=(SuperpixelFeatures*)malloc(sizeof(SuperpixelFeatures)*nbSuperpixels);
    for( i=0;i<nbSuperpixels;i++){
        (*data)[i].colorAtt1=0;
        (*data)[i].colorAtt2=0;
        (*data)[i].colorAtt3=0;
        (*data)[i].c=0;
        (*data)[i].r=0;
        (*data)[i].nbPixels=0;
        (*data)[i].classLabel=-1;
    }

    /* to update gimp progress bar */
    totalActions=height+nbSuperpixels;
    doneActions=0;

    for(r=0;r<height;++r){
        int c;
        for(c=0;c<width;++c){
            int label = superpixels->data[r][c];/* get pixel label */
            int i;

            /* get pixel RGB color in the image
              update superpixel data by adding color and position of the current pixel */
            (*data)[label].colorAtt1 += (double) image->red[r][c];
            (*data)[label].colorAtt2 += (double) image->green[r][c];
            (*data)[label].colorAtt3 += (double) image->blue[r][c];
            (*data)[label].r+=r;
            (*data)[label].c+=c;
            /* update number of pixels which are part of the superxiel */
            (*data)[label].nbPixels++;

            /* we now check if the pixel is a seed
               get the color of the pixel, in the seeds image */


            /* if the pixel is a seed, its color is in the classColors */
            for( i=0;i<algo->classColors->nbRow;++i){
                if(seeds->red[r][c]>0 || seeds->green[r][c] >0 || seeds->blue[r][c]>0){
                    if(seeds->red[r][c]==algo->classColors->red[i][0]&&
                            seeds->green[r][c]==algo->classColors->green[i][0]&&
                            seeds->blue[r][c]==algo->classColors->blue[i][0]){
                        /* the pixel is a seed */
                        if((*data)[label].classLabel==-1){
                            /* first pixel of the superpixel, wich is a seed
                               give its label to the superpixel */
                            (*data)[label].classLabel=i;/* save the class label of the pixel */
                        }else{
                            /* some seeds are already part of the superpixel */
                            if((*data)[label].classLabel != i){
                                /* two pixels have different label
                                   superpixel won't be used to training SVM */
                                (*data)[label].classLabel = -2;
                            }
                        }
                        break;
                    }
                }
            }
        }
        doneActions++;
        if(gimpPlugin) gimp_progress_update(doneActions/totalActions);
    }
    /* compute number of training superpixel */
    algo->nbTrainingData=0;
    for(i=0;i<nbSuperpixels;++i){

        if((*data)[i].classLabel >= 0 ){
            /* add superpixel to the training set
               if it has the label of a class */
            algo->nbTrainingData++;
        }

    }



    algo->trainingData_=(SuperpixelFeatures*)malloc(sizeof(SuperpixelFeatures)*algo->nbTrainingData);
    /* for each superpixels :
       compute its average color
       compute its average position
       add it to the traning set, if superpixel add seeds of only one class */
    idxTraining=0;
    for(i=0;i<nbSuperpixels;++i){
        /* compute average RGB color */
        (*data)[i].colorAtt1 /= (double) (*data)[i].nbPixels;
        (*data)[i].colorAtt2 /= (double) (*data)[i].nbPixels;
        (*data)[i].colorAtt3 /= (double) (*data)[i].nbPixels;

        /* normalize RGB color */
        (*data)[i].colorAtt1 = (double) (*data)[i].colorAtt1/255.;
        (*data)[i].colorAtt2 = (double) (*data)[i].colorAtt2/255;
        (*data)[i].colorAtt3 = (double) (*data)[i].colorAtt3/255.;

        /* compute average ordinate */
        (*data)[i].r /= (*data)[i].nbPixels;
        /* normalize ordinate */
        (*data)[i].r /= height -1;
        /* compute average abscissa */
        (*data)[i].c /= (*data)[i].nbPixels;
        /* normalize abscissa */
        (*data)[i].c /= width -1;

        if((*data)[i].classLabel >= 0 ){
            /* add superpixel to the training set
               if it has the label of a class */
            algo->trainingData_[idxTraining].colorAtt1=(*data)[i].colorAtt1;
            algo->trainingData_[idxTraining].colorAtt2=(*data)[i].colorAtt2;
            algo->trainingData_[idxTraining].colorAtt3=(*data)[i].colorAtt3;
            algo->trainingData_[idxTraining].r=(*data)[i].r;
            algo->trainingData_[idxTraining].c=(*data)[i].c;
            algo->trainingData_[idxTraining].nbPixels=(*data)[i].nbPixels;
            algo->trainingData_[idxTraining].classLabel=(*data)[i].classLabel;

            idxTraining++;
        }


        doneActions++;
        if(gimpPlugin) gimp_progress_update(doneActions/totalActions);

    }



}


/**
 * @brief clean up memory
 **/
void SVMAndSuperpixels_clearData(SVMAndSuperpixels* algo){
    free(algo->trainingData_);
    delete_Image(&(algo->classColors));
    svm_free_and_destroy_model(&(algo->model));
}

/**
 * @brief Main method of SCIS : do interactive segmentation
 **/
Image SVMAndSuperpixels_segment(SVMAndSuperpixels* algo,Image seeds, Image image, bool gimpPlugin) {

    int nbSuperpixels,nbData;
    int height,width;
    Matrix superpixels ;
    SuperpixelFeatures * data;
    if(gimpPlugin) gimp_progress_init("extract classes ...");
    /* store class colors
       extract all colors (with the exception of black) in seeds picture */
    algo->classColors = extractClassColors(seeds);
    if(algo->classColors==NULL) return;
    if(algo->classColors->nbRow<2){
        STATE=1;
        errror_handler();
        return NULL;
    }

    if(gimpPlugin) gimp_progress_end();


    /* compute superpixels */
    if(gimpPlugin) gimp_progress_init("oversegment ...");

    superpixels = felzenszwalbAlgorithm_segment(24,20,image,&nbSuperpixels,gimpPlugin);



    if(gimpPlugin) gimp_progress_end();
    height = image->nbRow;
    width = image->nbCol;

    /* Extract training and test data */
    if(gimpPlugin) gimp_progress_init("extract superpixels features..");

    SVMAndSuperpixels_extractTrainingAndTestData(algo,image,seeds,superpixels,nbSuperpixels,&data,&nbData,gimpPlugin);
    if(gimpPlugin) gimp_progress_end();



    /* Train SVM */
    if(gimpPlugin) gimp_progress_init("train classifier ...");
    algo->model = SVMAndSuperpixels_trainSVM(algo);
    if(algo->model==NULL){
        STATE=5;
        errror_handler();
        return;
    }
    if(gimpPlugin) gimp_progress_end();

    /* Use model to predict segmentation */
    if(gimpPlugin) gimp_progress_init("classify superpixels ...");
    {
        double predict_label ;
        svm_node* x=(svm_node*)malloc(sizeof(svm_node)*(algo->nbFeatures+1));
        Image res;


        /* Class for each superpixel */
        int* spClass=(int*)malloc(sizeof(int)*nbSuperpixels);
        int index;
        int i;
        for( i=0;i<nbSuperpixels;i++){
            spClass[i]=0;
        }

        for( index=0;index<nbData;++index){/* for each superpixel */
            int j=0;/* svm_node index used to store features of the current superpixel */
            int k=0;/* feature index */

            /* features equal to 0 are not stored */
            if(data[index].colorAtt1 !=0 ){
                x[j].index=k+1;
                x[j].value=data[index].colorAtt1;
                ++j;
            }
            k++;

            if(data[index].colorAtt2 !=0 ){
                x[j].index=k+1;
                x[j].value=data[index].colorAtt2;
                ++j;
            }
            k++;

            if(data[index].colorAtt3 !=0 ){
                x[j].index=k+1;
                x[j].value=data[index].colorAtt3;
                ++j;
            }
            k++;

            if(data[index].r !=0 ){
                x[j].index=k+1;
                x[j].value=data[index].r;
                ++j;
            }
            k++;

            if(data[index].c !=0 ){
                x[j].index=k+1;
                x[j].value=data[index].c;
                ++j;
            }
            k++;
            x[j].index=-1;/* end of features */
            /* test if superpixel is a seed */
            if(data[index].classLabel >=0  ){
                /* superpixel is a seed : do not predict its label
               label given by user is considered as better */
                spClass[index]=data[index].classLabel;
            }else{
                /* predict class of superpixel */
                predict_label = svm_predict(algo->model,x);
                /* store class label */
                spClass[index]=predict_label;
            }
            if(gimpPlugin) gimp_progress_update((double)(index+1)/(double)nbData);

        }
        if(gimpPlugin) gimp_progress_end();


        /* draw all pixel of the same class, with the color given by user
      for this class (thanks to seeds) */

        /* Image to store segmentation result */
        res = create_Image(height,width);
        if(gimpPlugin) gimp_progress_init("classify pixels ...");
        {
            int r;
            for(r=0;r<height;++r){
                int c;
                for( c=0;c<width;++c){
                    int l=(int)superpixels->data[r][c];
                    res->red[r][c]=algo->classColors->red[spClass[l]][0];
                    res->green[r][c]=algo->classColors->green[spClass[l]][0];
                    res->blue[r][c]=algo->classColors->blue[spClass[l]][0];
                }
                if(gimpPlugin) gimp_progress_update((double)r+1/(double)height);
            }
        }
        if(gimpPlugin) gimp_progress_end();

        free(spClass);
        delete_matrix(&superpixels);
        free(data);
        free(algo->x_space);
        free(algo->prob.x);
        free(algo->prob.y);
        SVMAndSuperpixels_clearData(algo);
        free(x);

        return res;
    }
}


/**
 *  @brief Store image from gimp drawable to more simle Image data structure.
 **/
Image convert_from_gimp_to_image(GimpDrawable *drawable){
    gint         i, j, k, channels;
    gint         x1, y1, x2, y2;
    GimpPixelRgn rgn_in;
    guchar      *row; /* to read drawable row by row (more efficient) */
    guchar      *outrow;
    int width,height;
    Image out;

    /* allows calculation of the filter's effect limits
    excluding any region that is not in the active selection */
    gimp_drawable_mask_bounds (drawable->drawable_id,
                               &x1, &y1,
                               &x2, &y2);
    /* returns the number of bytes per pixel (or the number of channels)
      for the specified drawable */
    channels = gimp_drawable_bpp (drawable->drawable_id);
    assert(channels>=3);

    /* to read pixels data in the active selection */
    gimp_pixel_rgn_init (&rgn_in,
                         drawable,
                         x1, y1,
                         x2 - x1, y2 - y1,
                         FALSE, FALSE);

    /* initialise enough memory for row1, row2, row3, outrow */
    row = g_new (guchar, channels * (x2 - x1));

    width = x2-x1;
    height = y2-y1;

    /* store data in Image */
    out = create_Image(height,width);
    if(channels==3){
        for (i = y1; i < y2; i++) {
            /* get row i */
            gimp_pixel_rgn_get_row (&rgn_in,
                                    row,
                                    x1, i,
                                    x2 - x1);

            for (j = x1; j < x2; j++)
            {
                out->red[i-y1][j-x1] = row[channels * (j - x1) ]  ;
                out->green[i-y1][j-x1] = row[channels * (j - x1) +1]  ;
                out->blue[i-y1][j-x1] = row[channels * (j - x1) +2]  ;
            }


            if (i % 10 == 0){
                gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
            }
        }
    }else{
        /* image has an alpha channel
          each pixel wich is not completly opaque is black */
        for (i = y1; i < y2; i++) {
            /* get row i */
            gimp_pixel_rgn_get_row (&rgn_in,
                                    row,
                                    x1, i,
                                    x2 - x1);


            for (j = x1; j < x2; j++)
            {
                if(row[channels * (j - x1) +3]>10){
                    out->red[i-y1][j-x1] = row[channels * (j - x1) ]  ;
                    out->green[i-y1][j-x1] = row[channels * (j - x1) +1]  ;
                    out->blue[i-y1][j-x1] = row[channels * (j - x1) +2]  ;
                }else{

                    out->red[i-y1][j-x1] = 0  ;
                    out->green[i-y1][j-x1] = 0  ;
                    out->blue[i-y1][j-x1] = 0  ;
                }
            }

            if (i % 10 == 0){
                gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
            }
        }
    }

    g_free (row);

    gimp_drawable_flush (drawable);
    gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
    gimp_drawable_update (drawable->drawable_id,
                          x1, y1,
                          x2 - x1, y2 - y1);
    return out;
}

/**
 * @brief store Image data in Gimp's drawable data structure
 **/
void convert_from_Image_to_Gimp(Image image, GimpDrawable *drawable){
    gint         i, j, k, channels;
    gint         x1, y1, x2, y2;
    GimpPixelRgn rgn_out;
    guchar      *row; /* to read drawable row by row (more efficient) */
    guchar      *outrow;
    int width,height;

    /* allows calculation of the filter's effect limits
       excluding any region that is not in the active selection */
    gimp_drawable_mask_bounds (drawable->drawable_id,
                               &x1, &y1,
                               &x2, &y2);
    /* returns the number of bytes per pixel (or the number of channels)
      for the specified drawable */
    channels = gimp_drawable_bpp (drawable->drawable_id);
    assert(channels>=3);

    /* to read pixels data in the active selection */
    gimp_pixel_rgn_init (&rgn_out,
                         drawable,
                         x1, y1,
                         x2 - x1, y2 - y1,
                         TRUE, TRUE);

    /* initialise enough memory for row1, row2, row3, outrow */
    row = g_new (guchar, channels * (x2 - x1));

    width = x2-x1;
    height = y2-y1;

    if(channels==3){
        for (i = y1; i < y2; i++) {
            /* get row i */
            gimp_pixel_rgn_get_row (&rgn_out,
                                    row,
                                    x1, i,
                                    x2 - x1);

            for (j = x1; j < x2; j++)
            {
                row[channels * (j - x1) ] =image->red[i-y1][j-x1]  ;
                row[channels * (j - x1) +1] = image->green[i-y1][j-x1]  ;
                row[channels * (j - x1) +2] = image->blue[i-y1][j-x1] ;
            }

            gimp_pixel_rgn_set_row (&rgn_out,
                                    row,
                                    x1, i,
                                    x2 - x1);
            if (i % 10 == 0){
                gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
            }
        }
    }else{
        for (i = y1; i < y2; i++) {
            /* get row i */
            gimp_pixel_rgn_get_row (&rgn_out,
                                    row,
                                    x1, i,
                                    x2 - x1);

            for (j = x1; j < x2; j++)
            {
                row[channels * (j - x1) ] =image->red[i-y1][j-x1]  ;
                row[channels * (j - x1) +1] = image->green[i-y1][j-x1]  ;
                row[channels * (j - x1) +2] = image->blue[i-y1][j-x1] ;
                row[channels * (j - x1) +3] = 255 ;
            }

            gimp_pixel_rgn_set_row (&rgn_out,
                                    row,
                                    x1, i,
                                    x2 - x1);
            if (i % 10 == 0){
                gimp_progress_update ((gdouble) (i - y1) / (gdouble) (y2 - y1));
            }
        }
    }

    g_free (row);

    gimp_drawable_flush (drawable);
    gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
    gimp_drawable_update (drawable->drawable_id,
                          x1, y1,
                          x2 - x1, y2 - y1);
}

/**
 * @brief run scis as gimp plugin
 **/
void runSCIS(const GimpParam *param,
             gint             *nreturn_vals,
             GimpParam       **return_vals){

    gint num_layers;
    gint32 image_ID;
    gint *layers ;
    gint32 imageLayer;
    gint32 seedsLayer;
    gint32 segmentationLayer;
    GimpDrawable     *image_drawable;
    GimpDrawable     *seeds_drawable;
    GimpDrawable     *segmentation_drawable;

    Image image;
    Image seeds;
    Image segmentation;
    SVMAndSuperpixels algo;

    static GimpParam  values[1];
    GimpPDBStatusType status = GIMP_PDB_SUCCESS;


    /* Setting mandatory output values */
    *nreturn_vals = 1;
    *return_vals  = values;

    values[0].type = GIMP_PDB_STATUS;
    values[0].data.d_status = status;

    image_ID = param[1].data.d_int32;
    layers = gimp_image_get_layers(image_ID, &num_layers);
    if(num_layers < 2  ){
        const GimpRunMode run_mode = (const GimpRunMode) param[0].data.d_int32;
        STATE=3;
        errror_handler();
        return;
    }
    if(num_layers==2){
        gboolean success;
        imageLayer = layers[1];
        seedsLayer = layers[0];
        /* get image drawable */
        image_drawable = gimp_drawable_get (imageLayer);
        /* get seeds drawable */
        seeds_drawable = gimp_drawable_get(seedsLayer);
        /* create layer with segmentation */
        segmentationLayer =  gimp_layer_new(image_ID,NULL,image_drawable->width,image_drawable->height,GIMP_RGB_IMAGE,50,GIMP_NORMAL_MODE);
        /* insert layer with segmentation */
        success = gimp_image_insert_layer(image_ID,segmentationLayer,0,0);
        /* get segmentation drawable */
        segmentation_drawable = gimp_drawable_get (segmentationLayer);
    }else{
        imageLayer = layers[num_layers-1];
        seedsLayer = layers[num_layers-2];
        segmentationLayer = layers[num_layers-3];
        /* get image drawable */
        image_drawable = gimp_drawable_get (imageLayer);
        /* get seeds drawable */
        seeds_drawable = gimp_drawable_get(seedsLayer);
        /* get segmentation drawable */
        segmentation_drawable = gimp_drawable_get (segmentationLayer);
    }

    gimp_progress_init("load original image...");
    image = convert_from_gimp_to_image(image_drawable);
    gimp_progress_end();
    gimp_progress_init("load seeds...");
    seeds = convert_from_gimp_to_image(seeds_drawable);
    gimp_progress_end();



    algo.C_=4;
    algo.minSize_felzenswalb_=20;
    algo.gamma_=4;
    algo.nbFeatures=5;

    /* run SCIS */
    segmentation =  SVMAndSuperpixels_segment(&algo,seeds,image,true);
    if(segmentation){
        gimp_progress_init("drawing segmentation ...");
        convert_from_Image_to_Gimp(segmentation,segmentation_drawable);
        /* set opacity  to result at 50 */
        gimp_layer_set_opacity(segmentationLayer,50);

        gimp_progress_end();

        delete_Image(&segmentation);
    }else{
        STATE=4;
        errror_handler();
    }


    delete_Image(&image);
    delete_Image(&seeds);
    gimp_displays_flush ();
    gimp_drawable_detach (image_drawable);
    gimp_drawable_detach (seeds_drawable);
    gimp_drawable_detach (segmentation_drawable);

    return;
}

/**
 * Gimp stuff
 *
 **/

static void query (void);
static void run   (const gchar      *name,
                   gint              nparams,
                   const GimpParam  *param,
                   gint             *nreturn_vals,
                   GimpParam       **return_vals);

GimpPlugInInfo PLUG_IN_INFO =
{
    NULL,
    NULL,
    query,
    run
};


MAIN()



static void
query (void)
{
    char * pluginName="plug-in-SCIS ";
    char * pluginBlurb="Segmentation interactive";
    char * pluginHelp="Multi-class segmentation using pixels labeled by user";
    char * pluginMenuLabel="_SCIS ...";
    static GimpParamDef args[] =
    {
        {
            GIMP_PDB_INT32,
            (gchar*)"run-mode",
            (gchar*)"Run mode"
        },
        {
            GIMP_PDB_IMAGE,
            (gchar*)"image",
            (gchar*)"Input image"
        },
        {
            GIMP_PDB_DRAWABLE,
            (gchar*)"drawable",
            (gchar*)"Input drawable"
        }
    };


    gimp_install_procedure (
                pluginName,
                pluginBlurb,
                pluginHelp,
                "Berengere Mathieu",
                "CeCILL-C ",
                "2016",
                pluginMenuLabel,
                "RGB*",
                GIMP_PLUGIN,
                G_N_ELEMENTS (args), 0,
                args, NULL);

    gimp_plugin_menu_register (pluginName,
                               "<Image>/Filters/Segmentation");
}

static void
run (const gchar      *name,
     gint              nparams,
     const GimpParam  *param,
     gint             *nreturn_vals,
     GimpParam       **return_vals)
{

    STATE=0;
    const GimpRunMode run_mode = (const GimpRunMode) param[0].data.d_int32;
    interactive_run_mode = (run_mode==GIMP_RUN_INTERACTIVE);

    runSCIS(param,nreturn_vals,return_vals);
}



