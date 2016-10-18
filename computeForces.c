#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

sem_t maxF, cheio, newPos;
static long seed = DEFAULT;
double dt, dt_old;        /* Alterado de static para global */
double max_f;
double sim_t;             /* Simulation time */
int npart, contador;  
int cnt;                  /* number of times in loop */
int tamanho_laco;         /* tamanho dos blocos de laço das threads*/
int numt;                 /* numero de threads */  

pthread_mutex_t barrier;  
pthread_cond_t go;

double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0) 
    seed = t;
  else 
    seed = t + MODULUS;
  return ((double) seed / MODULUS);
}

/*
 * End of the pRNG algorithm
 */

typedef struct {
    double x, y, z;
    double mass;
} Particle;

typedef struct {
    double xold, yold, zold;
    double fx, fy, fz;
} ParticleV;

//necessario para passar parametro na thread
typedef struct{ 
    int posicao; 
}ParametrosThread;

/*! para poder ser acessado no metodo computeForces
* tem que ficar localizado abaixo da struct
*/
Particle  * particles;   /* Particles */
ParticleV * pv;          /* Particle velocity */
//***************************************************

void *ComputeForces(void *parametro)
{
  ParametrosThread my_parametro = *((ParametrosThread *) parametro);
  int posicao = my_parametro.posicao;
  int loop = cnt;
  Particle * others = particles;
  Particle * myparticles = particles;

  int i = posicao * tamanho_laco;
  int backup = i;
  int aux = i + tamanho_laco;
  if((posicao+1) == numt-1){
    aux += npart%numt;
  }

  while(loop--){

    i = backup;
    //printf("%d ..... %d\n", i, aux);

    for(i; i < aux; i++){
      int j;
      double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
      rmin = 100.0;
      xi   = myparticles[i].x;
      yi   = myparticles[i].y;
      fx   = 0.0;
      fy   = 0.0;
      for (j=0; j<npart; j++) {
        rx = xi - others[j].x;
        ry = yi - others[j].y;
        mj = others[j].mass;
        r  = rx * rx + ry * ry;
        /* ignore overlap and same particle */
        if (r == 0.0) continue;
        if (r < rmin) rmin = r;
        r  = r * sqrt(r);
        fx -= mj * rx / r;
        fy -= mj * ry / r;
      }
      pv[i].fx += fx;
      pv[i].fy += fy;
      fx = sqrt(fx*fx + fy*fy)/rmin;
    
      sem_wait(&maxF);
      if (fx > max_f) max_f = fx;
      sem_post(&maxF);
    }

    pthread_mutex_lock(&barrier);
    contador++;
    if(contador < numt-1){
      pthread_cond_wait(&go, &barrier);
    }else{
      contador = 0;
      sem_post(&newPos);
      sem_wait(&cheio);
      pthread_cond_broadcast(&go);
    }
    pthread_mutex_unlock(&barrier);

  }
    pthread_exit(NULL);
}

void *ComputeNewPos()
{
  int loop = cnt;
  while(loop--){

    sem_wait(&newPos);
    //printf("NEW\n");

    int i;
    double a0, a1, a2;
    double dt_new;
    a0   = 2.0 / (dt * (dt + dt_old));
    a2   = 2.0 / (dt_old * (dt + dt_old));
    a1   = -(a0 + a2);
    for (i=0; i<npart; i++) {
     double xi, yi;
     xi             = particles[i].x;
     yi             = particles[i].y;
     particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
     particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
     pv[i].xold     = xi;
     pv[i].yold     = yi;
     pv[i].fx       = 0;
     pv[i].fy       = 0;
    }
    dt_new = 1.0/sqrt(max_f);
    /* Set a minimum: */
    if (dt_new < 1.0e-6) dt_new = 1.0e-6;
    /* Modify time step */
    if (dt_new < dt) {
     dt_old = dt;
     dt     = dt_new;
    }
    else if (dt_new > 4.0 * dt) {
     dt_old = dt;
     dt    *= 2.0;
    }

    sim_t += dt_old;

    sem_post(&cheio);
  }
  pthread_exit(NULL);
}

void InitParticles( Particle[], ParticleV [], int );

int main(int argc, char **argv)
{
    int         i, j;
    int tmp;
    if(argc != 4){
    printf("Wrong number of parameters.\nUsage: nbody num_bodies timesteps\n");
    exit(1);
  }
  npart = atoi(argv[1]);
  cnt = atoi(argv[2]);
  numt = atoi(argv[3]);
  dt = 0.001; 
  dt_old = 0.001;

    /* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
    
    /* Generate the initial values */
    InitParticles( particles, pv, npart);
    sim_t = 0.0;

    /* Inicializando semaforos*/
    sem_init(&maxF, 0, 1);
    sem_init(&newPos, 0, 0);
    sem_init(&cheio, 0, 0);

    int aux;
    if(npart < numt){
      aux = npart;
      tamanho_laco = 1;
    }else{
      aux = numt;
      tamanho_laco = (int)(npart / (numt-1));
      //printf("tamanho laco: %d\n", tamanho_laco);
    }

    ParametrosThread * parametros;
    parametros = (ParametrosThread *) malloc(sizeof(ParametrosThread)*aux);

    pthread_t thread[aux];
    
    for (long int i = 0; i < aux; i++){
      if(i < (aux-1)){
        parametros[i].posicao = i;
        /* Compute forces (2D only) */
        pthread_create(&thread[i], NULL, ComputeForces, (void*)&parametros[i]);
      }else{
         /* Once we have the forces, we compute the changes in position */
        pthread_create(&thread[i], NULL, ComputeNewPos, NULL);
      }
    }

    for (long int i = 0; i < aux+1; i++) {
      //identificador da thread
      //variavel que ira armazenar o valor retornado pela pthread_exit()
      pthread_join(thread[i], NULL);
    }
    
    //for (i=0; i<npart; i++)
    //  fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);
    
    sem_destroy(&maxF);
    sem_destroy(&newPos);
    sem_destroy(&cheio);

    return 0;
}

void InitParticles( Particle particles[], ParticleV pv[], int npart )
{
    int i;
    for (i=0; i<npart; i++) {
  particles[i].x    = Random();
  particles[i].y    = Random();
  particles[i].z    = Random();
  particles[i].mass = 1.0;
  pv[i].xold    = particles[i].x;
  pv[i].yold    = particles[i].y;
  pv[i].zold    = particles[i].z;
  pv[i].fx    = 0;
  pv[i].fy    = 0;
  pv[i].fz    = 0;
    }
}
