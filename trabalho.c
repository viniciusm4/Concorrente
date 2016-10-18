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

sem_t prod;
static long seed = DEFAULT;
double dt, dt_old;  /* Alterado de static para global */
int contador;       /* numero de aprticulas inicializadas */
int tamanho_laco;   /* tamanho dos blocos de laço das threads*/
int npart;

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


//necessario para passar mais de um parametro na thread
typedef struct{
    Particle particles; 
    ParticleV pv; 
    int posicao; 
}ParametrosThread;

/*! para poder ser acessado no metodo initParticles
* tem que ficar localizado abaixo da struct
*/
Particle  * particles;   /* Particles */
ParticleV * pv;          /* Particle velocity */
//***************************************************


void *InitParticles( void *parametro )
{
  ParametrosThread my_parametro = *((ParametrosThread *) parametro);

  int i = my_parametro.posicao * tamanho_laco;
  int aux = i + tamanho_laco;

  for(i; i < aux; i++){
    printf("%d\n", i);

    
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
  
    /* 
    my_parametro.particles[i].x    = Random();
    my_parametro.particles[i].y    = Random();
    my_parametro.particles[i].z    = Random();
    my_parametro.particles[i].mass = 1.0;
    my_parametro.pv[i].xold    = my_parametro.particles[i].x;
    my_parametro.pv[i].yold    = my_parametro.particles[i].y;
    my_parametro.pv[i].zold    = my_parametro.particles[i].z;
    my_parametro.pv[i].fx    = 0;
    my_parametro.pv[i].fy    = 0;
    my_parametro.pv[i].fz    = 0;

    particles[i] = my_parametro.particles[i];
    pv[i] = my_parametro.pv[i];
    */
    sem_wait(&prod);
    contador++;
    if(contador == npart){
        i = aux;
    }
    sem_post(&prod);
  }
}

double ComputeForces( Particle [], Particle [], ParticleV [], int );
double ComputeNewPos( Particle [], ParticleV [], int, double);

int main(int argc, char **argv)
{

    int    i, j;
    int    cnt;          /* number of times in loop */
    int    numt;         /* numero de threads */
    double sim_t;        /* Simulation time */
    int    tmp;
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

    sim_t = 0.0;

    //parte que cria as threads p executar initParticles
    //*************************************************
    //InitParticles( particles, pv, npart);
    ParametrosThread * parametros;
    parametros = (ParametrosThread *) malloc(sizeof(ParametrosThread)*npart);

    sem_init(&prod, 0, 1);
  
    int aux;
    if(npart < numt){
      aux = npart;
    }else{
      aux = numt;
    }

    tamanho_laco = (int)(npart / (numt-1));

    pthread_t thread[aux];
    /* Generate the initial values */
    for (long int i = 0; i < aux; i++){
      parametros[i].posicao = i;
      pthread_create(&thread[i], NULL, InitParticles, (void *)&parametros[i]);
    }

    for (long int i = 0; i < aux; i++) {
      //identificador da thread
      //variavel que ira armazenar o valor retornado pela pthread_exit()
      pthread_join(thread[i], NULL);
    }

    while (cnt--) {
      double max_f;
      /* Compute forces (2D only) */
      max_f = ComputeForces( particles, particles, pv, npart );
      /* Once we have the forces, we compute the changes in position */
      sim_t += ComputeNewPos( particles, pv, npart, max_f);
    }
    for (i=0; i<npart; i++)
      fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);
    return 0;

    sem_destroy(&prod);

    //finaliza a execucao de uma thread
    pthread_exit(NULL);
}

double ComputeForces( Particle myparticles[], Particle others[], ParticleV pv[], int npart )
{
  double max_f;
  int i;
  max_f = 0.0;
  for (i=0; i<npart; i++) {
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
    if (fx > max_f) max_f = fx;
  }
  return max_f;
}

double ComputeNewPos( Particle particles[], ParticleV pv[], int npart, double max_f)
{
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
  return dt_old;
}
