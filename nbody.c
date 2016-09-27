#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;
double dt, dt_old; /* Alterado de static para global */

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

void InitParticles( Particle [], ParticleV [], int );
double ComputeForces( Particle [], Particle [], ParticleV [], int );
double ComputeNewPos( Particle [], ParticleV [], int, double);
void SegundoFor(double xi, double yi, double mj, double others_x, double others_y, double &vetorFx[], double &vetorFy[], unsigned int &k);

int main(int argc, char **argv)
{
    Particle  * particles;   /* Particles */ // Partículas com somente a posição no espaco
    ParticleV * pv;          /* Particle velocity */ // Partículas com a posição no espaco e a força atuando sobre elas (velocidade)
    int         npart, i, j;
    int         cnt;         /* number of times in loop */
    double      sim_t;       /* Simulation time */
    int tmp;
    if(argc != 3){
		printf("Wrong number of parameters.\nUsage: nbody num_bodies timesteps\n");
		exit(1);
	}
    
    // npart recebe o primeiro parâmetro, contendo o número de partículas
    // cnt recebe o segundo parâmetro, contendo o número de iterações
	npart = atoi(argv[1]);
	cnt = atoi(argv[2]);

	dt = 0.001; 
	dt_old = 0.001;

    /* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
    
    /* Generate the initial values */
    InitParticles( particles, pv, npart);
    sim_t = 0.0;

    // Laço decrementa o timesteps e chama os métodos ComputeForces e ComputeNewPos, calculando o deslocamento das partículasd
    while (cnt--) {
      double max_f;
      /* Compute forces (2D only) */
      max_f = ComputeForces( particles, particles, pv, npart );
      /* Once we have the forces, we compute the changes in position */
      sim_t += ComputeNewPos( particles, pv, npart, max_f);
    }
    //for (i = 0; i < npart; i++)
      //fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);
    return 0;
}


// Inicializa as partículas

// Cada Partícula tem 3 atributos: x, y, z (posição da partícula no espaço).
// O For seta um número randômico para cada atributo da Partícula
// E seta a massa da Partícula como 1.

// Cada PartículaV tem 6 atributos: x_old, y_old, z_old (posição antiga da partícula), fx, fy, fz (forças da partícula).
// O For seta x_old, y_old e z_old com os atributos x, y e z da Partícula, respectivamente.
// E seta fx, fy e fz como 0.

void InitParticles( Particle particles[], ParticleV pv[], int npart )
{
    int i;
    for (i = 0; i < npart; i++) {
		particles[i].x	  = Random();  // Setar os atributos de cada partícula em 1 thread, e cada atributo de cada partícula em 1 thread
		particles[i].y	  = Random();
		particles[i].z	  = Random();
		particles[i].mass = 1.0;
		pv[i].xold	  = particles[i].x;
		pv[i].yold	  = particles[i].y;
		pv[i].zold	  = particles[i].z;
		pv[i].fx	  = 0;
		pv[i].fy	  = 0;
		pv[i].fz	  = 0;
    }
}

void SegundoFor(double xi, double yi, double mj, double others_x, double others_y, double &vetorFx[], double &vetorFy[], unsigned int &k) 
{
	double rx, ry, r, rmin;

	rmin = 100.0;
	rx = xi - others_x;
	ry = yi - others_y;
	r = rx * rx + ry * ry;

	if (r == 0.0)
		continue;
	if (r < rmin)
		rmin = r;

	r = r * sqrt(r);

	pthread_mutex_lock(&mutex);
	++k;
	pthread_mutex_unlock(&mutex);

	vetorFx[k] = mj * rx / r;
	vetorFy[k] = mj * ry / r;
	sem_wait(&semaphoreEmpty);

}

void ConsumidorSegundoFor(unsigned int &p, double &fx, double &vetorFx[], double &vetorFy[]) 
{
	fx -= vetorFx[p];
	fy -= vetorFy[p];
	++p;
	sem_post(&semaphoreEmpty);
}

// ComputeForces recebe dois arrays de partículas: "myparticles[]" e "others[]", e um array de ParticleV "pv[]",
// e um inteiro "npart" (número de partículas).
// Calcula o deslocamento das partículas

double ComputeForces( Particle myparticles[], Particle others[], ParticleV pv[], int npart )
{
  double max_f;
  int i;
  max_f = 0.0;

  // Iteração de 0 até o número de partículas
  for (i = 0; i < npart; i++) {
  	unsigned int *k = 0;
  	unsigned int *p = 0;

  	sem_t semaphoreEmpty;
  	pthreat_mutex_t mutex;
  	pthread_t vetorThreads[npart];

  	sem_init(&semaphoreEmpty, 0, npart);
  	pthread_mutex_init(&mutex, NULL);

  	double* vetorFx[npart];
  	double* vetorFy[npart];

    int j;
    double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin; 
    // xi, yi : posição das particulas I | mi : massa das partículas I
    // rx, ry : diferença da posição das partículas R | mj : massa das partículas Others
    // fx, fy : força atuando sobre as partículas | r, rmin (???)
    rmin = 100.0;

    //xi   = myparticles[i].x;
    //yi   = myparticles[i].y;

    fx   = 0.0;
    fy   = 0.0;


    // -----------INÍCIO DO CÁLCULO THREADS SECUNDÁRIAS---------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------
    for (j = 0; j < npart; j++) {
    	vetorThreads[j] = pthread_create(&vetorThreads[j], NULL, SegundoFor, myparticles[i].x, myparticles[i].y, others[j].mass, others[j].x, others[j].y, *vetorFx[], *vetorFy[], *k);
      // -----------------------------------------------------
      // // rx recebe a diferença da posição X entre as duas partículas
      // // ry recebe a diferença da posição Y entre as duas partículas 
      // rx = xi - others[j].x;
      // ry = yi - others[j].y;

      // // mj recebe a massa da partícula Others
      // mj = others[j].mass;

      // // r recebe a soma dos quadrados da diferença de posição das partículas
      // r  = rx * rx + ry * ry;

      // /* ignore overlap and same particle */
      // if (r == 0.0) continue;

      // // Se r for menor que o rmin, ele recebe o valor de rmin
      // if (r < rmin) rmin = r;

      // r  = r * sqrt(r);
      // -------------------------------------------------------------------------------------------------------------
      // fx -= mj * rx / r; // FAZER (CONSUMIDOR)
      // fy -= mj * ry / r;
    }

    for (int m = 0; m < npart; ++m) {
    	pthread_join(vetorThreads[m], NULL);
    }
    // --------------------------------------------------------------------------------------------------------------------------
    // ----------FIM DO CÁLCULO THREADS SECUNDÁRIAS----------------------------------------------------------------------------------------------------------------

    // Soma fx, fy atual com as novas fx, fy calculadas
    pv[i].fx += fx;
    pv[i].fy += fy;

    // Fx recebe a força gravitacional entre duas partículas
    fx = sqrt(fx*fx + fy*fy)/rmin;

    // Estabelece max_f como a maior força gravitacional entre duas partículas calculada no método
    if (fx > max_f) 
    	max_f = fx;

  }

  // Retorna a maior força gravitacional entre duas partículas
  return max_f;
}

// Recebe um array de partículas "particles[]", um array de partículasV "pv[]"
// um inteiro "npart" com o número de partículas e um double "max_f" com a maior força gravitacional entre duas partículas
double ComputeNewPos( Particle particles[], ParticleV pv[], int npart, double max_f)
{
  int i;
  double a0, a1, a2; // Aceleração nas direções x, y e z da partícula (ax = Fx/m, ay = Fy/m)
  double dt_new;

  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);

  for (i = 0; i < npart; i++) {
    double xi, yi;
    xi	           = particles[i].x;
    yi	           = particles[i].y;
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

