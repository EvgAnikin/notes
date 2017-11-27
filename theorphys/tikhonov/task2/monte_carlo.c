#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265359


double partial_action(double x, double xprev, double xnext, double alpha, double tau)
{
    return (x - xprev)*(x - xprev)/2./tau + 
           (x - xnext)*(x - xnext)/2./tau +
           (x*x/2 + alpha*x*x*x*x)*tau;
}


int metropolis_step(double *X, int N, double tau, double alpha)
{
    int i;
    i = rand() % N;

    double x_max, x_new, old_action, new_action;
    x_max = 10;
    x_new = -x_max + 2*x_max*drand48();
    
    old_action = partial_action(X[i],  X[(i-1+N)%N], X[(i+1)%N], alpha, tau);
    new_action = partial_action(x_new, X[(i-1+N)%N], X[(i+1)%N], alpha, tau);

    float p = drand48();
    if (p < exp(old_action - new_action))
    {
        X[i] = x_new;
        return 1;
    }
    return 0;
}


float myErfInv(float x)
{
    float tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0f : 1.0f;
    
    x = (1 - x)*(1 + x);        // x = 1 - x*x;
    lnx = logf(x);
    
    tt1 = 2/(PI*0.147) + 0.5f * lnx;
    tt2 = 1/(0.147) * lnx;
    
    return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}


int heatbath_step(double *X, int N, double tau)
{
    int i;
    double sigma;
    i = rand() % N;
    
    sigma = sqrt(tau/(2 + tau*tau));

    double gaussian_random = myErfInv(2*drand48() - 1);
    if (gaussian_random == NAN || 
        gaussian_random == -INFINITY ||
        gaussian_random == INFINITY)
    {
        gaussian_random = 0;
    }

    X[i] = (X[(i-1+N)%N] + X[(i+1)%N])/(2 + tau*tau) + 
           sqrt(2)*sigma*gaussian_random;

    return 1;
}


double *init_X(int N)
{
    double *X;
    X = (double *)malloc(sizeof(double)*N);

    int i;
    for(i = 0; i < N; i++)
    {
        X[i] = -0.5 + 2*drand48();
    } 
    return X;
}

double mean(double *X, int N)
{
    int i;
    double sum = 0;
    for(i = 0; i < N; i++)
    {
        sum += X[i];
    }
    return sum/N;
}

double mean_energy(double *X, int N, double alpha)
{
    int i;
    double sum_sq = 0;
    for(i = 0; i < N; i++)
    {
        sum_sq += X[i]*X[i] + 3*alpha*X[i]*X[i]*X[i]*X[i];
    }
    return sum_sq/N;
}


double error(double *X, int N)
{
    double mean_value = mean(X, N);
    double sum_sq_dev = 0;
    int i;

    for(i = 0; i < N; i++)
    {
        sum_sq_dev += (X[i] - mean_value)*(X[i] - mean_value);
    }

    return sqrt(sum_sq_dev/(N - 1)/N);
}


void set_args(char **argv, int *nsamples, int *nbetween, int *N, double *beta, double *alpha)
{
    sscanf(argv[1], "%d", nsamples);
    sscanf(argv[2], "%d", nbetween);
    sscanf(argv[3], "%d", N);
    sscanf(argv[4], "%lf", beta);
    sscanf(argv[5], "%lf", alpha);
}


/*double mean_nl_oscillator_energy(int nsamples, int nbetween, int N, double beta, double alpha)
{
    mean_energies = (double *)malloc(nsamples*sizeof(double));
    
    tau = beta/N;
    X = init_X(N);
    ninitial = 5*nbetween;
    for(i = 0; i < ninitial; i++)
    {
        nsteps += metropolis_step(X, N, tau, alpha);
    }

    for(i = 0; i < nsamples; i++)
    {
        for(j = 0; j < nbetween; j++)
        {
            nsteps += metropolis_step(X, N, tau, alpha);
        }
        mean_energies[i] = mean_energy(X, N, alpha);
    }
}*/


int main(int argc, char *argv[])
{
    srand(time(NULL));
    
    int i, j, N, ninitial, nsamples, nbetween;
    int nsteps = 0;
    double *X;
    double alpha, beta, tau;
    double *mean_energies;

    if (argc != 6)
    {   
        printf("Arguments missing\n");
        exit(1);
    }

    set_args(argv, &nsamples, &nbetween, &N, &beta, &alpha);
    mean_energies = (double *)malloc(nsamples*sizeof(double));
    
    tau = beta/N;
    X = init_X(N);
    ninitial = 5*nbetween;
    for(i = 0; i < ninitial; i++)
    {
        nsteps += metropolis_step(X, N, tau, alpha);
    }

    for(i = 0; i < nsamples; i++)
    {
        for(j = 0; j < nbetween; j++)
        {
            nsteps += metropolis_step(X, N, tau, alpha);
        }
        mean_energies[i] = mean_energy(X, N, alpha);
    }


    printf("fraction of steps: %f\n", ((float) nsteps)/nsamples/nbetween);
    printf("average of steps between measurements: %f\n", ((float) nsteps)/nsamples);
    
    printf("mean square: %lf\nerror: %lf\n", mean(mean_energies, nsamples), 
                                             error(mean_energies, nsamples));

    free(X);
    free(mean_energies);
    return 0;
}
