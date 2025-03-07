// This version uses pthreads to parallelize the simulation loop.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

// Struct to store the position and velocity of a body
typedef struct
{
    double x;
    double y;
} Vector2D;

// Creating a global ACCx and ACCy array to store the acceleration of all particles for each thread
double *GACCx;
double *GACCy;

// Function prototypes for the serial functions (not used by the parallel code)
// We now inline the work in the threads.
static double get_wall_seconds();

// Structure to pass thread arguments
typedef struct
{
    int tid;            // thread id
    int nthreads;       // number of threads
    int nsteps;         // number of timesteps
    int nstars;         // total number of stars
    double stepsize;    // size of each timestep
    double G;           // gravitational constant (scaled)
    double e0;          // softening factor
    Vector2D *position; // position of each star containing x and y
    Vector2D *velocity; // velocity of each star containing x and y
    double *mass;
    double *ACCx;
    double *ACCy;
} ThreadData;

// Global barrier for synchronization
pthread_barrier_t barrier;

void *simulation_thread(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    int nstars = data->nstars;
    double stepsize = data->stepsize;
    double G = data->G;
    double e0 = data->e0;
    int tid = data->tid;
    int nthr = data->nthreads;

    // --- Force Calculation ---
    for (int t = 0; t < data->nsteps; t++)
    {
        // Reset Global acceleration arrays to zero
        for (int i = 0; i < nstars; i++)
        {
            GACCx[i + (tid * nstars)] = 0;
            GACCy[i + (tid * nstars)] = 0;
        }

        // Synchronize to ensure all threads have reset acceleration arrays
        pthread_barrier_wait(&barrier);
        // Wait for all threads to reset acceleration arrays
        // Compute forces between stars
        for (int i = tid; i < nstars; i += nthr)
        {
            double posx = data->position[i].x;
            double posy = data->position[i].y;
            double mass_i = data->mass[i];
            // This time we just accumulate the forces because the force
            // on each star is updated by mutliple threads
            // This is because of Newton's 3rd law
            double accx = 0;
            double accy = 0;
            for (int j = i + 1; j < nstars; j++)
            {

                double dx = posx - data->position[j].x;
                double dy = posy - data->position[j].y;
                double rije0 = sqrt(dx * dx + dy * dy) + e0;
                double pow_rije0 = 1.0 / (rije0 * rije0 * rije0);

                double factor_i = data->mass[j] * pow_rije0;
                double factor_j = mass_i * pow_rije0;

                accx -= factor_i * dx;
                accy -= factor_i * dy;

                GACCx[j + (tid * nstars)] += factor_j * dx;
                GACCy[j + (tid * nstars)] += factor_j * dy;
            }

            GACCx[i + (tid * nstars)] += accx;
            GACCy[i + (tid * nstars)] += accy;
        }
        // Synchronize to ensure all threads have computed forces
        pthread_barrier_wait(&barrier);
        // Now multiply by G after all threads have computed forces
        if (tid == 0)
        {
            for (int i = 0; i < nstars; i++)
            {
                data->ACCx[i] = 0;
                data->ACCy[i] = 0;
                for (int j = 0; j < nthr; j++)
                {
                    data->ACCx[i] += GACCx[i + (j * nstars)];
                    data->ACCy[i] += GACCy[i + (j * nstars)];
                }
                data->ACCx[i] *= G;
                data->ACCy[i] *= G;
            }
        }

        // Synchronize to ensure all threads have computed forces

        pthread_barrier_wait(&barrier);

        // --- Update Velocities and Positions ---
        // Update velocities and positions of stars
        if (tid == 0)
        {
            for (int i = 0; i < nstars; i++)
            {
                data->velocity[i].x += stepsize * data->ACCx[i];
                data->velocity[i].y += stepsize * data->ACCy[i];
                data->position[i].x += stepsize * data->velocity[i].x;
                data->position[i].y += stepsize * data->velocity[i].y;
            }
        }
        // Synchronize to ensure all threads have finished updating before next timestep
        pthread_barrier_wait(&barrier);
        // This is the end of the timestep
    }
    return NULL;
}

int main(int argc, char *argv[])
{
    double start_time = get_wall_seconds();

    if (argc != 7)
    {
        fprintf(stderr, "Usage: %s <Number of stars> <input file> <Number of timesteps> <size of timesteps> <graphics> <Nthreads>\n", argv[0]);
        return 1;
    }

    // Read command line arguments (make const for potential optimizations)
    const int nstars = atoi(argv[1]);
    const char *input_file = argv[2];
    const int nsteps = atoi(argv[3]);
    const double stepsize = atof(argv[4]);
    const int graphics = atoi(argv[5]); // 0 or 1 for false or true
    int NUM_THREADS = atoi(argv[6]);
    if (graphics)
    {
        fprintf(stderr, "Graphics not supported\n");
        return 1;
    }

    // Limit number of threads to number of stars
    if (NUM_THREADS > nstars)
    {
        NUM_THREADS = nstars;
    }

    const double G = 100.0 / nstars;
    const double e0 = 1e-3; // softening factor

    // Allocate arrays for simulation
    // We pad the arrays to avoid false sharing
    int padding_size = nstars;
    Vector2D *position = (Vector2D *)malloc(padding_size * sizeof(Vector2D));
    double *mass = (double *)malloc(nstars * sizeof(double));
    Vector2D *velocity = (Vector2D *)malloc(padding_size * sizeof(Vector2D));
    double *brightness = (double *)malloc(nstars * sizeof(double));
    double *ACCx = (double *)malloc(padding_size * sizeof(double));
    double *ACCy = (double *)malloc(padding_size * sizeof(double));
    GACCx = (double *)malloc(padding_size * NUM_THREADS * sizeof(double));
    GACCy = (double *)malloc(padding_size * NUM_THREADS * sizeof(double));

    if (position == NULL || mass == NULL || velocity == NULL || brightness == NULL || ACCx == NULL || ACCy == NULL)
    {
        fprintf(stderr, "Error allocating memory\n");
        return 1;
    }

    // Initialize arrays to zero
    memset(position, 0, nstars * sizeof(Vector2D));
    memset(mass, 0, nstars * sizeof(double));
    memset(velocity, 0, nstars * sizeof(Vector2D));
    memset(brightness, 0, nstars * sizeof(double));
    memset(ACCx, 0, nstars * sizeof(double));
    memset(ACCy, 0, nstars * sizeof(double));

    // Read input binary file
    FILE *file = fopen(input_file, "rb");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", input_file);
        free(position);
        free(mass);
        free(velocity);
        free(brightness);
        free(ACCx);
        free(ACCy);
        return 1;
    }
    for (int i = 0; i < nstars; i++)
    {
        if (fread(&position[i], sizeof(Vector2D), 1, file) != 1 ||
            fread(&mass[i], sizeof(double), 1, file) != 1 ||
            fread(&velocity[i], sizeof(Vector2D), 1, file) != 1 ||
            fread(&brightness[i], sizeof(double), 1, file) != 1)
        {
            fprintf(stderr, "Error reading file\n");
            fclose(file);
            return 1;
        }
    }
    fclose(file);

    // Determine number of threads to use and distribute work
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];

    // Initialize barrier for NUM_THREADS threads
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    // Create worker threads
    for (int t = 0; t < NUM_THREADS; t++)
    {
        thread_data[t].tid = t;
        thread_data[t].nthreads = NUM_THREADS;
        // Distribute remainder stars to the last thread
        thread_data[t].nstars = nstars;
        thread_data[t].nsteps = nsteps;
        thread_data[t].stepsize = stepsize;
        thread_data[t].G = G;
        thread_data[t].e0 = e0;
        thread_data[t].position = position;
        thread_data[t].velocity = velocity;
        thread_data[t].mass = mass;
        thread_data[t].ACCx = ACCx;
        thread_data[t].ACCy = ACCy;

        if (pthread_create(&threads[t], NULL, simulation_thread, (void *)&thread_data[t]) != 0)
        {
            fprintf(stderr, "Error creating thread %d\n", t);
            exit(1);
        }
    }

    // Wait for all threads to finish
    for (int t = 0; t < NUM_THREADS; t++)
    {
        pthread_join(threads[t], NULL);
    }

    // Destroy barrier and mutexes
    pthread_barrier_destroy(&barrier);

    // Write out the results to a binary file
    FILE *output = fopen("result.gal", "wb");
    if (output == NULL)
    {
        fprintf(stderr, "Error opening output file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++)
    {
        fwrite(&position[i], sizeof(Vector2D), 1, output);
        fwrite(&mass[i], sizeof(double), 1, output);
        fwrite(&velocity[i], sizeof(Vector2D), 1, output);
        fwrite(&brightness[i], sizeof(double), 1, output);
    }
    fclose(output);

    // Free allocated memory
    free(ACCx);
    free(ACCy);
    free(position);
    free(mass);
    free(velocity);
    free(brightness);

    double end_time = get_wall_seconds();
    printf("Time taken: %f seconds\n", end_time - start_time);
    return 0;
}

static double get_wall_seconds()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec / 1000000;
}