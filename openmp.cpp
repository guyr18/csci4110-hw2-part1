#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "grid.cpp"
#include "omp.h"
#include <vector>
#include <iostream>

bool inBounds(const int i, const int j)
{

    return (i >= 0 && i < numBins && j >= 0 && j < numBins);
                        
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0, numthreads;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);
    
    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen (sumname, "a") : NULL;

    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);
    
    // Simulate a number of time steps.
    double simulation_time = read_timer();
	
    omp_set_num_threads(16);

    #pragma omp parallel private(dmin)
    {
	numthreads = omp_get_num_threads();

    for(int step = 0; step < NSTEPS; step++)
    {

	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;
    
        double density = .0005;
        double size = sqrt(density * n);
        double binSize = size / numBins;

        // Create a matrix grid to hold bins.
        grid_t grid;

        // Initialize grid.
        grid.init();
        
        #pragma openmp for
        // Fill matrix bins with particles.
        for(int i = 0; i < n; i++)
        {

            int binPosX = floor(particles[i].x / binSize);
            int binPosY = floor(particles[i].y / binSize);
            grid.collection[binPosY][binPosX].push_back(&particles[i]);

        }


        #pragma omp for reduction (+:navg) reduction (+:davg)
        // Calculate forces for local neighborhood(s).
        for (int curRow = 0; curRow < numBins; curRow++)
        {
            for (int curCol = 0; curCol < numBins; curCol++)
            {
            
                std::vector<bin_t> neighbors; // Cache neighbors.

                // Fill neighbors.
                for (int i = curRow - 1; i <= curRow + 1; i++) 
                {

                    for (int j = curCol - 1; j <= curCol + 1; j++)
                    {

                        // Check if in bounds.
                        if (inBounds(i, j))
                        {
    
                            neighbors.push_back(grid.collection[i][j]);

                        }
                    }
                }

                bin_t curBin = grid.collection[curRow][curCol];

                for (int x = 0; x < curBin.size(); x++) 
                {

                    // Stop current bin particle from accelerating.
                    curBin[x]->ax = 0;
                    curBin[x]->ay = 0;


                    // Loop across neighbor bins.
                    for (int y = 0; y < neighbors.size(); y++)
                    {

                        bin_t curNeighbor = neighbors[y];

                        // Loop across particles that are in range.
                        for (int z = 0; z < curNeighbor.size(); z++)
                        {

                            // Apply force.
                            apply_force(*curBin[x], *curNeighbor[z], &dmin, &davg, &navg);

                        }
                    }
                }
            }
        }

        #pragma omp for
        // Move particles.
        for( int i = 0; i < n; i++ )  
        {

            move(particles[i]);

        }

        if(find_option(argc, argv, "-no") == -1)
        {
         
	  #pragma omp master	
          // Computing statistical data.
          if (navg) 
          {

            absavg +=  davg / navg;
            nabsavg++;

          }
 
          if (dmin < absmin) absmin = dmin;
		
          #pragma omp master
          // Save if neccessary.
          if(fsave && (step % SAVEFREQ) == 0)
              save(fsave, n, particles);

        }
    }
 }

    simulation_time = read_timer() - simulation_time;
    
    printf("n = %d, threads = %d, simulation time = %g seconds", n, numthreads, simulation_time);

    if(find_option(argc, argv, "-no") == -1)
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf(", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    
    // Printing summary data.
    if(fsum) 
        fprintf(fsum, "%d %d %g\n", n, numthreads, simulation_time);
 
    // Clearing space.
    if(fsum)
        fclose(fsum);    
    free(particles);
    if(fsave)
        fclose(fsave);
    
    return 0;
}