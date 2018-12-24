
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

#define MASTER 0

#define LOW 0
#define HIGH 1

#define ABOVE 0
#define BELOW 1


void stencilwhole(const short nx, const short ny, float * restrict image, float * restrict  tmp_image);
void stenciltop(const short nx, const short bottom, float * restrict image, float * restrict tmp_image);
void stencilbottom(const short nx, const short top, const short bottom, float * restrict image, float * restrict tmp_image);
void stencilmiddle(const short nx, const short top, const short bottom, float * restrict image, float * restrict tmp_image);
void init_image(const short nx, const short ny, float * restrict  image, float * restrict  tmp_image);
void output_image(const char * file_name, const short nx, const short ny, float *image);
double wtime(void);
short calculateRowBoundary(const short ny, const short rank, const short noOfPartitions, short topOrBottom);
//short calculateAdjacentRank(const short direction, const short rank, const short noOfProcesses);

int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  short nx = atoi(argv[1]);
  short ny = atoi(argv[2]);
  short niters = atoi(argv[3]);

  // Allocate the image
  float * restrict image = malloc(sizeof(float)*nx*ny);
  float * restrict tmp_image = malloc(sizeof(float)*nx*ny);


  int rank;
  int size;
  int flag;
  int strlen;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;

  //initialise mpi
  MPI_Init( &argc, &argv );

  MPI_Initialized(&flag);
  if(flag!=1){
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  MPI_Get_processor_name(hostname,&strlen);
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  /* determine the RANK of the current process [0:SIZE-1] */
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  // Set the input image
  //if(rank == MASTER){
  init_image(nx, ny, image, tmp_image);
  output_image("startimage.pgm", nx, ny, image);
  //}

  // Call the stencil kernel
  double tic;
  double toc;
  // if(rank == MASTER){
  //   for (short t = 0; t < niters; ++t) {
  //     stencil(nx, ny, image, tmp_image);
  //     stencil(nx, ny, tmp_image, image);
  //   }
  // }
  // double toc = wtime();


  short first = calculateRowBoundary(ny, rank, size, LOW);
  short last = calculateRowBoundary(ny, rank, size, HIGH);

  if(size == 1){
    tic = wtime();
    for (short t = 0; t < niters; ++t) {

      stencilwhole(nx, ny, image, tmp_image);
      stencilwhole(nx, ny, tmp_image, image);
    }
    printf("\n\nsize only 1\n\n\n");
    toc = wtime();

    printf("------------------------------------\n");
    printf(" runtime size 1: %lf s\n", toc-tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, image);

  }
  else{
        if(rank == MASTER){
      tic = wtime();
      for(short t = 0; t < niters; t++){
        MPI_Send(&image[last*nx], nx, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);

        MPI_Recv(&image[(last+1)*nx], nx, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);


        stenciltop(nx, last, image, tmp_image);


        MPI_Send(&tmp_image[last*nx], nx, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);

        MPI_Recv(&tmp_image[(last+1)*nx], nx, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);
        stenciltop(nx, last, tmp_image, image);

      }
      for(short processRank = 1; processRank < size; processRank++){
        short processFirst = calculateRowBoundary(ny, processRank, size, LOW);
        short processLast = calculateRowBoundary(ny, processRank, size, HIGH);

        MPI_Recv(&image[processFirst * nx], nx*(processLast - processFirst + 1), MPI_FLOAT, processRank, 1, MPI_COMM_WORLD, &status);
      }
      toc = wtime();

      printf("------------------------------------\n");
      printf(" runtime size %d: %lf s\n", size, toc-tic);
      printf("------------------------------------\n");

      output_image(OUTPUT_FILE, nx, ny, image);


    }

    else {
      if(rank == size - 1){
        for(short t = 0; t < niters; t++){

          MPI_Recv(&image[(first-1)*nx], nx, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);

          MPI_Send(&image[first*nx], nx, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);


          stencilbottom(nx, first, last, image, tmp_image);

          MPI_Recv(&tmp_image[(first-1)*nx], nx, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);

          MPI_Send(&tmp_image[first*nx], nx, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);
          stencilbottom(nx, first, last, tmp_image, image);

        }
      }
      else{
        for(short t = 0; t < niters; t++){
          MPI_Sendrecv(&image[last*nx], nx, MPI_FLOAT, rank+1, 1, &image[(first-1)*nx], nx, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);

          MPI_Sendrecv(&image[first*nx], nx, MPI_FLOAT,rank-1, 1, &image[(last+1)*nx], nx, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);

          stencilmiddle(nx, first, last, image, tmp_image);

          MPI_Sendrecv(&tmp_image[last*nx], nx, MPI_FLOAT, rank+1, 1, &tmp_image[(first-1)*nx], nx, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);

          MPI_Sendrecv(&tmp_image[first*nx], nx, MPI_FLOAT,rank-1, 1, &tmp_image[(last+1)*nx], nx, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);

          stencilmiddle(nx, first, last, tmp_image, image);
          //send
        }
      }
      MPI_Send(&image[first*nx], nx*(last - first + 1), MPI_FLOAT, MASTER, 1, MPI_COMM_WORLD);

    }

  }



  free(image);

  //printf("Hello, world; from host %s: process %d of %d\n", hostname, rank, size);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
//used when there is only one node
void stencilwhole(const short nx, const short ny, float * restrict image, float * restrict tmp_image) {
  //when i=0
  //when j=0
  tmp_image[0] = image[0] * 0.6f;
  tmp_image[0] += image[nx] * 0.1f;
  tmp_image[0] += image[0] * 0.1f;
  //#pragma vector always
  for (int j = 1; j < ny-1; ++j) {
    tmp_image[j] = image[j] * 0.6f;
    tmp_image[j] += image[j  +nx] * 0.1f;
    tmp_image[j] += image[j-1] * 0.1f;
    tmp_image[j] += image[j+1] * 0.1f;
  }
  //when j=nx-1
  tmp_image[ny-1] = image[ny-1] * 0.6f;
  tmp_image[ny-1] += image[(ny-1) +nx] * 0.1f;
  tmp_image[ny-1] += image[(ny-1)-1] * 0.1f;

  //#pragma vector always
  for (int i = 1; i < nx-1; ++i) {
    //when j=0
    tmp_image[i*nx] = image[i*nx] * 0.6f;
    tmp_image[i*nx] += image[(i-1)*nx] * 0.1f;
    tmp_image[i*nx] += image[(i+1)*nx] * 0.1f;
    tmp_image[i*nx] += image[1+i*nx] * 0.1f;
    //#pragma vector always
    for (int j = 1; j < ny-1; ++j) {
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1f;
    }
    //when j=nx-1
    tmp_image[(ny-1)+i*nx] = image[(ny-1)+i*nx] * 0.6f;
    tmp_image[(ny-1)+i*nx] += image[(ny-1)  +(i-1)*nx] * 0.1f;
    tmp_image[(ny-1)+i*nx] += image[(ny-1)  +(i+1)*nx] * 0.1f;
    tmp_image[(ny-1)+i*nx] += image[(ny-1)-1+i*nx] * 0.1f;
  }
  //when i = ny-1
  //when j=0
  tmp_image[(nx-1)*nx] = image[(nx-1)*nx] * 0.6f;
  tmp_image[(nx-1)*nx] += image[((nx-1)-1)*nx] * 0.1f;
  tmp_image[(nx-1)*nx] += image[1+(nx-1)*nx] * 0.1f;
  //#pragma vector always
  for (int j = 1; j < ny-1; ++j) {
    tmp_image[j+(nx-1)*nx] = image[j+(nx-1)*nx] * 0.6f;
    tmp_image[j+(nx-1)*nx] += image[j  +((nx-1)-1)*nx] * 0.1f;
    tmp_image[j+(nx-1)*nx] += image[j-1+(nx-1)*nx] * 0.1f;
    tmp_image[j+(nx-1)*nx] += image[j+1+(nx-1)*nx] * 0.1f;
  }
  //when j=nx-1
  tmp_image[(ny-1)+(nx-1)*nx] = image[(ny-1)+(nx-1)*nx] * 0.6f;
  tmp_image[(ny-1)+(nx-1)*nx] += image[(ny-1)  +((nx-1)-1)*nx] * 0.1f;
  tmp_image[(ny-1)+(nx-1)*nx] += image[(ny-1)-1+(nx-1)*nx] * 0.1f;
}

//used when a partition is the top of the stencil
void stenciltop(const short nx, const short bottom, float * restrict image, float * restrict tmp_image){
  //when i=0
  //when j=0
  tmp_image[0] = image[0] * 0.6f;
  tmp_image[0] += image[nx] * 0.1f;
  tmp_image[0] += image[0] * 0.1f;
  //#pragma vector always
  for (int j = 1; j < nx-1; ++j) {
    tmp_image[j] = image[j] * 0.6f;
    tmp_image[j] += image[j  +nx] * 0.1f;
    tmp_image[j] += image[j-1] * 0.1f;
    tmp_image[j] += image[j+1] * 0.1f;
  }
  //when j=nx-1
  tmp_image[nx-1] = image[nx-1] * 0.6f;
  tmp_image[nx-1] += image[(nx-1) +nx] * 0.1f;
  tmp_image[nx-1] += image[(nx-1)-1] * 0.1f;



  for (int i = 1; i <= bottom; ++i) {
    //when j=0
    tmp_image[i*nx] = image[i*nx] * 0.6f;
    tmp_image[i*nx] += image[(i-1)*nx] * 0.1f;
    tmp_image[i*nx] += image[(i+1)*nx] * 0.1f;
    tmp_image[i*nx] += image[1+i*nx] * 0.1f;
    for (int j = 1; j < nx-1; ++j) {
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1f;
    }
    //when j=nx-1
    tmp_image[(nx-1)+i*nx] = image[(nx-1)+i*nx] * 0.6f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)  +(i-1)*nx] * 0.1f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)  +(i+1)*nx] * 0.1f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)-1+i*nx] * 0.1f;
  }

}

void stencilbottom(const short nx, const short top, const short bottom, float * restrict image, float * restrict tmp_image){
  for (int i = top; i < bottom; ++i) {
    //when j=0
    tmp_image[i*nx] = image[i*nx] * 0.6f;
    tmp_image[i*nx] += image[(i-1)*nx] * 0.1f;
    tmp_image[i*nx] += image[(i+1)*nx] * 0.1f;
    tmp_image[i*nx] += image[1+i*nx] * 0.1f;
    for (int j = 1; j < nx-1; ++j) {
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1f;
    }
    //when j=nx-1
    tmp_image[(nx-1)+i*nx] = image[(nx-1)+i*nx] * 0.6f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)  +(i-1)*nx] * 0.1f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)  +(i+1)*nx] * 0.1f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)-1+i*nx] * 0.1f;
  }


  //when j=0
  tmp_image[bottom*nx] = image[bottom*nx] * 0.6f;
  tmp_image[bottom*nx] += image[(bottom-1)*nx] * 0.1f;
  tmp_image[bottom*nx] += image[1+bottom*nx] * 0.1f;

  for (int j = 1; j < nx-1; ++j) {
    tmp_image[j+bottom*nx] = image[j+bottom*nx] * 0.6f;
    tmp_image[j+bottom*nx] += image[j  +(bottom-1)*nx] * 0.1f;
    tmp_image[j+bottom*nx] += image[j-1+bottom*nx] * 0.1f;
    tmp_image[j+bottom*nx] += image[j+1+bottom*nx] * 0.1f;
  }
  //when j = nx-1
  tmp_image[(nx-1)+bottom*nx] = image[(nx-1)+bottom*nx] * 0.6f;
  tmp_image[(nx-1)+bottom*nx] += image[(nx-1)  +(bottom-1)*nx] * 0.1f;
  tmp_image[(nx-1)+bottom*nx] += image[(nx-1)-1+bottom*nx] * 0.1f;
}



//middle partitions
void stencilmiddle(const short nx, const short top, const short bottom, float * restrict image, float * restrict tmp_image){
  for (int i = top; i <= bottom; ++i) {
    //when j=0
    tmp_image[i*nx] = image[i*nx] * 0.6f;
    tmp_image[i*nx] += image[(i-1)*nx] * 0.1f;
    tmp_image[i*nx] += image[(i+1)*nx] * 0.1f;
    tmp_image[i*nx] += image[1+i*nx] * 0.1f;
    for (int j = 1; j < nx-1; ++j) {
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1f;
    }
    //when j=nx-1
    tmp_image[(nx-1)+i*nx] = image[(nx-1)+i*nx] * 0.6f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)  +(i-1)*nx] * 0.1f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)  +(i+1)*nx] * 0.1f;
    tmp_image[(nx-1)+i*nx] += image[(nx-1)-1+i*nx] * 0.1f;
  }
}



// Create the input image
void init_image(const short nx, const short ny, float * restrict image, float * restrict tmp_image) {
  // Zero everything
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      image[j+i*ny] = 0.0f;
      tmp_image[j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int ii = j*nx/8; ii < (j+1)*nx/8; ++ii) {
        for (int jj = i*ny/8; jj < (i+1)*ny/8; ++jj) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const short nx, const short ny, float *image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (image[j+i*ny] > maximum)
      maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      fputc((char)(255.0f*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}
//calculate the boundaries of the given row
short calculateRowBoundary(const short ny, const short rank, const short noOfPartitions, short topOrBottom){
  short output;
  if(topOrBottom == LOW){
    output = (short)((ny * rank) / noOfPartitions);
  }
  if(topOrBottom == HIGH){
    output = (short)((ny * (rank+1)) / noOfPartitions) - 1;
  }
  return output;
}
/*
short calculateAdjacentRank(const short direction, const short rank, const short noOfProcesses){
short output;
if(direction == ABOVE){
if(rank == 0){
output = noOfProcesses - 1;
}
else{
output = rank - 1;
}
}
else if(direction == BELOW){
if(rank = noOfProcesses - 1){
output = 0;
}
else{
output = rank + 1;
}
}
return output;
}*/

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
