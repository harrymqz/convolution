//
//  convolution.c
//
//
//  Created by Josep Lluis Lerida on 11/03/15.
//
// This program calculates the convolution for PPM images.
// The program accepts an PPM image file, a text definition of the kernel matrix and the PPM file for storing the convolution results.
// The program allows to define image partitions for processing large images (>500MB)
// The 2D image is represented by 1D vector for chanel R, G and B. The convolution is applied to each chanel separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define MAX_THREADS 3

// Estructura per emmagatzemar el contingut d'una imatge.

struct imagenppm {
    int altura;
    int ancho;
    char *comentario;
    int maxcolor;
    int P;
    int *R;
    int *G;
    int *B;
};
typedef struct imagenppm* ImagenData;

// Estructura per emmagatzemar el contingut d'un kernel.

struct structkernel {
    int kernelX;
    int kernelY;
    float *vkern;
};
typedef struct structkernel* kernelData;

//Functions Definition
ImagenData initimage(char* nombre, FILE **fp, int partitions, int halo);
ImagenData duplicateImageData(ImagenData src, int partitions, int halo);

int readImage(ImagenData Img, FILE **fp, int dim, int halosize, long int *position);
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim);
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position);
int savingChunk(ImagenData img, FILE **fp, int dim, int offset);
int convolve2D(int* inbuf, int* outbuf, int sizeX, int sizeY, float* kernel, int ksizeX, int ksizeY);
void freeImagestructure(ImagenData *src);

//----------------------------- MPI----------------------------
int rank, nprocs;

void sendTask(ImagenData source, int chunksize, int partsize, int partitions, int offset, int dest);
void sendToMaster(ImagenData output, int chunksize, int partsize, int offset, int dest);
void receiveTask(ImagenData source, ImagenData output, int* chunksize, int* partsize, int* height, int* offset, int from);
void receiveFromMaster(ImagenData output, int* chunksize, int* partsize, int* offset, int from);
//-------------------------------------------------------------

//Open Image file and image struct initialization

ImagenData initimage(char* nombre, FILE **fp, int partitions, int halo) {
    char c;
    char comentario[300];
    int i = 0, chunk = 0;
    ImagenData img = NULL;

    /*Se habre el fichero ppm*/

    if ((*fp = fopen(nombre, "r")) == NULL) {
        perror("Error: ");
    } else {
        //Memory allocation
        img = (ImagenData) malloc(sizeof (struct imagenppm));

        //Reading the first line: Magical Number "P3"
        fscanf(*fp, "%c%d ", &c, &(img->P));

        //Reading the image comment
        while ((c = fgetc(*fp)) != '\n') {
            comentario[i] = c;
            i++;
        }
        comentario[i] = '\0';

        //Allocating information for the image comment
        img->comentario = calloc(strlen(comentario), sizeof (char));
        strcpy(img->comentario, comentario);

        //Reading image dimensions and color resolution
        fscanf(*fp, "%d %d %d", &img->ancho, &img->altura, &img->maxcolor);
        chunk = img->ancho * img->altura / partitions;

        //We need to read an extra row.
        chunk = chunk + img->ancho * halo;
        if ((img->R = calloc(chunk, sizeof (int))) == NULL) {
            return NULL;
        }
        if ((img->G = calloc(chunk, sizeof (int))) == NULL) {
            return NULL;
        }
        if ((img->B = calloc(chunk, sizeof (int))) == NULL) {
            return NULL;
        }
    }
    return img;
}

//Duplicate the Image struct for the resulting image

ImagenData duplicateImageData(ImagenData src, int partitions, int halo) {
    char c;
    char comentario[300];
    unsigned int imageX, imageY;
    int i = 0, chunk = 0;

    //Struct memory allocation
    ImagenData dst = (ImagenData) malloc(sizeof (struct imagenppm));

    //Copying the magic number
    dst->P = src->P;

    //Copying the string comment
    dst->comentario = calloc(strlen(src->comentario), sizeof (char));
    strcpy(dst->comentario, src->comentario);

    //Copying image dimensions and color resolution
    dst->ancho = src->ancho;
    dst->altura = src->altura;
    dst->maxcolor = src->maxcolor;
    chunk = dst->ancho * dst->altura / partitions;

    //We need to read an extra row.
    chunk = chunk + src->ancho * halo;
    if ((dst->R = calloc(chunk, sizeof (int))) == NULL) {
    }
    if ((dst->G = calloc(chunk, sizeof (int))) == NULL) {
        return NULL;
    }
    if ((dst->B = calloc(chunk, sizeof (int))) == NULL) {
        return NULL;
    }

    return dst;
}

//Read the corresponding chunk from the source Image

int readImage(ImagenData img, FILE **fp, int dim, int halosize, long *position) {
    // printf("\t dim: %d, halosize: %d, position %f\n", dim, halosize, *position);

    int i = 0, haloposition = 0;
    if (fseek(*fp, *position, SEEK_SET))
        perror("Error: ");

    haloposition = dim - (img->ancho * halosize * 2);

    // printf("\t dim: %d\n", dim);
    for (i = 0; i < dim; i++) {
        // printf("\t 1. Forsins --> %d --> %d and rank: %d\n", i, i, rank);

        // When start reading the halo store the position in the image file
        if (halosize != 0 && i == haloposition) *position = ftell(*fp);

        fscanf(*fp, "%d %d %d ", &img->R[i], &img->G[i], &img->B[i]);
    }
    // printf("\t Read 5. image\n");

    //    printf ("Readed = %d pixels, posicio=%lu\n",k,*position);
    return 0;
}

//Duplication of the  just readed source chunk to the destiny image struct chunk

int duplicateImageChunk(ImagenData src, ImagenData dst, int dim) {
    int i = 0;

    for (i = 0; i < dim; i++) {
        dst->R[i] = src->R[i];
        dst->G[i] = src->G[i];
        dst->B[i] = src->B[i];
    }

    //    printf ("Duplicated = %d pixels\n",i);
    return 0;
}

// Open kernel file and reading kernel matrix. The kernel matrix 2D is stored in 1D format.

kernelData leerKernel(char* nombre) {
    FILE *fp;
    int i = 0;
    kernelData kern = NULL;

    /*Opening the kernel file*/
    fp = fopen(nombre, "r");
    if (!fp) {
        perror("Error: ");
    } else {
        //Memory allocation
        kern = (kernelData) malloc(sizeof (struct structkernel));

        //Reading kernel matrix dimensions
        fscanf(fp, "%d,%d,", &kern->kernelX, &kern->kernelY);
        kern->vkern = (float *) malloc(kern->kernelX * kern->kernelY * sizeof (float));

        // Reading kernel matrix values
        for (i = 0; i < (kern->kernelX * kern->kernelY) - 1; i++) {
            fscanf(fp, "%f,", &kern->vkern[i]);
        }
        fscanf(fp, "%f", &kern->vkern[i]);
        fclose(fp);
    }

    return kern;
}

// Open the image file with the convolution results

int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position) {
    /*Se crea el fichero con la imagen resultante*/
    if ((*fp = fopen(nombre, "w")) == NULL) {
        perror("Error: ");
        return -1;
    }

    /*Writing Image Header*/
    fprintf(*fp, "P%d\n%s\n%d %d\n%d\n", img->P, img->comentario, img->ancho, img->altura, img->maxcolor);
    *position = ftell(*fp);
    return 0;
}

// Writing the image partition to the resulting file. dim is the exact size to write. offset is the displacement for avoid halos.

int savingChunk(ImagenData img, FILE **fp, int dim, int offset) {
    int i, k = 0;

    //Writing image partition
    for (i = offset; i < dim + offset; i++) {
        fprintf(*fp, "%d %d %d ", img->R[k], img->G[k], img->B[k]);
        //        if ((i+1)%6==0) fprintf(*fp,"\n");
        k++;
    }

    //    printf ("Writed = %d pixels, dim=%d, offset=%d\n",k,dim, offset);
    return 0;
}

// This function free the space allocated for the image structure.

void freeImagestructure(ImagenData *src) {

    free((*src)->comentario);
    free((*src)->R);
    free((*src)->G);
    free((*src)->B);

    free(*src);
}

///////////////////////////////////////////////////////////////////////////////
// 2D convolution
// 2D data are usually stored in computer memory as contiguous 1D array.
// So, we are using 1D array for 2D data.
// 2D convolution assumes the kernel is center originated, which means, if
// kernel size 3 then, k[-1], k[0], k[1]. The middle of index is always 0.
// The following programming logics are somewhat complicated because of using
// pointer indexing in order to minimize the number of multiplications.
//
//
// signed integer (32bit) version:
///////////////////////////////////////////////////////////////////////////////

int convolve2D(int* in, int* out, int dataSizeX, int dataSizeY,
        float* kernel, int kernelSizeX, int kernelSizeY) {
    int i, j, m, n;
    int *inPtr, *inPtr2, *outPtr;
    float *kPtr;
    int kCenterX, kCenterY;
    int rowMin, rowMax; // to check boundary of input array
    int colMin, colMax; //
    float sum; // temp accumulation buffer

    // check validity of params
    if (!in || !out || !kernel) return -1;
    if (dataSizeX <= 0 || kernelSizeX <= 0) return -1;

    // find center position of kernel (half of kernel size)
    kCenterX = (int) kernelSizeX / 2;
    kCenterY = (int) kernelSizeY / 2;

    // init working  pointers
    inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX]; // note that  it is shifted (kCenterX, kCenterY),
    outPtr = out;
    kPtr = kernel;

    // start convolution
    for (i = 0; i < dataSizeY; ++i) // number of rows
    {
        // compute the range of convolution, the current row of kernel should be between these
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;

        for (j = 0; j < dataSizeX; ++j) // number of columns
        {
            // compute the range of convolution, the current column of kernel should be between these
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;

            sum = 0; // set to 0 before accumulate

            // flip the kernel and traverse all the kernel values
            // multiply each kernel value with underlying input data
            for (m = 0; m < kernelSizeY; ++m) // kernel rows
            {
                // check if the index is out of bound of input array
                if (m <= rowMax && m > rowMin) {
                    for (n = 0; n < kernelSizeX; ++n) {
                        // check the boundary of array
                        if (n <= colMax && n > colMin)
                            sum += *(inPtr - n) * *kPtr;

                        ++kPtr; // next kernel
                    }
                } else
                    kPtr += kernelSizeX; // out of bound, move to next row of kernel

                inPtr -= dataSizeX; // move input data 1 raw up
            }

            // convert integer number
            if (sum >= 0) *outPtr = (int) (sum + 0.5f);
                //            else *outPtr = (int)(sum - 0.5f)*(-1);
                // For using with image editors like GIMP or others...
            else *outPtr = (int) (sum - 0.5f);
            // For using with a text editor that read ppm images like libreoffice or others...
            //            else *outPtr = 0;

            kPtr = kernel; // reset kernel to (0,0)
            inPtr = ++inPtr2; // next input
            ++outPtr; // next output
        }
    }

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
    //MPI
    //-------------------------------------------------------------------------
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //-------------------------------------------------------------------------

    int i = 0, j = 0, k = 0;
    //    int headstored=0, imagestored=0, stored;

    if (argc != 5) {
        printf("Usage: %s <image-file> <kernel-file> <result-file> <partitions>\n", argv[0]);

        printf("\n\nError, Missing parameters:\n");
        printf("format: ./serialconvolution image_file kernel_file result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n");
        printf("- partitions : Image partitions\n\n");
        return -1;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // READING IMAGE HEADERS, KERNEL Matrix, DUPLICATE IMAGE DATA, OPEN RESULTING IMAGE FILE
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // if (rank == 0)
    // {
    int imagesize, partitions, partsize, chunksize, halo, halosize;
    long position = 0;
    double start, tstart = 0, tend = 0, tread = 0, tcopy = 0, tconv = 0, tstore = 0, treadk = 0;
    struct timeval tim;
    FILE *fpsrc = NULL, *fpdst = NULL;
    ImagenData source = NULL, output = NULL;
    // }

    // Store number of partitions

    // if (rank == 0)
    // {
    partitions = atoi(argv[4]);
    // }

    //Reading kernel matrix
    start = MPI_Wtime();
    tstart = start;
    kernelData kern = NULL;

    if ((kern = leerKernel(argv[2])) == NULL) {
        //        free(source);
        //        free(output);
        return -1;
    }



    //The matrix kernel define the halo size to use with the image. The halo is zero when the image is not partitioned.
    if (partitions == 1) halo = 0;
    else halo = (kern->kernelY / 2)*2;

    treadk = treadk + (MPI_Wtime() - start);


    //Reading Image Header. Image properties: Magical number, comment, size and color resolution.
    if (rank == 0) {
        start = MPI_Wtime();

        //Memory allocation based on number of partitions and halo size.
        if ((source = initimage(argv[1], &fpsrc, partitions, halo)) == NULL) {
            return -1;
        }

        printf("IMAGE SIZE X = %d, Y = %d\n", source->ancho, source->altura);

        tread = tread + (MPI_Wtime() - start);

        //Duplicate the image struct.
        start = MPI_Wtime();

        if ((output = duplicateImageData(source, partitions, halo)) == NULL) {
            return -1;
        }

        tcopy = tcopy + (MPI_Wtime() - start);


        //Initialize Image Storing file. Open the file and store the image header.
        start = MPI_Wtime();
        if (initfilestore(output, &fpdst, argv[3], &position) != 0) {
            perror("Error: ");
            //        free(source);
            //        free(output);
            return -1;
        }
        tstore = tstore + (MPI_Wtime() - start);
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////
    // CHUNK READING
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if (rank == 0) {
        int c = 0, offset = 0;
        imagesize = source->altura * source->ancho;
        partsize = (source->altura * source->ancho) / partitions;


        while (c < partitions) {
            int i, mpi_partsize, mpi_partsize_module;
            for (i = 1; i < nprocs; i++) {
                ////////////////////////////////////////////////////////////////////////////////
                //Reading Next chunk.
                start = MPI_Wtime();
                if (c == 0) {
                    halosize = halo / 2;
                    chunksize = partsize + (source->ancho * halosize);
                    offset = 0;
                } else if (c < nprocs - 1) {
                    halosize = halo;
                    chunksize = partsize + (source->ancho * halosize);
                    offset = (source->ancho * halo / 2);
                } else {
                    halosize = halo / 2;
                    chunksize = partsize + (source->ancho * halosize);
                    offset = (source->ancho * halo / 2);
                }
                if (readImage(source, &fpsrc, chunksize, halo / 2, &position)) {
                    return -1;
                }
                tread = tread + (MPI_Wtime() - start);
                //Duplicate the image chunk
                start = MPI_Wtime();
                if (duplicateImageChunk(source, output, chunksize)) {
                    return -1;
                }

                tcopy = tcopy + (MPI_Wtime() - start);
                mpi_partsize = chunksize / (nprocs - 1);
                mpi_partsize_module = chunksize % (nprocs - 1);
                int height = (source->altura / partitions) + halosize;
                int begin = mpi_partsize * (i - 1);
                int array_length = mpi_partsize;

                // If is the last process we have to add the mudule to not lost anything
                if (i == nprocs - 1) {
                    array_length = mpi_partsize + mpi_partsize_module;
                }
                MPI_Status status;
                MPI_Send(&array_length, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&partsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&offset, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&height, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&source->ancho, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

                MPI_Send(&source->R[begin], array_length, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&source->G[begin], array_length, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&source->B[begin], array_length, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            tconv = MPI_Wtime();
            ImagenData output;
            int array_length = 0, partsize = 0, offset = 0;
            MPI_Status status;
            for (i = 1; i < nprocs; i++) {
                int begin = mpi_partsize * (i - 1);
                array_length = mpi_partsize;
                if (i == nprocs - 1) {
                    array_length = mpi_partsize + mpi_partsize_module;
                }
                MPI_Recv(&array_length, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                output = (ImagenData) malloc(sizeof (struct imagenppm));
                output->R = calloc(array_length, sizeof (int));
                output->G = calloc(array_length, sizeof (int));
                output->B = calloc(array_length, sizeof (int));
                MPI_Recv(&partsize, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&offset, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(output->R, array_length, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(output->G, array_length, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(output->B, array_length, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                gettimeofday(&tim, NULL);
                start = MPI_Wtime();
                if (savingChunk(output, &fpdst, array_length, begin)) {
                    perror("Error: ");
                    return -1;
                }
                gettimeofday(&tim, NULL);
                tstore = tstore + (MPI_Wtime() - start);
                //                }
            }
            c++;
        }
        tconv = tconv- MPI_Wtime();
        int i;
        int flag = -1;
        for (i = 1; i < nprocs; i++) {
            MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // CHUNK CONVOLUTION
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if (rank != 0) {
        ImagenData source, output;
        int array_length, height, partsize, offset;
        MPI_Status status;
        while (1) {
            MPI_Recv(&array_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            if (array_length == -1) {
                MPI_Finalize();
            }
            source = (ImagenData) malloc(sizeof (struct imagenppm));
            source->R = malloc((int) array_length  * sizeof (int));
            source->G = malloc((int) array_length  * sizeof (int));
            source->B = malloc((int) array_length  * sizeof (int));
            MPI_Recv(&partsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&source->ancho, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(source->R, array_length, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(source->G, array_length, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(source->B, array_length, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            output = (ImagenData) malloc(sizeof (struct imagenppm));
            output->R = malloc((int) array_length  * sizeof (int));
            output->G = malloc((int) array_length   * sizeof (int));
            output->B = malloc((int) array_length  * sizeof (int));
            #pragma omp parallel num_threads(MAX_THREADS)
            {
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        convolve2D(source->R, output->R, source->ancho, height/(nprocs-1), kern->vkern, kern->kernelX, kern->kernelY);
                    }
                    #pragma omp section
                    {
                        convolve2D(source->G, output->G, source->ancho, height/(nprocs-1), kern->vkern, kern->kernelX, kern->kernelY);
                    }
                    #pragma omp section
                    {
                        convolve2D(source->B, output->B, source->ancho, height/(nprocs-1), kern->vkern, kern->kernelX, kern->kernelY);
                    }
                }
            }
            MPI_Send(&array_length, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&partsize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&offset, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&output->R[0], array_length, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&output->G[0], array_length, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&output->B[0], array_length, MPI_INT, 0, 1, MPI_COMM_WORLD);
            free(source->R);
            free(source->G);
            free(source->B);
            free(source);
            free(output->R);
            free(output->G);
            free(output->B);
            free(output);
        }
        MPI_Finalize();
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // CHUNK SAVING
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Storing resulting image partition.

    if (rank == 0) {
        fclose(fpsrc);
        fclose(fpdst);
        gettimeofday(&tim, NULL);
        tend = MPI_Wtime();
        printf("Imatge: %s\n", argv[1]);
        printf("ISizeX : %d\n", source->ancho);
        printf("ISizeY : %d\n", source->altura);
        printf("kSizeX : %d\n", kern->kernelX);
        printf("kSizeY : %d\n", kern->kernelY);
        printf("%.6lf seconds elapsed for Reading image file.\n", tread);
        printf("%.6lf seconds elapsed for copying image structure.\n", tcopy);
        printf("%.6lf seconds elapsed for Reading kernel matrix.\n", treadk);
        printf("%.6lf seconds elapsed for make the convolution.\n", tconv);
        printf("%.6lf seconds elapsed for writing the resulting image.\n", tstore);
        printf("%.6lf seconds elapsed\n", tend - tstart);
        freeImagestructure(&source);
        freeImagestructure(&output);
    }
    return 0;
}


