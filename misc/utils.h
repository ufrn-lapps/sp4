#ifndef READ_2D_H
#define READ_2D_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * @brief  Read file into float array.
 * @note Expect metadata with array dimensions (size_x size_y)
 * @param  *filename: File name to be read.
 * @retval Array float where the data will be allocated.
 */
float *read_2d(const char *filename) {
    FILE *file;
    int size_x, size_y;
    float *arr;

    file = fopen(filename, "r");

    if (!fscanf(file, "%d", &size_y) || !fscanf(file, "%d", &size_x)) {
        printf("Error reading file %s.\n", filename);
    }

    // printf("%s\t size: %d %d\n", filename, size_y, size_x);

    arr = (float *)malloc(sizeof(float) * size_x * size_y);

    for (int j = 0; j < size_y; j++) {
        for (int i = 0; i < size_x; i++) {
            if (!fscanf(file, "%e", &arr[j * size_x + i])) {
                printf("Error reading file %s.\n", filename);
            }
        }
    }

    return arr;
}

/**
 * @brief  Write data structure into a file.
 * @note   If file exist it will be overwritten.
 * @param  *filename: Name of the file to be saved.
 * @param  *mydata: Data to be written in the file. Float array expected.
 * @param  metadata: struct with size of the data.
 * @retval None
 */
void array2d_to_file(char *filename, float *mydata, int size_x, int size_y) {
    FILE *fptr;
    fptr = fopen(filename, "w");

    if (fptr == NULL) {
        printf("Error!");
        exit(1);
    }

    // Write meta data into the file.
    printf("Filename=%s\tSize=%dx%d\n", filename, size_y, size_x);
    fprintf(fptr, "%d %d\n", size_x, size_y);

    // Write data.
    for (int j = 0; j < size_x; j++) {
        for (int i = 0; i < size_y; i++) {
            fprintf(fptr, "%e ", mydata[j * size_y + i]);
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    return;
}

#endif