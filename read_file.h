// read_file.h

#ifndef READ_FILE_H
#define READ_FILE_H

#include <stdlib.h>

// Define the struct
struct line {
    char** array_of_elements;
    int number_of_elements;
};

// Function prototypes
struct line split_line(char *line, char delimiter);
char* remove_trailing_newline(char *line);

#endif // READ_FILE_H
