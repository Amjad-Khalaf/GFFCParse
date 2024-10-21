#include "file.h"

void read_file(const char *filename) {
    // Declare file pointer
    FILE *file;

    // Open file in read mode
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Cannot open file!");
        return;
    }

    // Initial buffer size for the line
    size_t line_size = 128;

    // Allocate memory for the line buffer
    char *line = (char *)malloc(line_size);
    if (line == NULL) {
        perror("Failed to allocate memory!");
        fclose(file);  // Close the file before exiting
        return;
    }

    // Initialize to track how many characters are in the line
    int line_length = 0;

    // Maximum size of the chunk read in each iteration
    int chunk_length = 128;

    // Temporary buffer to store the chunk read from the file
    char chunk[128];

    // Read file in chunks and process
    while (fgets(chunk, chunk_length, file)) {
        // Check if the chunk fits in the allocated line buffer
        size_t chunk_len = strlen(chunk);
        if (chunk_len >= (line_size - line_length)) {
            line_size *= 2;  // Double the buffer size
            line = realloc(line, line_size);
            if (line == NULL) {
                perror("Failed to reallocate memory!");
                fclose(file);  // Close file before exiting
                return;
            }
        }

        // Append the chunk to the line buffer
        strcat(line, chunk);
        line_length += chunk_len;

        // If the line ends with a newline character
        if (line[line_length - 1] == '\n') {
            printf("%s", line);

            // Reset line_length to start fresh for the next line
            line[0] = '\0';  // Reset the line buffer to empty
            line_length = 0;
        }
    }

    // Free allocated memory
    free(line);

    // Close the file
    fclose(file);
}
