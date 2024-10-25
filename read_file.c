// read_file.c

#include "read_file.h"
#include <string.h>
#include <stdio.h>

// Function to split a line based on a given delimiter
struct line split_line(char *line, char delimiter) {
    // Declare array to store all the elements
    int number_of_elements = 0;
    char **split_line = NULL; // Initialize as NULL

    // Declare array to store each element
    int size_of_element = 1; // Accommodate for null terminator
    char *element = (char*)malloc(size_of_element * sizeof(char));
    
    // Iterate through line to find delimiter
    int position_in_line;
    int position_in_element = 0;
    for (position_in_line = 0; position_in_line < strlen(line); position_in_line++) {
        
        if (line[position_in_line] != delimiter) { // Append to element
            size_of_element++;
            element = realloc(element, size_of_element);
            element[position_in_element] = line[position_in_line];
            position_in_element++;
        }
        
        if (line[position_in_line] == delimiter) {
            element[position_in_element] = '\0'; // Add null terminator
            
            // Increase size of split_line to accommodate element
            number_of_elements++;
            split_line = realloc(split_line, number_of_elements * sizeof(char*));
            split_line[number_of_elements - 1] = strdup(element); // Append a copy of element
            
            // Reset position in individual element
            size_of_element = 1;
            element = (char*)malloc(size_of_element * sizeof(char));
            position_in_element = 0;
        }
    }

    // Handle the last element if there's no trailing delimiter
    if (position_in_element > 0) {
        element[position_in_element] = '\0'; // Null-terminate the last element
        number_of_elements++;
        split_line = realloc(split_line, number_of_elements * sizeof(char*));
        split_line[number_of_elements - 1] = strdup(element); // Append last element
    }

    free(element); // Free the last temporary element

    // Return result as line struct
    struct line result;
    result.array_of_elements = split_line;
    result.number_of_elements = number_of_elements;

    return result;
}

// Function to remove trailing newline
char* remove_trailing_newline(char *line) {
    line[strcspn(line, "\n")] = 0;
    return line;
}