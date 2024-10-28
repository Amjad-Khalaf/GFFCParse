// gene_count_hash.h

#ifndef GENE_COUNT_HASH_H
#define GENE_COUNT_HASH_H

#include <stdlib.h>  // For malloc
#include <string.h>  // For strdup

#define HASH_TABLE_SIZE 15000

struct gene_count {
    char *gene_name;
    int count;
    struct gene_count *next;  // Stores collisions in a linked list
};

// Global hash table to store gene counts
extern struct gene_count *gene_table[HASH_TABLE_SIZE];

// Function prototypes
unsigned int hash_function(const char *gene_name);
void increment_gene_count(const char *gene_name);
int get_gene_count(const char *gene_name);
int get_total_gene_count();

#endif  // GENE_COUNT_HASH_H
