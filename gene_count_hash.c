#include "gene_count_hash.h"

// Global hash table to store gene counts
struct gene_count *gene_table[HASH_TABLE_SIZE] = {NULL};

unsigned int hash_function(const char *gene_name) {
    unsigned int hash = 0;
    while (*gene_name) {
        hash = (hash << 5) + *gene_name++;
    }
    return hash % HASH_TABLE_SIZE;
}

void increment_gene_count(const char *gene_name) {
    unsigned int hash_index = hash_function(gene_name);
    struct gene_count *entry = gene_table[hash_index];

    // Check if the gene is already in the table
    while (entry != NULL) {
        if (strcmp(entry->gene_name, gene_name) == 0) {
            entry->count++;
            return;
        }
        entry = entry->next;
    }

    // If not found, add a new entry
    struct gene_count *new_entry = (struct gene_count *)malloc(sizeof(struct gene_count));
    new_entry->gene_name = strdup(gene_name); 
    new_entry->count = 1;
    new_entry->next = gene_table[hash_index];
    gene_table[hash_index] = new_entry;
}

int get_gene_count(const char *gene_name) {
    unsigned int index = hash_function(gene_name);
    struct gene_count *entry = gene_table[index];

    // Traverse the linked list for this hash index
    while (entry != NULL) {
        if (strcmp(entry->gene_name, gene_name) == 0) {
            return entry->count;  // Return the count if found
        }
        entry = entry->next;
    }

    return 0;  // If the gene is not found, return 0
}

int get_total_gene_count() {
    int total_gene_count = 0;

    // Loop through each index of the hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        struct gene_count *entry = gene_table[i];

        // Traverse the linked list at this index
        while (entry != NULL) {
            total_gene_count += entry->count;  // Add the count for this gene
            entry = entry->next;               // Move to the next node in the list
        }
    }

    return total_gene_count;
}
