#include <stdlib.h>           // For qsort
#include "exon_operations.h"   // For function prototypes and access to genome features

// Helper function for sorting, defined as static since it will only be used internally
static int compare_exons(const void *a, const void *b) {    
    struct exon *exon1 = *(struct exon **)a;
    struct exon *exon2 = *(struct exon **)b;
    return exon1->start - exon2->start;
}

void sort_exons(struct mRNA mRNA) {    
    qsort(mRNA.array_of_exons, mRNA.number_of_exons, sizeof(struct exon *), compare_exons);
}

void total_exon_length(struct mRNA *mRNA) {
    int total_length = 0;
    for (int i = 0; i < mRNA->number_of_exons; i++) {
        total_length += (mRNA->array_of_exons[i]->end - mRNA->array_of_exons[i]->start);
    }
    mRNA->total_exon_length = total_length; // Assign total length to mRNA struct
}

void calculate_exon_stats(struct mRNA *mRNA) {
    int total_length = 0; 
    int largest_exon_length = 0;

    for (int i = 0; i < mRNA->number_of_exons; i++) {
        int length = (mRNA->array_of_exons[i]->end - mRNA->array_of_exons[i]->start);
        total_length += length;
        if (length > largest_exon_length) {
            largest_exon_length = length;
        }
    }

    mRNA->total_exon_length = total_length;
    mRNA->average_exon_length = mRNA->number_of_exons > 0 ? total_length / mRNA->number_of_exons : 0;
    mRNA->largest_exon_length = largest_exon_length;
}