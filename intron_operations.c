#include <stdlib.h>           // For malloc and realloc
#include "intron_operations.h"  // For function prototypes and access to genome features

void identify_introns(struct mRNA *mRNA) {
    mRNA->number_of_introns = 0;
    struct intron **array_of_introns = NULL;

    // Identify first intron
    if ((mRNA->array_of_exons[0]->start - mRNA->start) > 0) {
        mRNA->number_of_introns++;
        array_of_introns = realloc(array_of_introns, sizeof(struct intron*) * mRNA->number_of_introns);
        struct intron *new_intron = malloc(sizeof(struct intron));
        new_intron->start = mRNA->start;
        new_intron->end = mRNA->array_of_exons[0]->start - 1;
        new_intron->length = (mRNA->array_of_exons[0]->start - 1) - mRNA->start;
        array_of_introns[mRNA->number_of_introns - 1] = new_intron;
    }

    // Identify inter-exonic introns
    for (int i = 0; i < mRNA->number_of_exons - 1; i++) {
        if ((mRNA->array_of_exons[i + 1]->start - mRNA->array_of_exons[i]->end) > 0) {
            mRNA->number_of_introns++;
            array_of_introns = realloc(array_of_introns, sizeof(struct intron*) * mRNA->number_of_introns);
            struct intron *new_intron = malloc(sizeof(struct intron));
            new_intron->start = mRNA->array_of_exons[i]->end + 1;
            new_intron->end = mRNA->array_of_exons[i + 1]->start - 1;
            new_intron->length = (mRNA->array_of_exons[i + 1]->start - 1) - (mRNA->array_of_exons[i]->end + 1);
            array_of_introns[mRNA->number_of_introns - 1] = new_intron;
        }
    }

    // Identify final intron
    if ((mRNA->end - mRNA->array_of_exons[mRNA->number_of_exons - 1]->end) > 0) {
        mRNA->number_of_introns++;
        array_of_introns = realloc(array_of_introns, sizeof(struct intron*) * mRNA->number_of_introns);
        struct intron *new_intron = malloc(sizeof(struct intron));
        new_intron->start = mRNA->array_of_exons[mRNA->number_of_exons - 1]->end + 1;
        new_intron->end = mRNA->end;
        new_intron->length = mRNA->end - (mRNA->array_of_exons[mRNA->number_of_exons - 1]->end + 1);
        array_of_introns[mRNA->number_of_introns - 1] = new_intron;
    }

    // Assign the intron array to the mRNA struct
    mRNA->array_of_introns = array_of_introns;
}

void total_intron_length(struct mRNA *mRNA) {
    int total_length = 0;
    for (int i = 0; i < mRNA->number_of_introns; i++) {
        total_length += (mRNA->array_of_introns[i]->end - mRNA->array_of_introns[i]->start);
    }
    mRNA->total_intron_length = total_length; // Assign total length to mRNA struct
}

void calculate_intron_stats(struct mRNA *mRNA) {
    int total_length = 0; 
    int largest_intron_length = 0;

    for (int i = 0; i < mRNA->number_of_introns; i++) {
        int length = (mRNA->array_of_introns[i]->end - mRNA->array_of_introns[i]->start);
        total_length += length;
        if (length > largest_intron_length) {
            largest_intron_length = length;
        }
    }

    mRNA->total_intron_length = total_length;
    mRNA->average_intron_length = mRNA->number_of_introns > 0 ? total_length / mRNA->number_of_introns : 0;
    mRNA->largest_intron_length = largest_intron_length;
}