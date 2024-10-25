// genome_features.h

#ifndef GENOME_FEATURES_H
#define GENOME_FEATURES_H

// Exon struct
struct exon {
    int start;
    int end;
    int length;
};

// Intron struct
struct intron {
    int start;
    int end;
    int length;
};

// mRNA struct
struct mRNA {
    char *name; // Generic names like mRNA_1, _2, _3 will be used
    int start;
    int end;
    int length;
    int number_of_exons;
    struct exon **array_of_exons;
    int number_of_introns;
    struct intron **array_of_introns;
    int total_exon_length;
    int total_intron_length;
    int average_exon_length;
    int largest_exon_length;
    int average_intron_length;
    int largest_intron_length;
};

// Gene struct
struct gene {
    char *name;
    int start;
    int end;
    int length;
    int number_of_mRNAs;
    struct mRNA **array_of_mRNAs;
};

#endif // GENOME_FEATURES_H
