// exon_operations.h

#ifndef EXON_OPERATIONS_H
#define EXON_OPERATIONS_H

#include "genome_features.h"  // For struct exon and struct mRNA

void sort_exons(struct mRNA mRNA);
void total_exon_length(struct mRNA *mRNA);
void calculate_exon_stats(struct mRNA *mRNA);

#endif  // EXON_OPERATIONS_H