// intron_operations.h

#ifndef INTRON_OPERATIONS_H
#define INTRON_OPERATIONS_H

#include "genome_features.h"  // For struct intron and struct mRNA

void identify_introns(struct mRNA *mRNA);
void total_intron_length(struct mRNA *mRNA);
void calculate_intron_stats(struct mRNA *mRNA);

#endif  // INTRON_OPERATIONS_H