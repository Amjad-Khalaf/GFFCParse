// gff_operations.h

#ifndef GFF_OPERATIONS_H
#define GFF_OPERATIONS_H
         
#include "read_file.h"
#include "intron_operations.h"  
#include "exon_operations.h"
#include "genome_features.h"
#include "gene_count_hash.h"

char* get_gene_name_from_gff_tag(char *tag);
void prepare_output_files(const char *mRNA_filename, const char *exon_filename, const char *intron_filename);
int count_genes_in_gff(const char *filename);
void parse_gff(const char *filename, const char *mRNA_filename, const char *exon_filename, const char *intron_filename);
void summarise_gff(const char *mRNA_filename, const char *summary_filename);

#endif  // GFF_OPERATIONS_H
