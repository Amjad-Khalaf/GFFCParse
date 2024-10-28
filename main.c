#include <stdio.h>
#include "gene_count_hash.h"     
#include "genome_features.h"      
#include "exon_operations.h"      
#include "intron_operations.h"
#include "gff_operations.h"    

int main(int argc, char *argv[]) {
    
    printf("\n");
    
    printf("░█▀▀░█▀▀░█▀▀░█▀▀░█▀█░█▀█░█▀▄░█▀▀░█▀▀\n");
    printf("░█░█░█▀▀░█▀▀░█░░░█▀▀░█▀█░█▀▄░▀▀█░█▀▀\n");
    printf("░▀▀▀░▀░░░▀░░░▀▀▀░▀░░░▀░▀░▀░▀░▀▀▀░▀▀▀\n");

    printf("\n");
    
    printf("          ▒▒     ▓▓▓▓▓              \n");
    printf("          ▓▓▓▒▒  ▒▓▓▓▓▓             \n");
    printf("          ▓▓▓▒▓▓▓▓▓▓▓▓▓             \n");
    printf("           ▓▓▓▓▓▓▓▓▓▓▓▓▓            \n");
    printf("           ▓▒▓▓▓▓▓▓▒▒▓▓▓▒▒▒         \n");
    printf("          █▓▓▓▓▓▓█▓▒▒▒▓▓▓▒▒▒        \n");
    printf("         ▓▓▓▒▒▒▓▓▓▒▒▓▒▒▒▒▒▒▒▒       \n");
    printf("         ▒▒▒▒▓▓▒▒▓▓▓▓▓▒▒▒▒▒▒▒▒      \n");
    printf("          ▒▒▒▒▒▒▓▒▒▒▒▒▒▒▒▒▒▒▒▒      \n");
    printf("           ▒▒▒▒▒▒▒▒▓▓▒▒▒▒▒▒▒▒▒      \n");
    printf("            ▒▒▒▒▒▒▒▒▓▓▒▒▒▒▒▒▒       \n");
    printf("            ▓▒▒▓▓▓▓▒▒▒▒ ██▓▓        \n");
    printf("            █▓▓▓ █▓▓▓               \n");
    
    printf("\n");
    
    printf(" GFFCParse is a C programme to parse\n"); 
    printf(" a GFF file and: 1. summarise it; 2.\n");
    printf("  extract intron, exon & mRNA sizes \n");
    printf(" locations, and counts; and 3. make \n");
    printf("      png & svg summary plots.      \n");

    printf("____________________________________\n");
    
    printf("         version 1.0.0\n");
    printf("\n");

                                                                                        

    // Check if the correct number of arguments is passed
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <gff_file> <output_prefix>\n", argv[0]);
        return 1;
    }

    // Assign the command line arguments to variables
    const char *gff_file = argv[1];
    const char *prefix = argv[2];

    #define MAX_FILENAME_LEN 256

    // Buffer to hold the full output filename
    char mRNA_output_filename[MAX_FILENAME_LEN];
    char exon_output_filename[MAX_FILENAME_LEN];
    char intron_output_filename[MAX_FILENAME_LEN];
    char summary_output_filename[MAX_FILENAME_LEN];

    // Construct the output filenames
    if (snprintf(mRNA_output_filename, MAX_FILENAME_LEN, "%s%s", prefix, "_mRNA.tsv") >= MAX_FILENAME_LEN) 
    {
        fprintf(stderr, "Output prefix chosen is too long\n");
        return 1;
    }
    if (snprintf(exon_output_filename, MAX_FILENAME_LEN, "%s%s", prefix, "_exon.tsv") >= MAX_FILENAME_LEN) 
    {
        fprintf(stderr, "Output prefix chosen is too long\n");
        return 1;
    }
    if (snprintf(intron_output_filename, MAX_FILENAME_LEN, "%s%s", prefix, "_intron.tsv") >= MAX_FILENAME_LEN) 
    {
        fprintf(stderr, "Output prefix chosen is too long\n");
        return 1;
    }
    if (snprintf(summary_output_filename, MAX_FILENAME_LEN, "%s%s", prefix, "_summary.tbl") >= MAX_FILENAME_LEN) 
    {
        fprintf(stderr, "Output prefix chosen is too long\n");
        return 1;
    }

    // Process the GFF file and generate output
    if (count_genes_in_gff(gff_file) != 0) {
        perror("Something is wrong with your gene names/gff tag, check again!");
        return -1;
    }

    prepare_output_files(mRNA_output_filename, exon_output_filename, intron_output_filename);
    parse_gff(gff_file, mRNA_output_filename, exon_output_filename, intron_output_filename);
    summarise_gff(mRNA_output_filename, summary_output_filename);

    char command[256];
    snprintf(command, sizeof(command), "Rscript plot.r %s", mRNA_output_filename);
    int result = system(command);

    if (result != 0) {
        printf("Error: plot.r failed to generate plots.\n");
    } else {
        printf("Enjoy your results :)\n");
    }

    return 0;
}
