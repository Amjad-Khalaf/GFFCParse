#include <stdio.h>
#include <stdlib.h>
#include <string.h>         
#include "read_file.h"
#include "intron_operations.h"  
#include "exon_operations.h"
#include "genome_features.h"
#include "gene_count_hash.h"

char* get_gene_name_from_gff_tag(char *tag) 
{
    char *gene_name = split_line(split_line(tag, ';').array_of_elements[0], '=').array_of_elements[1];
    return gene_name;
}

void prepare_output_files(const char *mRNA_filename, const char *exon_filename, const char *intron_filename) 
{
    // exon output file
    FILE *exon_output_file = fopen(exon_filename, "a");
    if (exon_output_file == NULL) {
        perror("Error opening exon output file");
        return;
    }
    fprintf(exon_output_file, "gene\tgene_start\tgene_end\tgene_len\tmRNA\tmRNA_start\tmRNA_end\tmRNA_len\texon\texon_start\texon_end\texon_len\n");
    fclose(exon_output_file);

    // Open intron file here
    FILE *intron_output_file = fopen(intron_filename, "a");
    if (intron_output_file == NULL) {
        perror("Error opening intron output file");
        return;
    }
    fprintf(intron_output_file, "gene\tgene_start\tgene_end\tgene_len\tmRNA\tmRNA_start\tmRNA_end\tmRNA_len\tintron\tintron_start\tintron_end\tintron_len\n");
    fclose(intron_output_file);

    // Open mRNA file here
    FILE *mRNA_output_file = fopen(mRNA_filename, "a");
    if (mRNA_output_file == NULL) {
        perror("Error opening mRNA output file");
        return;
    }
    fprintf(mRNA_output_file, "gene\tnum_copies\tgene_start\tgene_end\tgene_len\tnum_mRNAs\tmRNA\tmRNA_start\tmRNA_end\tmRNA_len\tnum_exons\tsum_exon_len\tavg_exon_len\tmax_exon_len\tnum_introns\tsum_intron_len\tavg_intron_len\tmax_intron_len\n");

    // Close mRNA file here
    fclose(mRNA_output_file);
}

int count_genes_in_gff(const char *filename) 
{
    FILE *file;

    // Open the GFF file
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Cannot open input file");
        return -1;
    }

    size_t line_size = 256;
    char *line = (char *)malloc(line_size);
    if (line == NULL) {
        perror("Failed to allocate memory");
        fclose(file);
        return -1;
    }

    // Read GFF file line by line
    while (fgets(line, line_size, file)) {
        char *trimmed_line = remove_trailing_newline(line);

        // Skip header lines
        if (trimmed_line[0] == '#') {
            continue;
        }

        // Split the line into columns
        struct line parsed_line = split_line(trimmed_line, '\t');
        if (parsed_line.number_of_elements < 9) {
            perror("Invalid format! Not enough columns.");
            free(line);
            fclose(file);
            return -1;
        }

        // Process gene entries
        if (strcmp(parsed_line.array_of_elements[2], "gene") == 0) {
            char *gene_name = get_gene_name_from_gff_tag(parsed_line.array_of_elements[8]);
            if (gene_name != NULL) {
                // Increment gene count
                increment_gene_count(gene_name);
                free(gene_name);  // Free the gene name after use
            }
        }
    }

    free(line);
    fclose(file);

    return 0;  // Success
}

void parse_gff(const char *filename, const char *mRNA_filename, const char *exon_filename, const char *intron_filename) 
{    
    // Declare file pointer
    FILE *file;

    // Open file in read mode
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Cannot open input file, sorry!");
        return;
    }

    // Initial buffer size for the line
    size_t line_size = 256;

    // Allocate memory for the line buffer
    char *line = (char *)malloc(line_size);
    if (line == NULL) {
        perror("Failed to allocate memory, sorry!");
        fclose(file);
        return;
    }

    // Declare reused variables & structs
    struct gene current_gene;
    int number_of_genes = 0;
    char *mRNA_prefix = "mRNA_";

    // Initialize gene's mRNA array
    current_gene.array_of_mRNAs = NULL;
    current_gene.number_of_mRNAs = 0;

    // Read file line by line
    while (fgets(line, line_size, file)) {
        char *trimmed_line = remove_trailing_newline(line);          
        
        // Handle header lines
        if (trimmed_line[0] == '#') {
            continue;
        }

        // Separate line
        struct line parsed_line = split_line(trimmed_line, '\t');  
        if (parsed_line.number_of_elements < 9) {
            perror("Your file is not tab-separated, sorry!");
            free(line);
            fclose(file);
            return;
        }

        // Gene
        if (strcmp(parsed_line.array_of_elements[2], "gene") == 0) {
            // Resolve previous gene
            if (number_of_genes != 0) {


                // Print all mRNAs for the previous gene
                for (int j = 0; j < current_gene.number_of_mRNAs; j++) {
                    
                    // Sort exons before output
                    sort_exons(*current_gene.array_of_mRNAs[j]);

                    // Open exon file here
                    FILE *exon_output_file = fopen(exon_filename, "a");
                    if (exon_output_file == NULL) {
                        perror("Error opening exon output file");
                        return;
                    }

                    for (int i = 0; i < current_gene.array_of_mRNAs[j]->number_of_exons; i++) {
                        fprintf(exon_output_file, "%s\t%d\t%d\t%d\t%s_%s\t%d\t%d\t%d\t%s_%s_exon_%d\t%d\t%d\t%d\n", current_gene.name, current_gene.start, current_gene.end, current_gene.length, current_gene.name, current_gene.array_of_mRNAs[j]->name, current_gene.array_of_mRNAs[j]->start, current_gene.array_of_mRNAs[j]->end, current_gene.array_of_mRNAs[j]->length, current_gene.name, current_gene.array_of_mRNAs[j]->name, i + 1, current_gene.array_of_mRNAs[j]->array_of_exons[i]->start, current_gene.array_of_mRNAs[j]->array_of_exons[i]->end, current_gene.array_of_mRNAs[j]->array_of_exons[i]->length);
                    }
                    
                    // Close exon file here
                    fclose(exon_output_file);

                    // Add intron information
                    identify_introns(current_gene.array_of_mRNAs[j]);
                    
                    // Open intron file here
                    FILE *intron_output_file = fopen(intron_filename, "a");
                    if (intron_output_file == NULL) {
                        perror("Error opening intron output file");
                        return;
                    }

                    for (int i = 0; i < current_gene.array_of_mRNAs[j]->number_of_introns; i++) {
                        fprintf(intron_output_file, "%s\t%d\t%d\t%d\t%s_%s\t%d\t%d\t%d\t%s_%s_intron_%d\t%d\t%d\t%d\n", current_gene.name, current_gene.start, current_gene.end, current_gene.length, current_gene.name, current_gene.array_of_mRNAs[j]->name, current_gene.array_of_mRNAs[j]->start, current_gene.array_of_mRNAs[j]->end, current_gene.array_of_mRNAs[j]->length, current_gene.name, current_gene.array_of_mRNAs[j]->name, i + 1, current_gene.array_of_mRNAs[j]->array_of_introns[i]->start, current_gene.array_of_mRNAs[j]->array_of_introns[i]->end, current_gene.array_of_mRNAs[j]->array_of_introns[i]->length);                    
                    }

                    // Close intron file here
                    fclose(intron_output_file);
                    
                    // Add mRNA summary stats
                    total_exon_length(current_gene.array_of_mRNAs[j]);
                    total_intron_length(current_gene.array_of_mRNAs[j]);
                    calculate_exon_stats(current_gene.array_of_mRNAs[j]);
                    calculate_intron_stats(current_gene.array_of_mRNAs[j]);
                    
                    // Get gene count
                    int count = get_gene_count(current_gene.name);

                    // Open mRNA file here
                    FILE *mRNA_output_file = fopen(mRNA_filename, "a");
                    if (mRNA_output_file == NULL) {
                        perror("Error opening mRNA output file");
                        return;
                    }
                    
                    // Gene, Start, End, Length, #mRNAs, mRNA, Start, End, Length, #Exons, Sum_Exon_Len, Avg_Exon_Len, Max_Exon_Len, #Introns, Sum_Intron_Len,Avg_Intron_Len, Max_Intron_Len
                    fprintf(mRNA_output_file, "%s\t%d\t%d\t%d\t%d\t%d\t%s_%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", current_gene.name, count, current_gene.start, current_gene.end, current_gene.length, current_gene.number_of_mRNAs, current_gene.name, current_gene.array_of_mRNAs[j]->name, current_gene.array_of_mRNAs[j]->start, current_gene.array_of_mRNAs[j]->end, current_gene.array_of_mRNAs[j]->length, current_gene.array_of_mRNAs[j]->number_of_exons, current_gene.array_of_mRNAs[j]->total_exon_length, current_gene.array_of_mRNAs[j]->average_exon_length, current_gene.array_of_mRNAs[j]->largest_exon_length, current_gene.array_of_mRNAs[j]->number_of_introns, current_gene.array_of_mRNAs[j]->total_intron_length, current_gene.array_of_mRNAs[j]->average_intron_length, current_gene.array_of_mRNAs[j]->largest_intron_length);

                    // Close mRNA file here
                    fclose(mRNA_output_file);
                }
            }

            number_of_genes++;
            char *gene_name = get_gene_name_from_gff_tag(parsed_line.array_of_elements[8]); 
            current_gene.name = strdup(gene_name);  // strdup allocates memory and copies the string
            current_gene.start = atoi(parsed_line.array_of_elements[3]);
            current_gene.end = atoi(parsed_line.array_of_elements[4]);
            current_gene.length = atoi(parsed_line.array_of_elements[4]) - atoi(parsed_line.array_of_elements[3]);
            current_gene.number_of_mRNAs = 0;
            current_gene.array_of_mRNAs = NULL;  // Reset for new gene
        }

        // mRNA
        if (strcmp(parsed_line.array_of_elements[2], "mRNA") == 0) {
            current_gene.number_of_mRNAs++;  // Increase number of mRNAs
            struct mRNA **temp_mRNAs = realloc(current_gene.array_of_mRNAs, sizeof(struct mRNA*) * current_gene.number_of_mRNAs);
            if (temp_mRNAs == NULL) {
                perror("Failed to reallocate memory for mRNAs");
                free(line);
                fclose(file);
                return;
            }
            current_gene.array_of_mRNAs = temp_mRNAs;

            // Allocate new mRNA for each iteration
            struct mRNA *new_mRNA = (struct mRNA*)malloc(sizeof(struct mRNA));
            if (new_mRNA == NULL) {
                perror("Failed to allocate memory for new mRNA");
                free(line);
                fclose(file);
                return;
            }

            // Create the mRNA name
            char mRNA_name[20];  // Increase buffer size if needed
            sprintf(mRNA_name, "%s%d", mRNA_prefix, current_gene.number_of_mRNAs);

            new_mRNA->name = strdup(mRNA_name);  // Allocate memory for mRNA name
            new_mRNA->start = atoi(parsed_line.array_of_elements[3]);
            new_mRNA->end = atoi(parsed_line.array_of_elements[4]);
            new_mRNA->length = atoi(parsed_line.array_of_elements[4]) - atoi(parsed_line.array_of_elements[3]);
            new_mRNA->number_of_exons = 0;
            new_mRNA->array_of_exons = NULL;  // Initialize exons array

            // Add new mRNA to gene's array
            current_gene.array_of_mRNAs[current_gene.number_of_mRNAs - 1] = new_mRNA;
        }

        // Exon
        if (strcmp(parsed_line.array_of_elements[2], "exon") == 0) {
            if (current_gene.number_of_mRNAs == 0) {
                fprintf(stderr, "No mRNA available to add exon\n");
                continue;
            }
            struct mRNA *last_mRNA = current_gene.array_of_mRNAs[current_gene.number_of_mRNAs - 1];
            last_mRNA->number_of_exons++;
            struct exon **temp_exons = realloc(last_mRNA->array_of_exons, sizeof(struct exon*) * last_mRNA->number_of_exons);
            if (temp_exons == NULL) {
                perror("Failed to reallocate memory for exons");
                free(line);
                fclose(file);
                return;
            }
            last_mRNA->array_of_exons = temp_exons;

            // Allocate memory for new exon
            struct exon *new_exon = (struct exon*)malloc(sizeof(struct exon));
            if (new_exon == NULL) {
                perror("Failed to allocate memory for new exon");
                free(line);
                fclose(file);
                return;
            }
            new_exon->start = atoi(parsed_line.array_of_elements[3]);
            new_exon->end = atoi(parsed_line.array_of_elements[4]);
            new_exon->length = atoi(parsed_line.array_of_elements[4]) - atoi(parsed_line.array_of_elements[3]);

            // Add new exon to mRNA's exon array
            last_mRNA->array_of_exons[last_mRNA->number_of_exons - 1] = new_exon;
        }
    }

    // Handle last gene
    for (int j = 0; j < current_gene.number_of_mRNAs; j++) 
    {
        // Open exon file here
        FILE *exon_output_file = fopen(exon_filename, "a");
        if (exon_output_file == NULL) {
            perror("Error opening exon output file");
            return;
        }
        for (int i = 0; i < current_gene.array_of_mRNAs[j]->number_of_exons; i++) 
        {
            fprintf(exon_output_file, "%s\t%d\t%d\t%d\t%s_%s\t%d\t%d\t%d\t%s_%s_exon_%d\t%d\t%d\t%d\n", current_gene.name, current_gene.start, current_gene.end, current_gene.length, current_gene.name, current_gene.array_of_mRNAs[j]->name, current_gene.array_of_mRNAs[j]->start, current_gene.array_of_mRNAs[j]->end, current_gene.array_of_mRNAs[j]->length, current_gene.name, current_gene.array_of_mRNAs[j]->name, i + 1, current_gene.array_of_mRNAs[j]->array_of_exons[i]->start, current_gene.array_of_mRNAs[j]->array_of_exons[i]->end, current_gene.array_of_mRNAs[j]->array_of_exons[i]->length);
        }
        fclose(exon_output_file);

        // Add intron information
        identify_introns(current_gene.array_of_mRNAs[j]);
                    
        // Open intron file here
        FILE *intron_output_file = fopen(intron_filename, "a");
        if (intron_output_file == NULL) {
            perror("Error opening intron output file");
            return;
        }

        for (int i = 0; i < current_gene.array_of_mRNAs[j]->number_of_introns; i++) {
            fprintf(intron_output_file, "%s\t%d\t%d\t%d\t%s_%s\t%d\t%d\t%d\t%s_%s_intron_%d\t%d\t%d\t%d\n", current_gene.name, current_gene.start, current_gene.end, current_gene.length, current_gene.name, current_gene.array_of_mRNAs[j]->name, current_gene.array_of_mRNAs[j]->start, current_gene.array_of_mRNAs[j]->end, current_gene.array_of_mRNAs[j]->length, current_gene.name, current_gene.array_of_mRNAs[j]->name, i + 1, current_gene.array_of_mRNAs[j]->array_of_introns[i]->start, current_gene.array_of_mRNAs[j]->array_of_introns[i]->end, current_gene.array_of_mRNAs[j]->array_of_introns[i]->length);                       
        }
        fclose(intron_output_file);

        // Add mRNA summary stats
        total_exon_length(current_gene.array_of_mRNAs[j]);
        total_intron_length(current_gene.array_of_mRNAs[j]);
        calculate_exon_stats(current_gene.array_of_mRNAs[j]);
        calculate_intron_stats(current_gene.array_of_mRNAs[j]);

        // Get gene count
        int count = get_gene_count(current_gene.name);

        // Open mRNA file here
        FILE *mRNA_output_file = fopen(mRNA_filename, "a");
        if (mRNA_output_file == NULL) {
            perror("Error opening mRNA output file");
            return;
        }
                    
        // Gene, Start, End, Length, #mRNAs, mRNA, Start, End, Length, #Exons, Sum_Exon_Len, Avg_Exon_Len, Max_Exon_Len, #Introns, Sum_Intron_Len,Avg_Intron_Len, Max_Intron_Len
        fprintf(mRNA_output_file, "%s\t%d\t%d\t%d\t%d\t%d\t%s_%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", current_gene.name, count, current_gene.start, current_gene.end, current_gene.length, current_gene.number_of_mRNAs, current_gene.name, current_gene.array_of_mRNAs[j]->name, current_gene.array_of_mRNAs[j]->start, current_gene.array_of_mRNAs[j]->end, current_gene.array_of_mRNAs[j]->length, current_gene.array_of_mRNAs[j]->number_of_exons, current_gene.array_of_mRNAs[j]->total_exon_length, current_gene.array_of_mRNAs[j]->average_exon_length, current_gene.array_of_mRNAs[j]->largest_exon_length, current_gene.array_of_mRNAs[j]->number_of_introns, current_gene.array_of_mRNAs[j]->total_intron_length, current_gene.array_of_mRNAs[j]->average_intron_length, current_gene.array_of_mRNAs[j]->largest_intron_length);        
        
        // Close mRNA file here
        fclose(mRNA_output_file);
    }

    // Free allocated memory
    free(line);
    fclose(file);

    // Free memory for mRNAs and exons
    for (int i = 0; i < current_gene.number_of_mRNAs; i++) {
        struct mRNA *m = current_gene.array_of_mRNAs[i];
        for (int k = 0; k < m->number_of_exons; k++) {
            free(m->array_of_exons[k]);  // Free each exon
        }
        free(m->array_of_exons);  // Free the array of exon pointers
        free(m->name);  // Free mRNA name
        free(m);  // Free mRNA struct
    }
    free(current_gene.array_of_mRNAs);  // Free the array of mRNA pointers
}

void summarise_gff(const char *mRNA_filename, const char *summary_filename) 
{
    // Declare file pointer
    FILE *file;

    // Open file in read mode
    file = fopen(mRNA_filename, "r");
    if (file == NULL) {
        perror("Cannot open input file, sorry!");
        return;
    }

    // Initial buffer size for the line
    size_t line_size = 256;

    // Allocate memory for the line buffer
    char *line = (char *)malloc(line_size);
    if (line == NULL) {
        perror("Failed to allocate memory, sorry!");
        fclose(file);
        return;
    }

    // Declare summary vars
    int number_of_single_isoform_genes = 0;
    int number_of_multi_isoform_genes = 0;

    int number_of_single_copy_isoforms = 0;
    int len_of_single_copy_isoforms = 0;
    int number_of_introns_single_copy_isoforms = 0;
    int len_of_introns_single_copy_isoforms = 0;
    int number_of_exons_single_copy_isoforms = 0;
    int len_of_exons_single_copy_isoforms = 0;  
    int number_of_mono_exonic_single_copy_isoforms = 0;
    int len_of_mono_exonic_single_copy_isoforms = 0;
    int number_of_multi_exonic_single_copy_isoforms = 0;
    int len_of_multi_exonic_single_copy_isoforms = 0;

    int number_of_multi_copy_isoforms = 0;
    int len_of_multi_copy_isoforms = 0;
    int number_of_introns_multi_copy_isoforms = 0;
    int len_of_introns_multi_copy_isoforms = 0;
    int number_of_exons_multi_copy_isoforms = 0;
    int len_of_exons_multi_copy_isoforms = 0;
    int number_of_mono_exonic_multi_copy_isoforms = 0;
    int len_of_mono_exonic_multi_copy_isoforms = 0;
    int number_of_multi_exonic_multi_copy_isoforms = 0;
    int len_of_multi_exonic_multi_copy_isoforms = 0;



    // Read file line by line
    while (fgets(line, line_size, file)) {
        char *trimmed_line = remove_trailing_newline(line);          

        // Separate line
        struct line parsed_line = split_line(trimmed_line, '\t');  
        if (parsed_line.number_of_elements < 9) {
            perror("Your file is not tab-separated, sorry!");
            free(line);
            fclose(file);
            return;
        }

        // Handle header line
        if (strcmp(parsed_line.array_of_elements[1], "num_gene_copies") == 0) {
            continue;
        }

        // Single Isoform Genes
        if (atoi(parsed_line.array_of_elements[5]) == 1)
        {
            number_of_single_isoform_genes++;
        }

        // Multi Isoform Genes
        if (atoi(parsed_line.array_of_elements[5]) > 1) 
        {
            number_of_multi_isoform_genes++;
        }
        // Single Copy Isoforms
        if (atoi(parsed_line.array_of_elements[1]) == 1)
        {
            number_of_single_copy_isoforms++;
            len_of_single_copy_isoforms = len_of_single_copy_isoforms + atoi(parsed_line.array_of_elements[9]);
            number_of_introns_single_copy_isoforms = number_of_introns_single_copy_isoforms + atoi(parsed_line.array_of_elements[14]);
            len_of_introns_single_copy_isoforms = len_of_introns_single_copy_isoforms + atoi(parsed_line.array_of_elements[15]);
            
            number_of_exons_single_copy_isoforms = number_of_exons_single_copy_isoforms + atoi(parsed_line.array_of_elements[10]);
            len_of_exons_single_copy_isoforms = len_of_exons_single_copy_isoforms + atoi(parsed_line.array_of_elements[11]);
            
            if (atoi(parsed_line.array_of_elements[10]) == 1) 
            {
                number_of_mono_exonic_single_copy_isoforms++;
                len_of_mono_exonic_single_copy_isoforms = len_of_mono_exonic_single_copy_isoforms + atoi(parsed_line.array_of_elements[9]);             
            } 
            else if (atoi(parsed_line.array_of_elements[10]) > 1) 
            {
                number_of_multi_exonic_single_copy_isoforms ++;
                len_of_multi_exonic_single_copy_isoforms = len_of_multi_exonic_single_copy_isoforms + atoi(parsed_line.array_of_elements[9]);
            }

        // Multi Copy Isoforms
        } else if (atoi(parsed_line.array_of_elements[1]) > 1) 
        {
            number_of_multi_copy_isoforms++;
            len_of_multi_copy_isoforms = len_of_multi_copy_isoforms + atoi(parsed_line.array_of_elements[9]);
            number_of_introns_multi_copy_isoforms = number_of_introns_multi_copy_isoforms + atoi(parsed_line.array_of_elements[14]);
            len_of_introns_multi_copy_isoforms = len_of_introns_multi_copy_isoforms + atoi(parsed_line.array_of_elements[15]);
            
            number_of_exons_multi_copy_isoforms = number_of_exons_multi_copy_isoforms + atoi(parsed_line.array_of_elements[10]);
            len_of_exons_multi_copy_isoforms = len_of_exons_multi_copy_isoforms + atoi(parsed_line.array_of_elements[11]);
            
            if (atoi(parsed_line.array_of_elements[10]) == 1) 
            {
                number_of_mono_exonic_multi_copy_isoforms++;
                len_of_mono_exonic_multi_copy_isoforms = len_of_mono_exonic_multi_copy_isoforms + atoi(parsed_line.array_of_elements[9]);             
            } 
            else if (atoi(parsed_line.array_of_elements[10]) > 1) 
            {
                number_of_multi_exonic_multi_copy_isoforms ++;
                len_of_multi_exonic_multi_copy_isoforms = len_of_multi_exonic_multi_copy_isoforms + atoi(parsed_line.array_of_elements[9]);
            }
        }
    }

    int total_gene_number = get_total_gene_count();

    // Open summary output file here
    FILE *summary_output_file = fopen(summary_filename, "a");
        if (summary_output_file == NULL) {
            perror("Error opening intron output file");
            return;
        }

    fprintf(summary_output_file, "GFF SUMMARY\n");
    fprintf(summary_output_file, "=================================================================\n");
    fprintf(summary_output_file, "%-50s %10s\n", "Description", "Count / Length / Percentage");
    fprintf(summary_output_file, "-----------------------------------------------------------------\n");

    // Overall summary
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of genes:", total_gene_number);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of all isoforms:", number_of_multi_copy_isoforms + number_of_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10.2f%%\n", "Percentage of single-isoform genes:", (number_of_single_isoform_genes / (float)(number_of_multi_copy_isoforms + number_of_single_copy_isoforms)) * 100);
    fprintf(summary_output_file, "%-50s %10.2f%%\n", "Percentage of multi-isoform genes:", (number_of_multi_isoform_genes / (float)(number_of_multi_copy_isoforms + number_of_single_copy_isoforms)) * 100);
    fprintf(summary_output_file, "%-50s %10.2f%%\n", "Percentage of mono-exonic isoforms (intronless):", ((number_of_mono_exonic_multi_copy_isoforms + number_of_mono_exonic_single_copy_isoforms) / (float)(number_of_multi_copy_isoforms + number_of_single_copy_isoforms)) * 100);
    fprintf(summary_output_file, "%-50s %10.2f%%\n", "Percentage of multi-exonic isoforms (with introns):", ((number_of_multi_exonic_multi_copy_isoforms + number_of_multi_exonic_single_copy_isoforms) / (float)(number_of_multi_copy_isoforms + number_of_single_copy_isoforms)) * 100);

    fprintf(summary_output_file, "-----------------------------------------------------------------\n");

    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of all isoforms:", len_of_multi_copy_isoforms + len_of_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of exons detected:", number_of_exons_multi_copy_isoforms + number_of_exons_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of exons detected:", len_of_exons_multi_copy_isoforms + len_of_exons_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of introns detected:", number_of_introns_multi_copy_isoforms + number_of_introns_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of introns detected:", len_of_introns_multi_copy_isoforms + len_of_introns_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of all mono-exonic isoforms:", number_of_mono_exonic_multi_copy_isoforms + number_of_mono_exonic_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of all mono-exonic isoforms:", len_of_mono_exonic_multi_copy_isoforms + len_of_mono_exonic_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of all multi-exonic isoforms:", number_of_multi_exonic_multi_copy_isoforms + number_of_multi_exonic_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of all multi-exonic isoforms:", len_of_multi_exonic_multi_copy_isoforms + len_of_multi_exonic_single_copy_isoforms);

    fprintf(summary_output_file, "=================================================================\n");

    // Single-copy summary
    fprintf(summary_output_file, "Single-Copy Summary\n");
    fprintf(summary_output_file, "-----------------------------------------------------------------\n");
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of single-copy isoforms:", number_of_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of single-copy isoforms:", len_of_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of exons in single-copy isoforms:", number_of_exons_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of exons in single-copy isoforms:", len_of_exons_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of introns in single-copy isoforms:", number_of_introns_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of introns in single-copy isoforms:", len_of_introns_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of mono-exonic single-copy isoforms:", number_of_mono_exonic_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of mono-exonic single-copy isoforms:", len_of_mono_exonic_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of multi-exonic single-copy isoforms:", number_of_multi_exonic_single_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of multi-exonic single-copy isoforms:", len_of_multi_exonic_single_copy_isoforms);

    fprintf(summary_output_file, "=================================================================\n");

    // Multi-copy summary
    fprintf(summary_output_file, "Multi-Copy Summary\n");
    fprintf(summary_output_file, "-----------------------------------------------------------------\n");
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of multi-copy isoforms:", number_of_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of multi-copy isoforms:", len_of_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of exons in multi-copy isoforms:", number_of_exons_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of exons in multi-copy isoforms:", len_of_exons_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of introns in multi-copy isoforms:", number_of_introns_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of introns in multi-copy isoforms:", len_of_introns_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of mono-exonic multi-copy isoforms:", number_of_mono_exonic_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of mono-exonic multi-copy isoforms:", len_of_mono_exonic_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10d\n", "Total number of multi-exonic multi-copy isoforms:", number_of_multi_exonic_multi_copy_isoforms);
    fprintf(summary_output_file, "%-50s %10dbp\n", "Total length of multi-exonic multi-copy isoforms:", len_of_multi_exonic_multi_copy_isoforms);

    fprintf(summary_output_file, "=================================================================\n");

    // Free allocated memory
    free(line);
    fclose(file);
}
