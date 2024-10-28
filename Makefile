# Define the compiler
CC = gcc

# Define the compiler flags
CFLAGS = -Wall -g

# Define the target executable
TARGET = GFFCParse

# Define the source files
SRCS = main.c gene_count_hash.c exon_operations.c intron_operations.c gff_operations.c read_file.c

# Define the object files (replace .c with .o)
OBJS = $(SRCS:.c=.o)

# Define the rule for building the target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# Rule for compiling source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJS) $(TARGET)
