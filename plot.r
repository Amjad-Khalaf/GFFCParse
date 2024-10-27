library(ggplot2)
library(dplyr)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("R script did not receive input file properly, check how it is being passed to main.c")
}

file_path <- args[1]

data <- read.delim(file_path, header = TRUE)

data <- data %>%
  mutate(
    Copy_Number_Category = ifelse(`num_copies` == 1, "Single-copy Transcripts", "Multi-copy Transcripts"),
    Exon_Category = ifelse(`num_exons` == 1, "Mono-exonic", "Multi-exonic"))

create_and_save_plot <- function(plot, filename_base, width = 8, height = 6) {
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height, dpi = 300)
  ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height)
}

# Plot 1: Exon Count
exon_count <- data %>%
  group_by(`num_exons`, Copy_Number_Category) %>%
  summarise(Transcript_Count = n(), .groups = 'drop')

ten_percent_line <- sum(exon_count$Transcript_Count) * 0.1

exon_plot <- ggplot(exon_count, aes(x = `num_exons`, y = Transcript_Count, fill = Copy_Number_Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Single-copy Transcripts" = "#00AFBB", "Multi-copy Transcripts" = "#E7B800")) +
  geom_hline(yintercept = ten_percent_line, linetype = "dashed", color = "red") +
  annotate("text", x = max(exon_count$num_exons), y = ten_percent_line, label = "10% of Transcripts", color = "red", vjust = -0.5, hjust = 1) +
  labs(x = "\nNumber of Exons", y = "Number of Transcripts\n") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) +
  scale_x_continuous(breaks = seq(0, max(exon_count$num_exons), by = 1))

create_and_save_plot(exon_plot, "exon_count")

# Plot 2: Intron Count
intron_count <- data %>%
  group_by(`num_introns`, Copy_Number_Category) %>%
  summarise(Transcript_Count = n(), .groups = 'drop')

ten_percent_line_intron <- sum(intron_count$Transcript_Count) * 0.1

intron_plot <- ggplot(intron_count, aes(x = `num_introns`, y = Transcript_Count, fill = Copy_Number_Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Single-copy Transcripts" = "#00AFBB", "Multi-copy Transcripts" = "#E7B800")) +
  geom_hline(yintercept = ten_percent_line_intron, linetype = "dashed", color = "red") +
  annotate("text", x = max(intron_count$num_introns), y = ten_percent_line_intron, label = "10% of Transcripts", color = "red", vjust = -0.5, hjust = 1) +
  labs(x = "\nNumber of Introns", y = "Number of Transcripts\n") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) +
  scale_x_continuous(breaks = seq(0, max(intron_count$num_introns), by = 1))

create_and_save_plot(intron_plot, "intron_count")

# Plot 3: mRNA Length Density
length_density_plot <- ggplot(data, aes(x = mRNA_len, color = Exon_Category, fill = Exon_Category)) +
  geom_density(linewidth = 0.5) +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("#CBCE91FF", "#EA738DFF")) +
  scale_y_continuous(labels = scales::scientific_format()) +
  scale_x_continuous(labels = comma) +
  labs(x = "\nTranscript Length / bp", y = "Density\n") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))

create_and_save_plot(length_density_plot, "mRNA_length")

# Plot 4: Number of mRNAs
mRNA_count_plot <- ggplot(data, aes(x = num_mRNAs)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#FAEBEFFF") +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
  labs(x = "\nNumber of Transcripts", y = "Number of Genes\n") +
  theme_bw() +
  theme(legend.position = "none")

create_and_save_plot(mRNA_count_plot, "mRNA_count")
