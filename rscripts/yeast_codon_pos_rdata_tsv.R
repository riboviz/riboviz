# Convert yeast_codon_pos_i200.RData RData to TSV file
# 
# yeast_codon_pos_i200.RData holds aa list of 61 matrixes.
# Each matrix is named after a codon ("TTT", "TTC", ... ,"GGG").
# Each matrix has two columns - gene and codon position - the position
# of that codon in the gene's sequence.
# Each matrix has N rows, one for each occurrence of the codon in the
# gene's sequence.
#
# The script operates as follows:
#
# 1. Load yeast_codon_pos_i200.RData.
# 2. Create a matrix with columns Gene, Codon, Pos
# 3. For each codon, add the rows of its matrix to the above matrix,
#    using the codon name for the values of the Codon column.
# 4. Convert matrix to data frame.
# 5. Convert factors to characters.
# 6. Convert Pos column to numberic.
# 7. Sort by Gene then by Pos.
# 8. Save data frame as TSV file.

suppressMessages(library(optparse, quietly = TRUE))
option_list <- list(
  make_option("--input-file",
              type = "character",
	      default = "data/yeast_codon_pos_i200.RData",
              help = "Yeast codon positions RData file"
  ),
  make_option("--output-file",
              type = "character",
	      default = "data/yeast_codon_pos_i200.tsv",
              help = "Yeast codon positions TSV file"
  )
)
opt <- optparse::parse_args(OptionParser(option_list = option_list),
                            convert_hyphens_to_underscores = TRUE)
attach(opt)

print(paste0("Loading ", input_file))
load(input_file)
print(paste0("class(codon_pos): ", class(codon_pos)))
print(paste0("typeof(codon_pos): ", typeof(codon_pos)))
print("Information about codon_pos entry for TTT:")
ttt = codon_pos[['TTT']]
print(paste0("class(ttt): ", class(ttt)))
print(paste0("typeof(ttt): ", typeof(ttt)))
head(ttt)

print("Creating a matrix with columns Gene, Codon, Pos")
all_genes <- matrix(, nrow = 0, ncol = 3)
colnames(all_genes) <- c("Gene", "Codon", "Pos")
print(all_genes)

print("For each codon, add the rows of its matrix to the above matrix, using the codon name for the values of the Codon column")
num_rows <- 0
codons <- attr(codon_pos, "name")
for (codon in codons) {
    print(codon)
    genes <- codon_pos[[codon]]
    num_rows <- num_rows + nrow(genes)
    genes_codon <- cbind(genes, codon)
    genes_codon <- genes_codon[, c(1, 3, 2)]
    all_genes <- rbind(all_genes, genes_codon)
}
print(paste0("Expected number of rows: ", num_rows))
print(paste0("Actual number of rows: ", nrow(all_genes)))
head(all_genes)

print("Converting matrix to data frame")
# Matrixes must have values of the same type, dataframes allow different types
all_genes_df <- as.data.frame(all_genes)
sapply(all_genes_df, typeof)
sapply(all_genes_df, class)

print("Converting factors to characters")
all_genes_df <- data.frame(lapply(all_genes_df, as.character), stringsAsFactors=FALSE)
sapply(all_genes_df, typeof)
sapply(all_genes_df, class)

print("Converting Pos to numeric")
all_genes_df <- transform(all_genes_df, Pos = as.numeric(Pos))
sapply(all_genes_df, typeof)
sapply(all_genes_df, class)
head(all_genes_df)
nrow(all_genes_df)

print("Sorting by Gene then by Pos")
all_genes_df <- all_genes_df[with(all_genes_df, order(Gene, Pos)), ]
# all_genes_df <- all_genes_df[order(all_genes_df$Gene, all_genes_df$Pos),]
head(all_genes_df)
nrow(all_genes_df)

print(paste0("Save data frame as TSV file to ", output_file))
write.table(all_genes_df, output_file, row.names=FALSE, sep="\t", quote=FALSE)
