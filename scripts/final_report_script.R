##################################################
## Project: EJP RD Istanbul Workshop
## Purpose: Create final report
## Date: August 2022
## Author: Berk Gurdamar
##################################################

args = commandArgs(trailingOnly=TRUE)

intervar_out <- data.table::fread(args[1])

intervar_out$acmg_classification <- sapply(intervar_out$`InterVar: InterVar and Evidence`, function(x) strsplit(strsplit(x, "InterVar: ")[[1]][2], " PVS")[[1]][1])

intervar_out$idx <- paste0(intervar_out$`#Chr`, "_",
                           intervar_out$Start, "_",
                           intervar_out$End, "_",
                           intervar_out$Ref, "_",
                           intervar_out$Alt)

intervar_out <- intervar_out[,35:36]
intervar_out$idx <- paste0("chr", intervar_out$idx)


gemini_out <- data.table::fread(args[2])

gemini_out$start <- gemini_out$start + 1

gemini_out$idx <- paste0(gemini_out$chrom, "_",
                           gemini_out$start, "_",
                           gemini_out$end, "_",
                           gemini_out$ref, "_",
                           gemini_out$alt)

comb_df <- merge(x=gemini_out, y=intervar_out, by="idx", all.x=TRUE)

predictions <- data.table::fread(args[3])

predictions$idx <- paste0(predictions$Chr, "_",
                         predictions$Start, "_",
                         predictions$End, "_",
                         predictions$Ref, "_",
                         predictions$Alt)

predictions <- predictions[,c(14:117,150)]

comb_df <- merge(x=comb_df, y=predictions, by="idx", all.x=TRUE)
comb_df <- comb_df[,-1]

data.table::fwrite(comb_df, args[4], sep = "\t")
