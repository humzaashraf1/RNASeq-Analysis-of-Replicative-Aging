from pytximport import tximport
from pytximport.utils import create_transcript_gene_map
import os
import pandas as pd

quant_files = []
salmon_path = '/path/salmon_output'
for folder in os.listdir(salmon_path):
    quant_path = salmon_path + '/' + folder + '/' + 'quant.sf'
    quant_files.append(quant_path)

transcript_gene_map = create_transcript_gene_map(species="human")

results = tximport(
    file_paths=quant_files,
    data_type="salmon",
    transcript_gene_map=transcript_gene_map
)

expression_df = pd.DataFrame(
    results.X.toarray() if hasattr(results.X, "toarray") else results.X,
    index=results.obs.index, 
    columns=results.var.index
)
combined_df = pd.concat([results.obs, expression_df], axis=1)