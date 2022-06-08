import pandas as pd 
from pathlib import Path

HERE = Path(__file__).parent.resolve()

df  = pd.read_csv(HERE.parent.joinpath("data/wikidata_GO_terms.tsv"), sep="\t")

df = df[["sitelink", "itemLabel", "gene_symbol"]]

df = df.drop_duplicates()

full_gmt = ""
for i in df.groupby("itemLabel"):
  itemLabel = i[0]
  sitelink = i[1]["sitelink"].values[0]
  genes = i[1]["gene_symbol"].values
  full_gmt += f"""{itemLabel}	{sitelink}	{"	".join(genes)}
"""

HERE.parent.joinpath("datasets/wikigoa_gene_symbol.gmt").write_text(full_gmt)