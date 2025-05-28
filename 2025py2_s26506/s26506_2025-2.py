#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = "BioScriptEx10"

    def search_taxid(self, taxid):
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        print(f"Organism: {records[0]['ScientificName']} (TaxID: {taxid})")

        handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y")
        result = Entrez.read(handle)
        self.count = int(result["Count"])
        self.webenv = result["WebEnv"]
        self.query_key = result["QueryKey"]
        print(f"Found {self.count} records.")
        return self.count

    def fetch_filtered_records(self, min_len, max_len, max_fetch=100):
        batch_size = min(max_fetch, 500)
        handle = Entrez.efetch(
            db="nucleotide", rettype="gb", retmode="text",
            retmax=batch_size, webenv=self.webenv, query_key=self.query_key
        )
        records = list(SeqIO.parse(handle, "genbank"))
        filtered = [r for r in records if min_len <= len(r.seq) <= max_len]
        print(f"Filtered {len(filtered)} records in range [{min_len}, {max_len}].")
        return filtered

def save_csv(records, filename="output.csv"):
    data = [{
        "Accession": r.id,
        "Length": len(r.seq),
        "Description": r.description
    } for r in records]
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"Saved CSV to {filename}")
    return df

def plot_lengths(df, filename="plot.png"):
    df_sorted = df.sort_values(by="Length", ascending=False)
    plt.figure(figsize=(10, 5))
    plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o')
    plt.xticks(rotation=90)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title("GenBank Record Lengths")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot to {filename}")

def main():
    email = input("Enter email: ")
    api_key = input("Enter api_key: ")
    taxid = input("Enter TaxID: ")
    min_len = int(input("Minimum sequence length: "))
    max_len = int(input("Maximum sequence length: "))

    retriever = NCBIRetriever(email, api_key)
    if retriever.search_taxid(taxid):
        records = retriever.fetch_filtered_records(min_len, max_len)
        if records:
            df = save_csv(records)
            plot_lengths(df)
        else:
            print("No records matched the length criteria.")

if __name__ == "__main__":
    main()
