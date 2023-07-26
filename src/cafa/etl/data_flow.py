import joblib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

from .extraction import (
    dataframe_from_tsv,
    gene_ontology_reader,
    protein_sequence_reader,
)
from .transformations import (
    GO_counts_by_protein,
    add_relationship_count_to_go,
    add_subontologies_to_go,
    aspect_distributions,
    generate_amino_count,
    generate_subontologies,
    protein_sequence_length,
)

TRAIN_TERMS_FILE = "../data/cafa-5-protein-function-prediction/Train/train_terms.tsv"
TRAIN_TAX_FILE = "../data/cafa-5-protein-function-prediction/Train/train_taxonomy.tsv"
PROTEIN_SEQUENCE_FILE = "../data/cafa-5-protein-function-prediction/Train/train_sequences.fasta"
GENE_ONTOLOGY_FILE = "../data/cafa-5-protein-function-prediction/Train/go-basic.obo"


def store_data(data, filepath):
    joblib.dump(data, filepath)


def data_flow(go_filepath, protein_filepath):
    # Extraction
    train_terms = dataframe_from_tsv(TRAIN_TERMS_FILE)
    train_tax = dataframe_from_tsv(TRAIN_TAX_FILE)
    train_sequences = protein_sequence_reader(PROTEIN_SEQUENCE_FILE)
    go_connections = gene_ontology_reader(GENE_ONTOLOGY_FILE, train_terms)

    # GO connection data processing
    subontologies = generate_subontologies(train_terms)
    go_connections = add_subontologies_to_go(go_connections, subontologies)
    go_connections = add_relationship_count_to_go(go_connections, "is_a")
    go_connections = add_relationship_count_to_go(go_connections, "part_of")
    go_connections = add_relationship_count_to_go(go_connections, "regulates")

    # Protein data processing
    train_proteins = train_sequences.join(train_tax.set_index("EntryID"), on="EntryID", how="left")
    train_proteins = protein_sequence_length(train_proteins)
    train_proteins = GO_counts_by_protein(train_proteins, train_terms)
    train_proteins = generate_amino_count(train_proteins)
    train_proteins = aspect_distributions(train_proteins, train_terms)

    # Storing (Load)
    store_data(train_proteins, protein_filepath)
    store_data(go_connections, go_filepath)
