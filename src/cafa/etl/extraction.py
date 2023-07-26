import re

import obonet
import pandas as pd
from Bio import SeqIO


def protein_sequence_reader(filename):
    """Expects a fasta filetype. Returns dataframe."""
    sequences = SeqIO.parse(filename, "fasta")
    # read in fasta data protein by protein
    ## create list of protein ids, descriptions & sequences
    protein_ids = []
    protein_descriptions = []
    protein_sequences = []
    for seq in sequences:
        protein_ids.append(seq.id)
        protein_descriptions.append(seq.description)
        protein_sequences.append(seq.seq)
    # assemble dataframe for train sequences
    d = {"EntryID": protein_ids, "description": protein_descriptions, "sequence": protein_sequences}
    train_sequences = pd.DataFrame(data=d)
    return train_sequences


def extract_rel(term_list, graph, rel="is_a"):
    """Extract relationships from graph df for each gene ontology."""
    extraction_regex = {
        "part_of": (r"part_of GO:\d{7}", r"GO:\d{7}"),
        "is_a": {"is_a": "is_a"},
        "regulates": (r"regulates GO:\d{7}", r"GO:\d{7}"),
    }
    rel_list = []
    for term in term_list:
        if rel == "is_a":
            rel_list.append(graph.nodes[term][rel]) if rel in graph.nodes[
                term
            ] else rel_list.append(None)
        elif rel == "part_of" or rel == "regulates":
            mapper = extraction_regex[rel]
            full_finder = re.compile(mapper[0])
            term_finder = re.compile(mapper[1])
            if "relationship" in graph.nodes[term]:
                rel_list.append(
                    set(
                        [
                            term_finder.search(x)[0] if full_finder.search(x) else None
                            for x in graph.nodes[term]["relationship"]
                        ]
                    )
                )
            else:
                rel_list.append("N/A")
    return rel_list


def gene_ontology_reader(filename, train_terms):
    """Expects an obonet filetype. Returns dataframe."""
    graph = obonet.read_obo(filename)
    # get unique gene ontology terms from train_terms df
    unique_GO_terms = list(train_terms.term.unique())
    # extract part_of, regulates, & is_a relationships for each unique_GO_term
    is_a = extract_rel(unique_GO_terms, graph, rel="is_a")
    part_of = extract_rel(unique_GO_terms, graph, rel="part_of")
    regulates = extract_rel(unique_GO_terms, graph, rel="regulates")
    # create df for go_term, is_a, part_of & regulates
    d = {"term": unique_GO_terms, "is_a": is_a, "part_of": part_of, "regulates": regulates}
    go_connections = pd.DataFrame(d)
    go_connections.sort_values(by="term", inplace=True)
    return go_connections


def dataframe_from_tsv(filepath):
    return pd.read_csv(filepath, sep="\t")
