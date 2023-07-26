import joblib
import re
import pandas as pd
import matplotlib.pyplot as plt
import obonet
import seaborn as sb
from Bio import SeqIO

TRAIN_TERMS_FILE = '../data/cafa-5-protein-function-prediction/Train/train_terms.tsv'
TRAIN_TAX_FILE = '../data/cafa-5-protein-function-prediction/Train/train_taxonomy.tsv'
PROTEIN_SEQUENCE_FILE = '../data/cafa-5-protein-function-prediction/Train/train_sequences.fasta'

def protein_sequence_reader(filename):
    """Expects a fasta filetype. Returns dataframe.
    """
    sequence = SeqIO.parse(filename, "fasta")
    #read in fasta data protein by protein
    ## create list of protein ids, descriptions & sequences
    protein_ids=[]
    protein_descriptions=[]
    protein_sequences=[]
    for seq in sequences:
        protein_ids.append(seq.id)
        protein_descriptions.append(seq.description)
        protein_sequences.append(seq.seq)
    # assemble dataframe for train sequences
    d = {'EntryID': protein_ids,
         'description': protein_descriptions,
         'sequence': protein_sequences}
    train_sequences = pd.DataFrame(data=d)
    return train_sequences

def extract_rel(term_list, graph, rel="is_a"):
    """ Expects list of gene ontologies, obo graph, & relationship. 
    Extracts relationship information from obo graph. """
    extraction_regex = {
    "part_of": ("part_of GO:\d{7}", "GO:\d{7}"),
    "is_a": {"is_a": "is_a"},
    "regulates": ("regulates GO:\d{7}", "GO:\d{7}")
    }
    terms =[]
    rel_list = []
    for term in term_list:
        if rel == "is_a":
            rel_list.append(graph.nodes[term][rel]) if rel in graph.nodes[term] else rel_list.append(None)
        elif rel == "part_of" or rel == "regulates":
            mapper = extraction_regex[rel]
            full_finder = re.compile(mapper[0])
            term_finder = re.compile(mapper[1])
            if "relationship" in graph.nodes[term]:
                rel_list.append(set([term_finder.search(x)[0] if full_finder.search(x) else None for x in graph.nodes[term]["relationship"]]))
            else:
                rel_list.append("N/A")  
    return rel_list

def gene_ontology_reader(filename):
    """Expects a obo filetype. Returns dataframe.
    """
    #load in gene ontology data
graph = obonet.read_obo('../data/cafa-5-protein-function-prediction/Train/go-basic.obo')

# get unique go_terms from train_terms
unique_GO_terms = list(train_terms.term.unique())

#extract part_of, regulates, & is_a relationships for each unique_GO_term
a = extract_rel(unique_GO_terms, graph, rel='part_of')
b = extract_rel(unique_GO_terms, graph, rel='regulates')
c = extract_rel(unique_GO_terms, graph, rel='is_a')

# create df for go_term, is_a, part_of & regulates
d = {'term': unique_GO_terms, 'is_a': c, 'part_of': a, 'regulates': b}
go_connections = pd.DataFrame(d)
go_connections.sort_values(by='term', inplace=True)

################ Write similiar function for graph (obonet file). #####################

# Extraction
train_terms = pd.read_csv(TRAIN_TERMS_FILE, sep="\t")
train_tax = pd.read_csv(TRAIN_TAX_FILE, sep="\t")
train_sequences = protein_sequence_reader(PROTEIN_SEQUENCE_FILE)


# Transformation

# Storing (Load)


#pull subontologies for each unique GO from train_terms
subontologies = train_terms[['term','aspect']].groupby('term').aspect.unique().apply(lambda x: ''.join(x))
subontologies = pd.DataFrame(subontologies).reset_index()

#add subontologies to go_connections by joining on terms
go_connections = go_connections.join(subontologies.set_index('term'), on='term', how='left')

def len_of_list(x):
    """Take in list/series entry & return len of the list/series.
    Nonetypes should be marked as length 0."""
    if x!={None} and x!="N/A" and x is not None:
        return len(x)
    else:
        return 0

#add count columns for is_a, part_of, and regulates entries
go_connections['is_a_count'] = go_connections['is_a'].apply(len_of_list)
go_connections['part_of_count'] = go_connections['part_of'].apply(len_of_list)
go_connections['regulates_count'] = go_connections['regulates'].apply(len_of_list)

# join train_tax to train_sequences by EntryID
train_proteins = train_sequences.join(train_tax.set_index('EntryID'), on='EntryID', how='left')

#for each protein, get len of sequence id & number of associated GOs
train_proteins['sequence_length'] = train_proteins.sequence.apply(lambda x: len(x))

# GO counts for each protein from train_terms
GO_counts_by_protein = pd.DataFrame(train_terms.EntryID.value_counts().reindex()).reset_index()
# join counts column to train_proteins
train_proteins = train_proteins.join(GO_counts_by_protein.set_index('EntryID'), on='EntryID', how='left')
#rename count column to say term_count (for GO "terms")
train_proteins.rename(columns={'count': 'term_count'}, inplace=True);

#get unique # of amino acids in the protein sequence
train_proteins['amino_count'] = train_proteins.sequence.apply(lambda x: len(set(x)))

# Group train terms by EntryID
# get distribution of aspect proportions for each Entry ID
groups = train_terms.groupby("EntryID")
aspect_dist = groups.aspect.value_counts(normalize=True)

# convert aspect_dist into dataframe
aspect_dist = aspect_dist.reset_index()
aspect_dist = aspect_dist.set_index(['EntryID', 'aspect']).proportion.unstack(fill_value = 0).reset_index()

# add suffix: _prop to BPO, CCO, and MFO
cols = ['BPO', 'CCO', 'MFO']
aspect_dist = aspect_dist.rename(columns={c: c+'_prop' for c in aspect_dist.columns if c in cols})

# join aspect_dist to train_proteins by EntryID
train_proteins = train_proteins.join(aspect_dist.set_index('EntryID'), on='EntryID', how='left')