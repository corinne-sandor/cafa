import pandas as pd


def len_of_list(x):
    """Take in list/series entry & return len of the list/series.
    Nonetypes should be marked as length 0."""
    if x!={None} and x!="N/A" and x is not None:
        return len(x)
    else:
        return 0


def generate_subontologies(train_terms):
    subontologies = train_terms[['term','aspect']].groupby('term').aspect.unique().apply(lambda x: ''.join(x))
    return pd.DataFrame(subontologies).reset_index()


def add_subontologies_to_go(go_connections, subontologies):
    return go_connections.join(subontologies.set_index('term'), on='term', how='left')


def add_relationship_count_to_go(go_connections, relationship="is_a"):
    rel_count = f"{relationship}_count"
    go_connections[rel_count] = go_connections[relationship].apply(len_of_list)
    return go_connections


def protein_sequence_length(train_proteins):
    train_proteins["sequence_length"] = train_proteins.sequence.apply(lambda x: len(x))
    return train_proteins


def GO_counts_by_protein(train_proteins, train_terms):
    GO_counts_protein = pd.DataFrame(train_terms.EntryID.value_counts().reindex()).reset_index()
    train_proteins = train_proteins.join(GO_counts_protein.set_index('EntryID'), on='EntryID', how='left')
    train_proteins.rename(columns={'count': 'term_count'}, inplace=True)
    return train_proteins


def generate_amino_count(train_proteins):
    train_proteins["amino_count"] = train_proteins.sequence.apply(lambda x: len(set(x)))
    return train_proteins


def aspect_distributions(train_proteins, train_terms):
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
    return train_proteins