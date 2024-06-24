# Summary

## Goal of Project
The goal of this competition is to predict the gene ontologies for a set of proteins.

---

## Data Provided
Training Data:
* train_sequences.fasta - amino acid sequences for proteins in training set
* train_terms.tsv - the training set of proteins and corresponding annotated GO terms
* train_taxonomy.tsv - taxon ID for proteins in training set
* go-basic.obo - ontology graph structure

For reference:
* An OBO file is text-based file used by OBO-Edit, the open source, platform-independent application for viewing and editing ontologies.
* Fasta files are text-based file for representing either nucleotide sequences or amino acid (protein) sequences


## Training Data
* testsuperset.fasta - amino acid sequences for proteins on which the predictions should be made
* testsuperset-taxon-list.tsv - taxon ID for proteins in test superset (Note: you may need to use encoding="ISO-8859-1" to read this file in pandas)


## Submission Files
* IA.txt - Information Accretion for each term. This is used to weight precision and recall
* sample_submission.csv - a sample submission file in the correct format

---

## Data Handling

### Loading Data
We loaded in the training data provided by the competition. Of the 4 datasets, I was unfamiliar with two of the filetypes: .fasta and .obo

In order to read in the fasta file, we had to 
- read in the protein data sequence by sequence using a third-party package (SeqIO)
- stored the data (id, description, sequence) for each protein sequence into separate lists
- created a dataframe out of the lists.

In order to read in the obo file, we had to
- read in the gene ontology data using a third-party package (obonet)
- created a function to extract the “is_a”, “part_of”, and “regulates” relationship information for each gene ontology within the file
    - Challenges: “part_of” and “regulates” were both values in the “relationship” key, so had to use look regex mapper within the relationship key to extract the specific information being saught. Also needed to specify “N/A” when information wasn’t present.
- stored each relationship type within a list
- created a dataframe using the gene ontology terms and the relationship lists

Feature Engineering & Merging DFs
We did some feature engineering  for the the Gene Ontology Dataframe (go_connections):
    - Extracted subontologies (“aspect”) for each GO from train_terms & append to go_connections
    - Created a custom len_of_list function to take in list/series & return the length of the list/seriesgiven nuances of what an “empty” input could look like.
    - Used len_of_list function to take create a count column for the “is_a”, “part_of”, and “regulates”entries for each go_term in the go_connections df

We also decided to merge two protein dataframes together into one, called train_proteins, and engineereda few additional features to add to it:
    - Created a column with the sequence length for each protein sequence
    - Got gene ontology counts for each protein from train_terms
    - Got unique number of amino acids in each protein sequence
    - Get subontology (aspect) proportions for each protein from train_terms, and add as 3 aspect proportion columns

Export Dataframes
We exported the finalized dataframes to use for data exploration using a third party package (joblib).

---

## Data Exploration

### Loading Data
We loaded in the data we cleaned up, merged, and did feature engineering on during data handling.
- go_df: contains go terms, is_a relationships, part_of relationships, regulates relationships, aspect (subontology), is_a count, part_of count & regulates count
- proteins_df: contains protein ID, protein description, protein sequence, taxonomy ID, protein sequence length, go term counts per protein, & counts of unique amino acids, & aspect proportions (BPO, CCO, and MFO)
- terms_df: contains protein ID, go term, and aspect (subontology)

### Univariate Exploration
#### Exploration 1
- **Question**: What are the most common gene ontologies?
- **Code**: I pulled the top 50 go-terms with their value counts, then plotted a bar chart with the top 15 in light blue and the next 35 in blue. Then I calculated what proportion of the data the top 15 go-terms account for, and how many more go-terms it would take to double this proportion.
- **Answer**: There are a clear top 5 gene ontologies, but there is also a natural break at the top 15 most common gene ontologies (as seen in light blue). Each of the top 15 gene ontologies are associated with anywhere from ~34,000 - 93,000 out of 142246 proteins. Collectively, these top 15 gene ontologies are associated with 18% of the protein/gene ontology relationships. For relative measure, it takes the next 53 gene ontologies to account for the next 18% of the data!

#### Exploration 2
- **Question**: What is the subontology distribution for gene ontologies? ...for the top 15 most common gene ontologies?
- **Code**: I created a common go-term dataframe using the top 15 go-terms and their corresponding data from go_df. Then I created a pie chart from the aspect (subontology) normalized value counts for all go terms and just the top 15 go-terms.
- **Answer**: BPO (Biological Process Ontology) is by far the the most common subonotology, making up 67.6% of the gene ontologies. Next is MFO (Molecular Functions Ontology) making up 23% of the gene ontologies. Last is CCO (Cellular Components Ontology) making up a measly 6% of the gene ontologies. Interestingly, when you look at just the top 15 most common gene ontologies, the subontology breakdown for BPO and CCO drastically changes. BPO drops 40% and CCO goes from being in only 6% of the data to half of the gene ontologies.

#### Exploration 3 
- **Question**: What is the distribution of is_a_count for gene ontologies? … for the top 15 most common gene ontologies?
- **Code**: I created a function called plt_count_bar() to plot bar charts for the relationship counts for is_a, part_of and regulates. It takes in the specified dataframe and relationship, get the value counts in ascending order, and plots the data with specified labels using the .format() method. I also had to specify the xticks in the function to go by integer values up to the max count, otherwise they were grouping strangely. I then used plt_count_bar() to plot is_a count bar chart for all go-terms and the top 15 go_terms.
- **Answer**: Gene ontologies look to have anywhere from 0-10 is_a relationships with other gene ontologies. The distribution of gene ontology is_a counts is right skewed. The highest frequency at 1 is_a relationship, almost double any other is_a count. There are very few gene ontologies with 0 & 4-10 is_a relationships.When looking at the top 15 most common gene ontologies, they only have 0-2 is_a relationships. Unsurprisingly, the majority of them have 1 is_a relationship.

#### Exploration 4
- **Question**: What is the distribution of part_of_count for gene ontologies? … for the top 15 most common gene ontologies?
- **Code**: I used plt_count_bar() to plot part_of count bar chart for all go-terms and the top 15 go_terms.
- **Answer**: Gene ontologies look to have anywhere from 0-3 part_of relationships with other gene ontologies. The distribution of gene ontology part_of counts is right skewed. The highest frequency is at 0 part_of relationships, over 5x the next highest frequency. The frequency of gene ontologies with part_of relationships decay exponentially as the number of relationships increase. When looking at the top 15 most common gene ontologies, they only have 0 or 1 part_of relationships. Again, following the trend of the rest of the gene ontologies, the majority of them have 0 part_of relationships.

#### Exploration 5
- **Question**: What is the distribution of regulates_count for gene ontologies? … for the top 15 most common gene ontologies?
- **Code**: I used plt_count_bar() to plot regulates count bar chart for all go-terms and the top 15 go_terms.
- **Answer**: Gene ontologies look to have anywhere from 0-2 regulates relationships with other gene ontologies. The distribution of gene ontology regulates counts is right skewed. The highest frequency is at 0 regulates relationships, roughly 4x the next highest frequency. There are very few gene ontologies with 2 regulates relationships. When looking at the top 15 most common gene ontologies, they only have 0 or 1 regulates relationships with a majority of them having 0 regulates relationships.

#### Exploration 6
- **Question**: What are the most common protein IDs?
- **Code**: I pulled the top 50 proteins with their value counts, then plotted a bar chart with the top 15 in light blue and the next 35 in blue. Then I calculated what proportion of the data the top 15 proteins account for.
- **Answer**: There are a clear top 5 proteins, but we'll consider the top 15 most common proteins (as seen in light blue). Each of the top proteins have anywhere from ~540 - 810 out of 31466 gene ontologies. Collectively, these top 15 proteins make up 0.1% of the data, so even though they are the most common, they make up a very small portion of the data.

#### Exploration 7
- **Question**: What is the distribution of protein sequence length?
    - **Code**: I created a function called plt_count_hist() to plot a histogram. It takes in the dataframe, column name, number of bins wanted, and the distribution name for the title. It calculates the min and max values from the column, and will calculate a bin width using those and the number of bins wanted. Then it plots the histogram along with labels/titles using information from the function input and .format() method. I used the plt_count_hist() to plot the distribution of protein sequenc_length. I also had to rescale the y-axis in thousands.
- **Answer**: Protein sequence lengths are very much skewed to the right. The majority of the sequence lengths lie between 0-5000 amino acids. The bin with the highest frequency are proteins with sequence lengths between 0 and 1200.

#### Exploration 8
- **Question**: What is the distribution of gene ontology counts for proteins?
- **Code**: I used  plt_count_hist() to plot a histogram of protein gene ontology term_count.
- **Answer**:  Proteins gene ontology counts are very much skewed to the right. The majority of the proteins have between 0-100 gene ontologies associated with them. The bin with the highest frequency are proteins with 0-30 gene ontologies.

#### Exploration 9
- **Question**: What is the distribution of amino counts for proteins?
- **Code**: I used plt_count_bar to plot a bar chart of protein amino_counts. 
- **Answer**: Protein Sequence unique amino acid counts are very much skewed to the left. The majority of the proteins have between 18-20 unique amino acids in thier protein sequence. The bin with the highest frequency are proteins with 20 unique amino acids.

---

## Model
### Loading Data
We loaded in the data provided by the competition & two datasets provided by a competitor; Sergei Fironov.

### Training Data Modification
We started with over 31466 unique labels. In order to build a model, we chose torestrict to the top 1500 labels. The reason for this was partially due to competition limitations, as well as computational resource limitations.

In restricting the labels, we had to do some data modification. 
- We created a 2D array where the number of columns is 1500, and the numberof rows is equal to the number of training observations (i.e. proteins EntryIDs).
- We initialized everything to zero first, and then for each protein EntryID,we assign a 1 to a column if that ID is related to the current GO-term.


### Training
- We specified the neural network architecture aka the model structure. 
- We then fit the model with the training data in batch. 
- Finally, we looked at the model’s performance over the training data on a per batch basis. Generic metrics included were cross-entropy loss and binary accuracy.

### Submission
We follow instructions on Kaggle to submit model predictions on a test dataset. Our model had a score of 0.39399 on the test dataset landing us at 1430/1675 teams. For reference, the winning team had a score of 0.65143.

This was my first time ever building out a model and submitting to Kaggle. 

---

## Restructure of Project Code
After wrapping up the project in Jupiter notebooks, we restructured the data handling code in Python files. We did this to practice better coding.

We created a file with functions for..
- general data flow
- extraction of data from .fasta and .obonet files
- transformations for feature engineering.Restructure of Project Code
After wrapping up the project in Jupiter notebooks, we restructured the data handling code in Python files. We did this to practice better coding.

We created a file with functions for..
- general data flow
- extraction of data from .fasta and .obonet files
- Transformations for feature engineering
