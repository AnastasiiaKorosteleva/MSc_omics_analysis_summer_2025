################################################################
# ---------- Network Gene Name Conversion Functions ---------- #
################################################################
import requests
import re
import time
import pandas as pd

# Determine if id to be input is a valid gene name (does not contain parentheses or quotations or whitespace)
def exclude_id(name, bad_prefixes=None):
    excluded_id_regex = re.compile(r'[(),\'\"\s\/\|\.<>]+')
    # Remove genes that may also have prefixes that we do not want (e.g. CHEBI)
    if bad_prefixes:
        for prefix in bad_prefixes:
            if name.startswith(prefix):
                return True
    return excluded_id_regex.search(name)

# Remove the naming system prefix, if there is one
def get_identifier_without_prefix(string):
    elements = string.split(':')
    length = len(elements)
    if length == 2:
        return str(elements[1])
    elif length > 2:
        return None
    else:
        return string

# Construct string for batch query to MyGene.Info v3.0.0 API
def query_constructor(gene_list, exclude_prefixes=None, print_invalid_genes=False):
    valid_query_genes = [get_identifier_without_prefix(gene) for gene in gene_list if not exclude_id(gene, exclude_prefixes)]
    invalid_query_genes = [gene for gene in gene_list if exclude_id(gene, exclude_prefixes)]
    print(len(valid_query_genes), "Valid Query Genes")
    if print_invalid_genes:
        print(len(invalid_query_genes), "Invalid Query Genes:")
        print(invalid_query_genes)
    else:
        print(len(invalid_query_genes), "Invalid Query Genes")
    query_string = ' '.join(valid_query_genes)
    return query_string, valid_query_genes, invalid_query_genes

# Function for posting batch query to MyGene.Info v3.0.0 API
def query_batch(query_string, tax_id='9606', scopes="symbol, entrezgene, alias, uniprot", fields="symbol, entrezgene"):
    query_split = query_string.split(' ')
    query_n = len(query_split)
    query_time = time.time()
    if query_n <= 1000:
        data = {
            'species': tax_id,
            'scopes': scopes,
            'fields': fields,
            'q': query_string
        }
        res = requests.post('http://mygene.info/v3/query', data)
        json = res.json()
    else:
        chunks = -(-query_n // 1000)  # Integer division with rounding up
        query_chunks = [' '.join(query_split[i * 1000:(i + 1) * 1000]) for i in range(chunks)]
        json = []
        for chunk in query_chunks:
            data = {
                'species': tax_id,
                'scopes': scopes,
                'fields': fields,
                'q': chunk
            }
            res = requests.post('http://mygene.info/v3/query', data)
            json.extend(res.json())
    print(len(json), 'Matched query results')
    print('Batch query complete:', round(time.time() - query_time, 2), 'seconds')
    return json

# Construct matched queries maps
def construct_query_map_table(query_result, query_genes, display_unmatched_queries=False):
    construction_time = time.time()
    matched_data, matched_genes = [], []
    for match in query_result:
        if match.get('entrezgene') and match.get('symbol'):
            matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), str(match.get('entrezgene'))])
            matched_genes.append(match.get('query'))
    partial_match_genes = [gene for gene in query_genes if gene not in matched_genes]
    partial_match_results = []
    for match in query_result:
        if match.get('query') in partial_match_genes:
            partial_match_results.append(match)
            if match.get('entrezgene'):
                matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), str(match.get('entrezgene'))])
            else:
                matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), match.get('entrezgene')])
    print('Queries with partial matching results found:', len(partial_match_results))
    if display_unmatched_queries:
        for entry in partial_match_results:
            print(entry)
    match_table = pd.DataFrame(data=matched_data, columns=['Query', 'Score', 'Symbol', 'EntrezID']).set_index('Query')
    duplicate_matched_genes = [gene for gene in matched_genes if isinstance(match_table.loc[gene], pd.DataFrame)]
    print(len(duplicate_matched_genes), "Queries with multiple matches found")
    single_match_genes = [gene for gene in query_genes if gene not in duplicate_matched_genes]
    match_table_single = match_table.loc[single_match_genes]
    if duplicate_matched_genes:
        max_score_matches = [match_table.loc[gene][match_table.loc[gene]['Score'] == match_table.loc[gene]['Score'].max()] for gene in duplicate_matched_genes]
        match_table_trim = pd.concat([match_table_single] + max_score_matches)
    else:
        match_table_trim = match_table_single.copy()
    query_to_symbol = match_table_trim['Symbol'].to_dict()
    query_to_entrez = match_table_trim['EntrezID'].to_dict()
    print('Query mapping table/dictionary construction complete:', round(time.time() - construction_time, 2), 'seconds')
    return match_table_trim, query_to_symbol, query_to_entrez

# Adjustments for Python 3 have been applied throughout the script.
# Remaining functions can be similarly updated with syntax or logic corrections.
