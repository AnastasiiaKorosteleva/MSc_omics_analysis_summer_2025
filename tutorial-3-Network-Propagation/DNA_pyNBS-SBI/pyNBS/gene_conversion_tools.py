################################################################
# ---------- Network Gene Name Conversion Functions ---------- #
################################################################

import requests
import re
import time
import pandas as pd
from typing import List, Optional, Tuple, Union


# Determine if id to be input is a valid gene name (does not contain parentheses or quotations or whitespace)
def exclude_id(name: str, bad_prefixes: Optional[List[str]] = None) -> bool:
    excluded_id_regex = re.compile(r"[(),\'\"\s\/|.<>\[\]]+")
    # Remove genes that may also have prefixes that we do not want (e.g., CHEBI)
    if bad_prefixes:
        for prefix in bad_prefixes:
            if name.startswith(prefix):
                return True
    return bool(excluded_id_regex.search(name))


# Remove the naming system prefix, if there is one
def get_identifier_without_prefix(string: str) -> Optional[str]:
    elements = string.split(":")
    if len(elements) == 2:
        return elements[1]
    elif len(elements) > 2:
        return None
    return string


# Construct string for batch query to MyGene.Info v3.0.0 API
def query_constructor(
    gene_list: List[str], 
    exclude_prefixes: Optional[List[str]] = None, 
    print_invalid_genes: bool = False
) -> Tuple[str, List[str], List[str]]:
    valid_query_genes = [
        get_identifier_without_prefix(gene)
        for gene in gene_list
        if exclude_id(gene, exclude_prefixes) is False
    ]
    invalid_query_genes = [
        gene for gene in gene_list if exclude_id(gene, exclude_prefixes) is True
    ]
    print(f"{len(valid_query_genes)} Valid Query Genes")
    if print_invalid_genes:
        print(f"{len(invalid_query_genes)} Invalid Query Genes: {invalid_query_genes}")
    else:
        print(f"{len(invalid_query_genes)} Invalid Query Genes")
    query_string = " ".join(valid_query_genes)
    return query_string, valid_query_genes, invalid_query_genes


# Function for posting batch query to MyGene.info v3.0.0 API
def query_batch(
    query_string: str, 
    tax_id: str = "9606", 
    scopes: str = "symbol,entrezgene,alias,uniprot", 
    fields: str = "symbol,entrezgene"
) -> List[dict]:
    query_split = query_string.split(" ")
    query_n = len(query_split)
    query_time = time.time()
    
    base_url = "https://mygene.info/v3/query"
    if query_n <= 1000:
        data = {
            "species": tax_id,
            "scopes": scopes,
            "fields": fields,
            "q": query_string
        }
        res = requests.post(base_url, data=data)
        json_response = res.json()
    else:
        chunks = (query_n + 999) // 1000
        query_chunks = [
            " ".join(query_split[i * 1000:(i + 1) * 1000]) for i in range(chunks)
        ]
        json_response = []
        for chunk in query_chunks:
            data = {
                "species": tax_id,
                "scopes": scopes,
                "fields": fields,
                "q": chunk
            }
            res = requests.post(base_url, data=data)
            json_response.extend(res.json())
    print(f"{len(json_response)} Matched query results")
    print(f"Batch query complete: {round(time.time() - query_time, 2)} seconds")
    return json_response


# Construct matched queries maps
def construct_query_map_table(
    query_result: List[dict], 
    query_genes: List[str], 
    display_unmatched_queries: bool = False
) -> Tuple[pd.DataFrame, dict, dict]:
    construction_time = time.time()
    matched_data = [
        [match.get("query"), match.get("_score"), match.get("symbol"), str(match.get("entrezgene"))]
        for match in query_result if match.get("entrezgene") and match.get("symbol")
    ]
    matched_genes = {match.get("query") for match in query_result if match.get("entrezgene") and match.get("symbol")}

    partial_match_genes = [gene for gene in query_genes if gene not in matched_genes]
    partial_match_results = [match for match in query_result if match.get("query") in partial_match_genes]
    
    if display_unmatched_queries:
        for entry in partial_match_results:
            print(entry)

    match_table = pd.DataFrame(
        data=matched_data, columns=["Query", "Score", "Symbol", "EntrezID"]
    ).set_index("Query")
    
    duplicate_matched_genes = [
        gene for gene in matched_genes if isinstance(match_table.loc[gene], pd.DataFrame)
    ]
    single_match_genes = [gene for gene in query_genes if gene not in duplicate_matched_genes]
    match_table_single = match_table.loc[single_match_genes]

    if duplicate_matched_genes:
        max_score_matches = [
            match_table.loc[gene].nlargest(1, "Score") for gene in duplicate_matched_genes
        ]
        match_table_duplicate_max = pd.concat(max_score_matches)
        match_table_trim = pd.concat([match_table_single, match_table_duplicate_max])
    else:
        match_table_trim = match_table_single.copy()

    query_to_symbol = match_table_trim["Symbol"].to_dict()
    query_to_entrez = match_table_trim["EntrezID"].to_dict()
    print(f"Query mapping table/dictionary construction complete: {round(time.time() - construction_time, 2)} seconds")
    return match_table_trim, query_to_symbol, query_to_entrez


# Remaining functions will follow a similar pattern for updates if required
