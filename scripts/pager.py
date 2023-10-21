"""
This module provides a class interface to the PAGER API, which offers functionalities
to perform hypergeometric tests, retrieve enriched PAGs, and access various other data 
related to the PAGER database.
"""

import requests
import pandas as pd


class PAGER:
    def __init__(self):
        self.params = {}

    def run_pager(self, genes, **kwargs):
        """
        Connects to the PAGER API and performs a hypergeometric test to retrieve enriched PAGs.

        Parameters:
        - genes (list of str): A list of gene names to be analyzed.
        - **kwargs: Optional keyword arguments to override default parameters.
            Supported parameters include:
            - source (list of str): The sources for the PAGs. Default is all available sources.
            - type (str): Type of PAGs. Default is "All".
            - min_size (int): Minimum size for PAGs. Default is 1.
            - max_size (int): Maximum size for PAGs. Default is 2000.
            - similarity (int): Similarity threshold. Default is 0.
            - overlap (str): Overlap threshold. Default is "1".
            - organism (str): Organism for the analysis. Default is "All".
            - n_coco (str): Cohesion threshold. Default is "0".
            - p_value (float): P-value threshold. Default is 0.05.
            - fdr (float): False discovery rate threshold. Default is 0.05.

        Returns:
        - pd.DataFrame: A DataFrame containing the enriched PAGs returned by the PAGER API.

        Note:
        - This method constructs a payload by merging the default parameters with the provided kwargs and sends a POST request to the PAGER API.
        - Any provided kwargs will override the default parameters.
        """

        all_sources = [
            "BioCarta",
            "DSigDB",
            "GAD",
            "GeneSigDB",
            "GOA",
            "GOA_EXCL",
            "GTEx",
            "GWAS Catalog",
            "Isozyme",
            "KEGG_2021_HUMAN",
            "Microcosm Targets",
            "mirTARbase",
            "MSigDB",
            "NCI-Nature Curated",
            "NGS Catalog",
            "Pfam",
            "PharmGKB",
            "PheWAS",
            "Protein Lounge",
            "Reactome_2021",
            "Spike",
            "TargetScan",
            "HPA-normProtein",
            "HPA-PathologyAtlas",
            "HPA-CellAtlas",
            "HPA-RNAcon",
            "HPA-normRNA",
            "HPA-GTEx",
            "HPA-FANTOM5",
            "HPA-TCGA",
            "The Genes Reported in Articles Published by Cell",
            "I2D database, version 2.9",
            "GeoMx Cancer Transcriptome Atlas",
            "WikiPathway_2021",
            "CellMarker",
        ]

        # Default parameters
        default_params = {
            "source": all_sources,
            "type": "All",
            "min_size": 1,
            "max_size": 2000,
            "similarity": 0,
            "overlap": "1",
            "organism": "All",
            "n_coco": "0",
            "p_value": 0.05,
            "fdr": 0.05,
        }
        # Handle case-insensitivity for kwargs
        kwargs_lower = {k.lower(): v for k, v in kwargs.items()}

        # Override defaults with provided kwargs
        for key in default_params.keys():
            default_params[key] = kwargs_lower.get(key, default_params[key])

        # Create the payload
        params = {
            "genes": "%20".join(genes),
            "source": "%20".join(default_params["source"]),
            "type": default_params["type"],
            "ge": default_params["min_size"],
            "le": default_params["max_size"],
            "sim": str(default_params["similarity"]),
            "olap": str(default_params["overlap"]),
            "organism": default_params["organism"],
            "cohesion": str(default_params["n_coco"]),
            "pvalue": default_params["p_value"],
            "FDR": default_params["fdr"],
        }

        response = requests.post(
            "http://discovery.informatics.uab.edu/PAGER/index.php/geneset/pagerapi",
            data=params,
        )
        return pd.DataFrame(response.json())

    def path_member(self, PAG_IDs):
        """Connected to PAGER API to retrieve the membership of PAGs."""
        params = {"pag": ",".join(PAG_IDs)}
        response = requests.post(
            "http://discovery.informatics.uab.edu/PAGER/index.php/geneset/get_members_by_ids/",
            data=params,
        )
        return pd.DataFrame(response.json()["data"])

    def path_int(self, PAG_IDs):
        """Connected to PAGER API to retrieve the m-type relationships of PAGs."""
        params = {"pag": ",".join(PAG_IDs)}
        response = requests.post(
            "http://discovery.informatics.uab.edu/PAGER/index.php/pag_pag/inter_network_int_api/",
            data=params,
        )
        return pd.DataFrame(response.json()["data"])

    def path_reg(self, PAG_IDs):
        """Connected to PAGER API to retrieve the r-type relationships of PAGs."""
        params = {"pag": ",".join(PAG_IDs)}
        response = requests.post(
            "http://discovery.informatics.uab.edu/PAGER/index.php/pag_pag/inter_network_reg_api/",
            data=params,
        )
        return pd.DataFrame(response.json()["data"])

    def pag_ranked_gene(self, PAG_id):
        """Connected to PAGER API to retrieve RP-ranked genes with RP-score."""
        response = requests.get(
            f"http://discovery.informatics.uab.edu/PAGER/index.php/genesinPAG/viewgenes/{PAG_id}"
        )
        return pd.DataFrame(response.json()["gene"])

    def pag_gene_int(self, PAG_id):
        """Connected to PAGER API to retrieve gene interaction network."""
        response = requests.get(
            f"http://discovery.informatics.uab.edu/PAGER/index.php/pag_mol_mol_map/interactions/{PAG_id}"
        )
        return pd.DataFrame(response.json()["data"])

    def pag_gene_reg(self, PAG_id):
        """Connected to PAGER API to retrieve gene regulatory network."""
        response = requests.get(
            f"http://discovery.informatics.uab.edu/PAGER/index.php/pag_mol_mol_map/regulations/{PAG_id}"
        )
        return pd.DataFrame(response.json()["data"])

    def path_ngsea(self, genes, PAG_member):
        """Connected to PAGER API to generate the network-based GSEA result."""
        gene_exp_str = "\\t\\t".join(
            [
                f"{row[0]}\\t\\t{row[1]}\\t\\t\\t"
                for row in genes.itertuples(index=False)
            ]
        )
        pag_sets_str = "\\t\\t".join(
            [
                f"{row[0]}\\t\\t{row[1]}\\t\\t\\t"
                for row in PAG_member.itertuples(index=False)
            ]
        )

        params = {"geneExpStr": gene_exp_str, "PAGsetsStr": pag_sets_str}
        response = requests.post(
            "http://discovery.informatics.uab.edu/PAGER/index.php/geneset/ngseaapi/",
            data=params,
        )
        return pd.DataFrame(response.json()["data"])
