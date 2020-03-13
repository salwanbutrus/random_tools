import numpy as np
import scanpy as sc
from anndata import AnnData

def module_score(
    adata:AnnData,
    genes_use: list,
    score_name: Optional[str] = None,
    verbose: bool = True):
    
    """\
   
    Compute module scores for all genes in adata as described in methods of RGC-dev paper.
    
    
    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    genes_use
        list of genes in module of interest
    score_name
        Name endowed to the module score to be computed
        e.g. "Mod1"
    verbose
        Inform user of fraction of module genes that are in adata
        
    Returns
    -------
    adata with a new .obs called score_name
    
    """
    
    if (score_name==None):
        score_name = str(input("Provide a name for this score (no spaces): "))
        
    genes_use0 = genes_use
    genes_use = list(set(genes_use).intersection(adata.var_names))#genes that are both in module and `adata`
    
     
    if (len(genes_use) == 0):
        raise ValueError("Error : Must provide a list of genes that are present in the data")
        
    
    if (verbose):
        if(len(genes_use0) > len(genes_use)):
            n = len(genes_use0) - len(genes_use)
            print("Note that", n, "genes (",n*100/len(genes_use0), " %) in your module do not exist in the data set." )
    
    
    
    adata_score = adata.copy()
    adata_score = adata[:,genes_use]
    
    counts_modgenes = adata_score.X.toarray() #all cells, module genes
    counts_all = adata.X.toarray() #all cells, all genes
    scores = np.mean(counts_modgenes, axis=1) - np.mean(counts_all, axis=1) #(row means of counts_modgenes ) - (row means of counts_all)
    
    adata.obs[score_name] = scores

    
    

