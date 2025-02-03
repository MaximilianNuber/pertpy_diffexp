from collections.abc import Sequence

import numpy as np
import pandas as pd
from lamin_utils import logger
from scipy.sparse import issparse

from ._base import LinearModelBase
from ._checks import check_is_integer_matrix
from ._conversion import _py_to_r, _r_to_py, _ad_to_rmat, _ad_to_dge


class EdgeR(LinearModelBase):
    """Differential expression test using EdgeR"""

    def _check_counts(self):
        check_is_integer_matrix(self.data)

    def fit(self, **kwargs):  # adata, design, mask, layer
        """Fit model using edgeR.

        Note: this creates its own AnnData object for downstream.

        Args:
            **kwargs: Keyword arguments specific to glmQLFit()
        """
        # For running in notebook
        # pandas2ri.activate()
        # rpy2.robjects.numpy2ri.activate()
        try:
            import rpy2.robjects.numpy2ri
            import rpy2.robjects.pandas2ri
            from rpy2 import robjects as ro
            from rpy2.robjects import numpy2ri, pandas2ri
            from rpy2.robjects.conversion import localconverter
            from rpy2.robjects.packages import importr
            rbase = importr("base")


        except ImportError:
            raise ImportError("edger requires rpy2 to be installed.") from None

        try:
            edger = importr("edgeR")
        except ImportError as e:
            raise ImportError(
                "edgeR requires a valid R installation with the following packages:\n"
                "edgeR, BiocParallel, RhpcBLASctl"
            ) from e

        # Convert dataframe
        '''
        with localconverter(ro.default_converter + numpy2ri.converter):
            expr = self.adata.X if self.layer is None else self.adata.layers[self.layer]
            if issparse(expr):
                expr = expr.T.toarray()
            else:
                expr = expr.T
        '''
        if "layer" in kwargs.keys():
            layer = kwargs["layer"]
            expr_r = _ad_to_rmat(self.adata, layer = layer)
            kwargs.pop("layer")
        # expr_r = ro.conversion.py2rpy(pd.DataFrame(expr, index=self.adata.var_names, columns=self.adata.obs_names))
        expr_r = _ad_to_rmat(self.adata)

        dge = edger.DGEList(counts=expr_r, samples=_py_to_r(self.adata.obs))

        # logger.info("Calculating NormFactors")
        dge = edger.calcNormFactors(dge)

        
        # logger.info("Estimating Dispersions")
        
        r_design = _py_to_r(self.design)
        r_design = rbase.as_matrix(r_design)

        dge = edger.estimateDisp(dge, design=r_design)

        # logger.info("Fitting linear model")
        fit = edger.glmQLFit(dge, design=r_design, **kwargs)

        # ro.globalenv["fit"] = fit
        self.fit = fit

    def _test_single_contrast(self, contrast: Sequence[float], **kwargs) -> pd.DataFrame:
        """Conduct test for each contrast and return a data frame

        Args:
            contrast: numpy array of integars indicating contrast i.e. [-1, 0, 1, 0, 0]
        """
        ## -- Check installations
        # For running in notebook
        # pandas2ri.activate()
        # rpy2.robjects.numpy2ri.activate()

        # ToDo:
        #  parse **kwargs to R function
        #  Fix mask for .fit()

        try:
            import rpy2.robjects.numpy2ri
            import rpy2.robjects.pandas2ri
            from rpy2 import robjects as ro
            from rpy2.robjects import numpy2ri, pandas2ri
            from rpy2.robjects.conversion import localconverter
            from rpy2.robjects.packages import importr
            from rpy2.robjects.packages import STAP

        except ImportError:
            raise ImportError("edger requires rpy2 to be installed.") from None

        try:
            edger = importr("edgeR")
            r_base = importr("base")
        except ImportError:
            raise ImportError(
                "edgeR requires a valid R installation with the following packages: " "edgeR, BiocParallel, RhpcBLASctl"
            ) from None

        # Convert vector to R, which drops a category like `self.design_matrix` to use the intercept for the left out.
        contrast_vec_r = _py_to_r(np.asarray(contrast))

        test = edger.glmQLFTest(self.fit, contrast = contrast_vec_r)
        de_res = edger.topTags(test, n = _py_to_r(np.inf), adjust_method = _py_to_r("BH"))
        de_res = r_base.as_data_frame(de_res)

        ### avoid STAP
        # r_string = '''
        # edger_test <- function(fit, contrast) {
        #     test = edgeR::glmQLFTest(fit, contrast=contrast)
        #     de_res = edgeR::topTags(test, n=Inf, adjust.method="BH", sort.by="PValue")$table
        #     return(de_res)
        # }
        # '''

        # test_pkg = STAP(r_string, "test_pkg")

       #  de_res = test_pkg.edger_test(self.fit, contrast_vec_r)
        de_res = _r_to_py(de_res)
        
        '''
        ro.globalenv["contrast_vec"] = contrast_vec_r
        # Test contrast with R
        ro.r(
            """
            test = edgeR::glmQLFTest(fit, contrast=contrast_vec)
            de_res =  edgeR::topTags(test, n=Inf, adjust.method="BH", sort.by="PValue")$table
            """
        )
        '''

        # Convert results to pandas
        # de_res = _r_to_py(ro.globalenv["de_res"])
        de_res.index.name = "variable"
        de_res = de_res.reset_index()

        return de_res.rename(columns={"PValue": "p_value", "logFC": "log_fc", "FDR": "adj_p_value"})