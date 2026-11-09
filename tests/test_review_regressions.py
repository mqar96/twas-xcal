import numpy as np
import pandas as pd
import pytest
from xcal.core import run_xcal, xcal_gene, xcal_from_moments, harmonize, estimate_c

def test_coverage_excludes_panel_monomorphic_weight_and_applies_threshold():
    rng=np.random.default_rng(108)
    G=rng.binomial(2,.3,(300,4)).astype(float);G[:,0]=0
    snps=np.array(['s0','s1','s2','s3'])
    weights=pd.DataFrame({'gene':['A','A'],'snp':['s0','s1'],'weight':[10.,1.]})
    gwas=pd.DataFrame({'snp':snps,'z':[50.,2.,0.,0.],'eaf':[.3]*4})
    row=run_xcal(weights,gwas,G,snps,np.arange(4),c=.1).iloc[0]
    assert row.n_snps_used==1 and row.n_monomorphic_in_panel==1
    assert np.isclose(row.weight_coverage,1/101)
    filtered=run_xcal(weights,gwas,G,snps,np.arange(4),c=.1,min_coverage=.5).iloc[0]
    assert filtered.status=='insufficient_coverage'

def test_nonfinite_denominator_is_rejected():
    with pytest.raises(ValueError,match='finite and positive'):
        xcal_from_moments(1.,np.inf,1.,c=.1)

def test_direct_gene_interface_warns_about_panel_sd_fallback():
    G=np.random.default_rng(109).binomial(2,.3,(100,2))
    with pytest.warns(RuntimeWarning,match='panel SDs'):
        xcal_gene([1.,1.],[1.,1.],G,c=.1)

def test_harmonize_keeps_effect_frequency_aligned_and_is_stable():
    gwas=pd.DataFrame({'snp':['s0'],'effect_allele':['G'],'other_allele':['A'],'z':[2.],'eaf':[.2]})
    panel=pd.DataFrame({'snp':['s0'],'a1':['A'],'a2':['G']})
    out,_=harmonize(gwas,panel,value_col='z')
    assert out.z.iloc[0]==-2 and np.isclose(out.eaf.iloc[0],.8)
    assert np.isclose(out.original_eaf.iloc[0],.2)
    again,_=harmonize(out,panel,value_col='z')
    assert np.isclose(again.eaf.iloc[0],.8) and again.z.iloc[0]==-2

def test_ldsc_cannot_report_an_se_from_one_or_invalid_blocks():
    rng=np.random.default_rng(110)
    for bad in [1,True,-1,1.5]:
        with pytest.raises(ValueError,match='n_blocks'):
            estimate_c(rng.normal(size=2000),rng.uniform(1,10,2000),n_blocks=bad)
