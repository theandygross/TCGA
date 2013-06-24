In[1]:

.. code:: python

    cd ../src

.. parsed-literal::

    /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/src


In[2]:

.. code:: python

    import os as os
    import pickle as pickle
    import pandas as pd
    
    from Reports.NotebookTools import *
    from Processing.Tests import *
    from Reports.Figures import *
    
    pd.set_option('precision',3)
    pd.set_option('display.line_width', 100)
    pd.set_option('display.width', 300)

In[3]:

.. code:: python

    result_path = '/scratch/TCGA/Firehose__2012_01_16/ucsd_analyses'
    run = sorted(os.listdir(result_path))[0]
    run = pickle.load(open('/'.join([result_path, run, 'RunObject.p']), 'rb'))

In[9]:

.. code:: python

    cancer = run.load_cancer('OV')
    clinical = cancer.load_clinical()
    
    mut = cancer.load_data('MAF')
    mut.uncompress()
    rna = cancer.load_data('mRNASeq')

In[13]:

.. code:: python

    draw_survival_curves(mut.features.ix['BIOCARTA_P38MAPK_PATHWAY'], clinical.survival.survival, show=True, ann='p')

Out[13]:

.. parsed-literal::

    <Reports.Figures.Show at 0x7b4add0>

In[14]:

.. code:: python

    draw_survival_curves(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'], clinical.survival.survival, show=True, ann='p')

Out[14]:

.. parsed-literal::

    <Reports.Figures.Show at 0x79f9190>

In[ ]:

.. code:: python

    draw_survival_curves(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'], clinical.survival.survival, show=True, ann='p')

In[65]:

.. code:: python

    cancer = run.load_cancer('LUSC')
    clinical = cancer.load_clinical()
    
    mut = cancer.load_data('MAF')
    mut.uncompress()
    
    draw_survival_curves(mut.features.ix['BIOCARTA_P38MAPK_PATHWAY'], clinical.survival.survival, show=True, ann='p')

Out[65]:

.. parsed-literal::

    <Reports.Figures.Show at 0x7a57d50>

In[74]:

.. code:: python

    pd.crosstab(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'],
                clinical.clinical.tumor_stage).plot(kind='bar')

Out[74]:

.. parsed-literal::

    <matplotlib.axes.AxesSubplot at 0x715d650>

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_00.png

In[73]:

.. code:: python

    draw_survival_curves(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'],
                         clinical.survival.survival, show=True, ann='p')

Out[73]:

.. parsed-literal::

    <Reports.Figures.Show at 0x715d950>

In[76]:

.. code:: python

    meth  = cancer.load_data('Methylation')

In[68]:

.. code:: python

    draw_survival_curves(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'], 
                         clinical.survival.survival, clinical.clinical.tumor_stage, show=True, ann='p')

Out[68]:

.. parsed-literal::

    <Reports.Figures.Show at 0x73b1650>

In[21]:

.. code:: python

    rna = cancer.load_data('mRNASeq')

In[23]:

.. code:: python

    from Processing.Helpers import *

In[39]:

.. code:: python

    r = extract_pc(rna.df.ix[run.gene_sets['REACTOME_METABOLISM_OF_NITRIC_OXIDE']].dropna(how='all', axis=1))

In[42]:

.. code:: python

    draw_survival_curves(r['pat_vec'], clinical.survival.survival, show=True, show_legend=False, std=1)

Out[42]:

.. parsed-literal::

    <Reports.Figures.Show at 0x75a9e10>

In[12]:

.. code:: python

    pathway_plot(mut.df.ix[run.gene_sets['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION']])

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_01.png

In[10]:

.. code:: python

    cancer = run.load_cancer('LUSC')
    clinical = cancer.load_clinical()
    
    mut = cancer.load_data('MAF')
    mut.uncompress()
    
    draw_survival_curves(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'], clinical.survival.survival, show=True, ann='p')

Out[10]:

.. parsed-literal::

    <Reports.Figures.Show at 0x4064110>

In[16]:

.. code:: python

    cancer = run.load_cancer('KIRC')
    clinical = cancer.load_clinical()
    
    mut = cancer.load_data('MAF')
    mut.uncompress()

In[25]:

.. code:: python

    global_vars = cancer.load_global_vars()

In[88]:

.. code:: python

    venn_pandas(mut.features.ix['REACTOME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'], clinical.clinical.metastasis=='m1')

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_02.png

In[22]:

.. code:: python

    draw_survival_curves(mut.features.ix['REACTO
    ME_SEMA3A_PLEXIN_REPULSION_SIGNALING_BY_INHIBITING_INTEGRIN_ADHESION'], clinical.survival.event_free_survival, show=True, ann='p')

Out[22]:

.. parsed-literal::

    <Reports.Figures.Show at 0x814ea90>

In[87]:

.. code:: python

    draw_survival_curves(mut.features.ix['BIOCARTA_P38MAPK_PATHWAY'], clinical.survival.survival, show=True, ann='p')

Out[87]:

.. parsed-literal::

    <Reports.Figures.Show at 0x78ab5d0>

In[23]:

.. code:: python

    emt = pd.read_csv('EMT_sig.csv', index_col=0, squeeze=True)
    emt.name = 'mirna'
    
    args = [run.data_path] + [cancer.name]
    f = '{}/stddata/{}/mirnaseq/illuminahiseq_mirnaseq/bcgsc_ca/Level_3/miR_gene_expression/data/data.txt'.format(*args)
    mirna = pd.read_table(f, index_col=0, header=None)
    mirna = mirna.T.set_index(['miRNA_ID', 'Hybridization REF'])
    mirna = mirna.sortlevel(level=0).ix['reads_per_million_miRNA_mapped']
    mirna = np.log2(mirna.astype(float)).replace(-np.inf, -3.) #close enough to 0
    mirna = mirna.T
    mirna.columns = pd.MultiIndex.from_tuples([(i[:12], i[13:15]) for i 
                                               in mirna.columns])
    mirna = mirna.T.xs('01', level=1).T #pandas bug

In[26]:

.. code:: python

    survival_test = 'survival_5y'
    covariates = ['age']
    cov_df = global_vars.join(clinical.clinical, how='outer')
    cov_df = cov_df[covariates]
    surv = clinical.survival[survival_test]
    test = SurvivalTest(surv, cov_df)
    test.name = survival_test

In[28]:

.. code:: python

    mir_res = run_feature_matrix(mirna, test, fp_cutoff=1.)

In[89]:

.. code:: python

    mir_res.head(10)

Out[89]:

.. parsed-literal::

                       Full      Full                                               Full  Univariate                    
                         LR      LR_q                                               fmla     hazzard         p         q
    hsa-mir-21     2.21e-07  2.31e-04                Surv(days, event) ~ feature + age\n        2.69  1.43e-06  7.49e-04
    hsa-mir-153-1   1.5e-06  7.82e-04                Surv(days, event) ~ feature + age\n        1.42  3.28e-06  8.87e-04
    hsa-mir-676    4.54e-06  1.58e-03                Surv(days, event) ~ feature + age\n        0.63  3.39e-06  8.87e-04
    hsa-mir-149    1.27e-05  3.31e-03  Surv(days, event) ~ feature + age + age:feature\n        1.55  5.24e-06  1.10e-03
    hsa-mir-101-1  2.78e-05  5.81e-03                Surv(days, event) ~ feature + age\n        0.42  3.94e-05  4.58e-03
    hsa-mir-130b   4.77e-05  8.20e-03                Surv(days, event) ~ feature + age\n        2.05  7.24e-06  1.15e-03
    hsa-mir-138-1  5.49e-05  8.20e-03                Surv(days, event) ~ feature + age\n        1.34  1.03e-04  6.19e-03
    hsa-mir-153-2  7.61e-05  9.95e-03                Surv(days, event) ~ feature + age\n        1.42  5.11e-05  4.86e-03
    hsa-mir-365-2  9.15e-05  1.06e-02                Surv(days, event) ~ feature + age\n        1.86  1.21e-04  6.19e-03
    hsa-mir-3651   0.000103  1.08e-02  Surv(days, event) ~ feature + age + age:feature\n        1.43  4.76e-05  4.86e-03

In[44]:

.. code:: python

    from Data.Containers import Dataset
    from Processing.Helpers import *

In[45]:

.. code:: python

    df = mirna.ix[mir_res[mir_res.Full.LR_q < .01].index].dropna(how='all', axis=1)

In[46]:

.. code:: python

    pc = extract_pc(df)
    pct_var, mirna_vec, pc1 = pc['pct_var'], pc['gene_vec'], pc['pat_vec']
    pc1.name = 'mirna expression'

In[82]:

.. code:: python

    mirna_vec.order().plot(kind='bar')

Out[82]:

.. parsed-literal::

    <matplotlib.axes.AxesSubplot at 0x7dcfa90>

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_03.png

In[49]:

.. code:: python

    test.full_test(pc1)

Out[49]:

.. parsed-literal::

    LR                                      2.46e-11
    feature_p                               5.68e-12
    fmla         Surv(days, event) ~ feature + age\n
    hazzard                                 2.72e+05
    dtype: object

In[51]:

.. code:: python

    tt = lambda a,b: kruskal_pandas(b,a)

In[70]:

.. code:: python

    rna = cancer.load_data('mRNASeq')

In[80]:

.. code:: python

    violin_plot_pandas(clinical.clinical.metastasis, pc1);

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_04.png

In[91]:

.. code:: python

    draw_survival_curves(to_quants(pc1, std=1), clinical.survival.survival_5y, show=True, ann='p')

Out[91]:

.. parsed-literal::

    <Reports.Figures.Show at 0x7741110>

In[124]:

.. code:: python

    draw_survival_curves(to_quants(pc1, q=.5), clinical.survival.survival_5y, clinical.clinical.tumor_stage, show=True, ann='p')

Out[124]:

.. parsed-literal::

    <Reports.Figures.Show at 0x94e9c10>

In[92]:

.. code:: python

    violin_plot_pandas(clinical.clinical.tumor_stage, pc1)

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_05.png

In[119]:

.. code:: python

    draw_survival_curves(to_quants(rna.df.ix['FN1'], std=1), surv, clinical.clinical.tumor_stage, show=True)

Out[119]:

.. parsed-literal::

    <Reports.Figures.Show at 0x93b2450>

In[111]:

.. code:: python

    draw_survival_curves(to_quants(rna.df.ix['VIM'], q=.25), surv, clinical.clinical.tumor_stage, show=True)

Out[111]:

.. parsed-literal::

    <Reports.Figures.Show at 0x8fdba50>

In[122]:

.. code:: python

    draw_survival_curves(to_quants(rna.df.ix['CDH1'], std=1)==-1, surv, clinical.clinical.tumor_stage, show=True)

Out[122]:

.. parsed-literal::

    <Reports.Figures.Show at 0x94c20d0>

In[96]:

.. code:: python

    tt = lambda a,b: kruskal_pandas(b,a)

In[116]:

.. code:: python

    s = screen_feature(pc1, tt, rna.features)

In[117]:

.. code:: python

    s.head()

Out[117]:

.. parsed-literal::

                                                                                    H     p     q
    BIOCARTA_ACE2_PATHWAY                                                         203  0.49  0.49
    REACTOME_LOSS_OF_NLP_FROM_MITOTIC_CENTROSOMES                                 203  0.49  0.49
    REACTOME_LAGGING_STRAND_SYNTHESIS                                             203  0.49  0.49
    REACTOME_JNK_PHOSPHORYLATION_AND_ACTIVATION_MEDIATED_BY_ACTIVATED_HUMAN_TAK1  203  0.49  0.49
    REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS                                   203  0.49  0.49

In[29]:

.. code:: python

    mir_res.head()

Out[29]:

.. parsed-literal::

                       Full      Full                                               Full  Univariate                    
                         LR      LR_q                                               fmla     hazzard         p         q
    hsa-mir-21     2.21e-07  2.31e-04                Surv(days, event) ~ feature + age\n        2.69  1.43e-06  7.49e-04
    hsa-mir-153-1   1.5e-06  7.82e-04                Surv(days, event) ~ feature + age\n        1.42  3.28e-06  8.87e-04
    hsa-mir-676    4.54e-06  1.58e-03                Surv(days, event) ~ feature + age\n        0.63  3.39e-06  8.87e-04
    hsa-mir-149    1.27e-05  3.31e-03  Surv(days, event) ~ feature + age + age:feature\n        1.55  5.24e-06  1.10e-03
    hsa-mir-101-1  2.78e-05  5.81e-03                Surv(days, event) ~ feature + age\n        0.42  3.94e-05  4.58e-03

In[192]:

.. code:: python

    p53_mut = mut.df.ix['TP53'].clip_upper(1.)
    p53_mut.name = 'p53_mut'
    p53_meth = meth.df.ix['TP53']
    p53_meth.name = 'p53_meth'
    p53_rna = rna.df.ix['TP53']
    p53_rna.name = 'p53_rna'

In[102]:

.. code:: python

    surv = clinical.survival.event_free_survival_5y
    age = clinical.clinical.age

In[103]:

.. code:: python

    import Data.Firehose as FH
    from Processing.Helpers import *
    from Figures.Survival import draw_survival_curve

In[104]:

.. code:: python

    gistic = FH.get_gistic_gene_matrix(run.data_path, cancer.name)

In[164]:

.. code:: python

    survival_test = 'survival_5y'
    covariates = ['age']
    cov_df = global_vars.join(clinical.clinical, how='outer')
    cov_df = cov_df[covariates]
    #cov_df[('cna', 'chrom_instability')] = global_vars[('cna', 'chrom_instability')]
    remerge = lambda s: '__'.join(s) if type(s) != str else s
    cov_df = cov_df.rename(columns=remerge)
    surv = clinical.survival[survival_test]
    test = SurvivalTest(surv, cov_df)
    test.name = survival_test
    #test.check_feature = lambda s: (len(s.unique()) == 2) and ((pd.crosstab(s, two_hit).stack() > 9).sum() > 2)
    test.check_feature = lambda s: s.value_counts()[0] < (len(s) - 10)
        
    def fp(feature):
        return get_cox_ph_ms(test.surv, feature, covariates=cov_df, return_val='p_haz',
                             formula='Surv(days, event) ~ ' + ' + '.join(list(cov_df.columns) + ['feature']))
    test.first_pass = fp

In[106]:

.. code:: python

    gv = run_feature_matrix(global_vars.T, test, fp_cutoff=1.)

In[75]:

.. code:: python

    gv.head()

Out[75]:

.. parsed-literal::

                                Full  Full                           Full  Univariate            
                                  LR  LR_q                           fmla     hazzard     p     q
    mutation rate_non         0.0332  0.45  Surv(days, event) ~ feature\n        0.00  0.04  0.41
             C->(G_A)         0.0576  0.45  Surv(days, event) ~ feature\n        0.25  0.04  0.41
             *CpG->T          0.0932  0.45  Surv(days, event) ~ feature\n        3.75  0.06  0.41
    cna      lesion_amp_high   0.105  0.45  Surv(days, event) ~ feature\n        1.02  0.15  0.48
    mutation rate_dbsnp         0.12  0.45  Surv(days, event) ~ feature\n         inf  0.09  0.45

In[76]:

.. code:: python

    s = pd.DataFrame({f: fisher_exact_test(global_vars['mutation']['rate_non'].dropna() == 0, feature) for f,feature in mut.features.iterrows()}).T
    s.sort(columns='p').head()

::

    ---------------------------------------------------------------------------
    ValueError                                Traceback (most recent call last)
    <ipython-input-76-9caa87c374a4> in <module>()
    ----> 1 s = pd.DataFrame({f: fisher_exact_test(global_vars['mutation']['rate_non'].dropna() == 0, feature) for f,feature in mut.features.iterrows()}).T
          2 s.sort(columns='p').head()
    
    /cellar/users/agross/semillon/epd/lib/python2.7/site-packages/pandas-0.11.0rc1-py2.7-linux-x86_64.egg/pandas/core/frame.pyc in __init__(self, data, index, columns, dtype, copy)
        392             mgr = self._init_mgr(data, index, columns, dtype=dtype, copy=copy)
        393         elif isinstance(data, dict):
    --> 394             mgr = self._init_dict(data, index, columns, dtype=dtype)
        395         elif isinstance(data, ma.MaskedArray):
        396             mask = ma.getmaskarray(data)
    
    /cellar/users/agross/semillon/epd/lib/python2.7/site-packages/pandas-0.11.0rc1-py2.7-linux-x86_64.egg/pandas/core/frame.pyc in _init_dict(self, data, index, columns, dtype)
        523 
        524         return _arrays_to_mgr(arrays, data_names, index, columns,
    --> 525                               dtype=dtype)
        526 
        527     def _init_ndarray(self, values, index, columns, dtype=None,
    
    /cellar/users/agross/semillon/epd/lib/python2.7/site-packages/pandas-0.11.0rc1-py2.7-linux-x86_64.egg/pandas/core/frame.pyc in _arrays_to_mgr(arrays, arr_names, index, columns, dtype)
       5304     # figure out the index, if necessary
       5305     if index is None:
    -> 5306         index = extract_index(arrays)
       5307     else:
       5308         index = _ensure_index(index)
    
    /cellar/users/agross/semillon/epd/lib/python2.7/site-packages/pandas-0.11.0rc1-py2.7-linux-x86_64.egg/pandas/core/frame.pyc in extract_index(data)
       5342 
       5343         if not indexes and not raw_lengths:
    -> 5344             raise ValueError('If use all scalar values, must pass index')
       5345 
       5346         if have_series or have_dicts:
    
    ValueError: If use all scalar values, must pass index

In[77]:

.. code:: python

    draw_survival_curves(global_vars['mutation']['rate_non'].dropna(), surv, p53_mut, show=True, std=.5)

Out[77]:

.. parsed-literal::

    <Reports.Figures.Show at 0x9fc7110>

In[96]:

.. code:: python

    fig, axs = subplots(1,2, figsize=(10,4))
    f = clinical.clinical.age.dropna().order()
    aa = pd.DataFrame({x: get_cox_ph_ms(surv, (f >= x)*1., return_val='p_haz') for x in f[1:-1:2]}).T
    aa.dropna().hazzard.plot(ax=axs[0])
    axs[0].set_ylabel('Hazzard Ratio')
    axs[0].set_xlabel('Cutoff')
    draw_survival_curve((f < aa.p.idxmin()).map({False: 'high', True: 'low'}), surv, ax=axs[1])

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_06.png

In[94]:

.. code:: python

    fig, axs = subplots(1,2, figsize=(10,4))
    f = global_vars['mutation']['rate_non'].dropna()
    aa = pd.DataFrame({x: get_cox_ph_ms(surv, (f >= x)*1., return_val='p_haz') for x in f[1:-1:2]}).T
    aa.dropna().hazzard.plot(ax=axs[0])
    axs[0].set_ylabel('Hazzard Ratio')
    axs[0].set_xlabel('Cutoff')
    draw_survival_curve((f < aa.p.idxmin()).map({False: 'high', True: 'low'}), surv, ax=axs[1])

.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_07.png

In[79]:

.. code:: python

    mut_high = (f < aa.p.idxmin()).map({False: 'high', True: 'low'})

In[80]:

.. code:: python

    [f for f in mut.features.index if 'TGF' in f ]

Out[80]:

.. parsed-literal::

    ['BIOCARTA_TGFB_PATHWAY', 'KEGG_TGF_BETA_SIGNALING_PATHWAY']

In[89]:

.. code:: python

    clinical.artificially_censor(8)

In[92]:

.. code:: python

    draw_survival_curves(mut.features.ix['BIOCARTA_TGFB_PATHWAY'], clinical.survival.event_free_survival_8y, mut_high, show=True)

Out[92]:

.. parsed-literal::

    <Reports.Figures.Show at 0x6a93f10>

In[48]:

.. code:: python

    fisher_exact_test(p53_mut, del_3p)

Out[48]:

.. parsed-literal::

    odds_ratio    1.34
    p             0.58
    dtype: float64

In[116]:

.. code:: python

    fhit = gistic.xs('FHIT', level=2).ix[0]
    del_3p = fhit < 0
    del_3p.name = 'del_3p'

In[28]:

.. code:: python

    surv = clinical.survival.survival_5y
    combo = ((del_3p==True)*1.).add((p53_mut > 0)*2.)
    #combo = combo[hpv[hpv = 'HPV-'].index][age[age < 85].index]
    #combo = combo.ix[combo.index.diff(true_index(hpv=='HPV+').union(true_index(age >= 85)))]
    rate = global_vars['mutation']['rate_non'].order().dropna()
    combo = combo[rate[to_quants(rate, std=3) == 0].index]
    combo = combo.dropna()
    combo = combo.map({0:'WT',1:'3p',2:'p53',3:'3p+p53'})
    combo.name = 'molecular_status'
    m = get_cox_ph_ms(surv, feature=None, covariates=pd.concat([combo], axis=1), interactions=True, return_val='model_desc')
    draw_survival_curve(combo, surv, colors=array(rcParams['axes.color_cycle'])[[0,4,1,2]])

.. parsed-literal::

    Call:  function (formula, data, weights, subset, na.action, init, control, 
        ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
        robust = FALSE, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
        ...) 
    {
        ties <- match.arg(ties)
        Call <- match.call()
        indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
            names(Call), nomatch = 0)
        if (indx[1] == 0) 
            stop("A formula argument is required")
        temp <- Call[c(1, indx)]
        temp[[1]] <- as.name("model.frame")
        special <- c("strata", "cluster", "tt")
        temp$formula <- if (missing(data)) 
            terms(formula, special)
        else terms(formula, special, data = data)
        if (is.R()) 
            m <- eval(temp, parent.frame())
        else m <- eval(temp, sys.parent())
        if (nrow(m) == 0) 
            stop("No (non-missing) observations")
        Terms <- terms(m)
        if (missing(control)) 
            control <- coxph.control(...)
        Y <- model.extract(m, "response")
        if (!inherits(Y, "Surv")) 
            stop("Response must be a survival object")
        type <- attr(Y, "type")
        if (type != "right" && type != "counting") 
            stop(paste("Cox model doesn't support \"", type, "\" survival data", 
                sep = ""))
        weights <- model.weights(m)
        data.n <- nrow(Y)
        strats <- attr(Terms, "specials")$strata
        if (length(strats)) {
            stemp <- untangle.specials(Terms, "strata", 1)
            if (length(stemp$terms) > 0) 
                Terms2 <- Terms[-stemp$terms]
            else Terms2 <- Terms
            if (length(stemp$vars) == 1) 
                strata.keep <- m[[stemp$vars]]
            else strata.keep <- strata(m[, stemp$vars], shortlabel = TRUE)
            strats <- as.numeric(strata.keep)
        }
        else Terms2 <- Terms
        timetrans <- attr(Terms, "specials")$tt
        if (length(timetrans)) {
            timetrans <- untangle.specials(Terms, "tt")
            ntrans <- length(timetrans$terms)
            if (missing(tt) || is.null(tt)) {
                tt <- function(x, time, riskset, weights) {
                    obrien <- function(x) {
                      r <- rank(x)
                      (r - 0.5)/(0.5 + length(r) - r)
                    }
                    unlist(tapply(x, riskset, obrien))
                }
            }
            if (is.function(tt)) 
                tt <- list(tt)
            if (is.list(tt)) {
                if (any(!sapply(tt, is.function))) 
                    stop("The tt argument must contain function or list of functions")
                if (length(tt) != ntrans) {
                    if (length(tt) == 1) {
                      temp <- vector("list", ntrans)
                      for (i in 1:ntrans) temp[[i]] <- tt[[1]]
                      tt <- temp
                    }
                    else stop("Wrong length for tt argument")
                }
            }
            else stop("The tt argument must contain function or list of functions")
            if (ncol(Y) == 2) {
                if (length(strats) == 0) {
                    sorted <- order(-Y[, 1], Y[, 2])
                    newstrat <- rep.int(0L, nrow(Y))
                    newstrat[1] <- 1L
                }
                else {
                    sorted <- order(strats, -Y[, 1], Y[, 2])
                    newstrat <- as.integer(c(1, 1 * (diff(strats[sorted]) != 
                      0)))
                }
                if (storage.mode(Y) != "double") 
                    storage.mode(Y) <- "double"
                counts <- .Call("coxcount1", Y[sorted, ], as.integer(newstrat))
                tindex <- sorted[counts$index]
            }
            else {
                if (length(strats) == 0) {
                    sort.end <- order(-Y[, 2], Y[, 3])
                    sort.start <- order(-Y[, 1])
                    newstrat <- c(1L, rep(0, nrow(Y) - 1))
                }
                else {
                    sort.end <- order(strats, -Y[, 2], Y[, 3])
                    sort.start <- order(strata, -Y[, 1])
                    newstrat <- c(1L, as.integer(diff(strats[sort.end]) != 
                      0))
                }
                if (storage.mode(Y) != "double") 
                    storage.mode(Y) <- "double"
                counts <- .Call("coxcount2", Y, as.integer(sort.start - 
                    1L), as.integer(sort.end - 1L), as.integer(newstrat))
                tindex <- counts$index
            }
            m <- m[tindex, ]
            Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
            type <- "right"
            strats <- factor(rep(1:length(counts$nrisk), counts$nrisk))
            weights <- model.weights(m)
            for (i in 1:ntrans) m[[timetrans$var[i]]] <- (tt[[i]])(m[[timetrans$var[i]]], 
                Y[, 1], strats, weights)
        }
        offset <- model.offset(m)
        if (is.null(offset) | all(offset == 0)) 
            offset <- rep(0, nrow(m))
        cluster <- attr(Terms, "specials")$cluster
        if (length(cluster)) {
            robust <- TRUE
            tempc <- untangle.specials(Terms2, "cluster", 1:10)
            ord <- attr(Terms2, "order")[tempc$terms]
            if (any(ord > 1)) 
                stop("Cluster can not be used in an interaction")
            cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
            Terms2 <- Terms2[-tempc$terms]
        }
        else {
            if (!missing(robust)) 
                warning("The robust option is depricated")
            else robust <- FALSE
        }
        attr(Terms2, "intercept") <- 1
        X <- model.matrix(Terms2, m)
        Xatt <- attributes(X)
        if (is.R()) {
            assign <- lapply(attrassign(X, Terms2)[-1], function(x) x - 
                1)
            xlevels <- .getXlevels(Terms2, m)
            contr.save <- attr(X, "contrasts")
        }
        else {
            assign <- lapply(attr(X, "assign")[-1], function(x) x - 
                1)
            xvars <- as.character(attr(Terms2, "variables"))
            xvars <- xvars[-attr(Terms2, "response")]
            if (length(xvars) > 0) {
                xlevels <- lapply(m[xvars], levels)
                xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
                if (length(xlevels) == 0) 
                    xlevels <- NULL
            }
            else xlevels <- NULL
            contr.save <- attr(X, "contrasts")
        }
        X <- X[, -1, drop = F]
        if (missing(init)) 
            init <- NULL
        pterms <- sapply(m, inherits, "coxph.penalty")
        if (any(pterms)) {
            pattr <- lapply(m[pterms], attributes)
            pname <- names(pterms)[pterms]
            ord <- attr(Terms, "order")[match(pname, attr(Terms, 
                "term.labels"))]
            if (any(ord > 1)) 
                stop("Penalty terms cannot be in an interaction")
            pcols <- assign[match(pname, names(assign))]
            fit <- coxpenal.fit(X, Y, strats, offset, init = init, 
                control, weights = weights, method = method, row.names(m), 
                pcols, pattr, assign)
        }
        else {
            if (method == "breslow" || method == "efron") {
                if (type == "right") 
                    fitter <- get("coxph.fit")
                else fitter <- get("agreg.fit")
            }
            else if (method == "exact") {
                if (type == "right") 
                    fitter <- get("coxexact.fit")
                else fitter <- get("agexact.fit")
            }
            else stop(paste("Unknown method", method))
            fit <- fitter(X, Y, strats, offset, init, control, weights = weights, 
                method = method, row.names(m))
        }
        if (is.character(fit)) {
            fit <- list(fail = fit)
            if (is.R()) 
                class(fit) <- "coxph"
            else oldClass(fit) <- "coxph"
        }
        else {
            if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
                vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
                msg <- paste("X matrix deemed to be singular; variable", 
                    paste(vars, collapse = " "))
                if (singular.ok) 
                    warning(msg)
                else stop(msg)
            }
            fit$n <- data.n
            fit$nevent <- sum(Y[, ncol(Y)])
            fit$terms <- Terms
            fit$assign <- assign
            if (is.R()) 
                class(fit) <- fit$method
            else oldClass(fit) <- fit$method[1]
            if (robust) {
                fit$naive.var <- fit$var
                fit$method <- method
                fit2 <- c(fit, list(x = X, y = Y, weights = weights))
                if (length(strats)) 
                    fit2$strata <- strats
                if (length(cluster)) {
                    temp <- residuals.coxph(fit2, type = "dfbeta", 
                      collapse = cluster, weighted = TRUE)
                    if (is.null(init)) 
                      fit2$linear.predictors <- 0 * fit$linear.predictors
                    else fit2$linear.predictors <- c(X %*% init)
                    temp0 <- residuals.coxph(fit2, type = "score", 
                      collapse = cluster, weighted = TRUE)
                }
                else {
                    temp <- residuals.coxph(fit2, type = "dfbeta", 
                      weighted = TRUE)
                    fit2$linear.predictors <- 0 * fit$linear.predictors
                    temp0 <- residuals.coxph(fit2, type = "score", 
                      weighted = TRUE)
                }
                fit$var <- t(temp) %*% temp
                u <- apply(as.matrix(temp0), 2, sum)
                fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u, 
                    control$toler.chol)$test
            }
            if (length(fit$coefficients) && is.null(fit$wald.test)) {
                nabeta <- !is.na(fit$coefficients)
                if (is.null(init)) 
                    temp <- fit$coefficients[nabeta]
                else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
                fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta], 
                    temp, control$toler.chol)$test
            }
            na.action <- attr(m, "na.action")
            if (length(na.action)) 
                fit$na.action <- na.action
            if (model) {
                if (length(timetrans)) {
                    m[[".surv."]] <- Y
                    m[[".strata."]] <- strats
                    stop("Time transform + model frame: code incomplete")
                }
                fit$model <- m
            }
            if (x) {
                Xatt$dim <- attr(X, "dim")
                Xatt$dimnames <- attr(X, "dimnames")
                Xatt$assign <- Xatt$assign[-1]
                attributes(X) <- Xatt
                fit$x <- X
                if (length(strats)) {
                    if (length(timetrans)) 
                      fit$strata <- strats
                    else fit$strata <- strata.keep
                }
            }
            if (y) 
                fit$y <- Y
        }
        if (!is.null(weights) && any(weights != 1)) 
            fit$weights <- weights
        fit$formula <- formula(Terms)
        if (length(xlevels) > 0) 
            fit$xlevels <- xlevels
        fit$contrasts <- contr.save
        if (any(offset != 0)) 
            fit$offset <- offset
        fit$call <- Call
        fit$method <- method
        fit
    }(formula = Surv(days, event) ~ 1, data = structure(list(molecular_status = structure(c("p53", 
    "p53", "p53", "p53", "3p+p53", "WT", "3p+p53", "3p+p53", "3p+p53", 
    "p53", "p53", "p53", "3p", "p53", "3p+p53", "p53", "3p+p53", 
    "p53", "3p+p53", "p53", "3p+p53", "p53", "3p+p53", "3p+p53", 
    "3p+p53", "p53", "p53", "p53", "p53", "p53", "p53", "p53", "p53", 
    "3p+p53", "3p+p53", "3p+p53", "3p+p53", "WT", "p53", "WT", "3p", 
    "p53", "p53", "p53", "p53", "3p+p53", "3p+p53", "p53", "WT", 
    "p53", "WT", "p53", "p53", "WT", "3p+p53", "p53", "WT", "p53", 
    "p53", "p53", "3p+p53", "3p+p53", "3p+p53", "p53", "p53", "p53", 
    "p53", "p53", "3p+p53", "p53", "p53", "p53", "3p+p53", "p53", 
    "p53", "3p", "p53", "p53", "3p+p53", "p53", "3p+p53", "p53", 
    "p53", "p53", "p53", "p53", "p53", "p53", "3p+p53", "3p+p53", 
    "p53", "p53", "p53", "p53", "3p", "WT", "3p+p53", "p53", "3p+p53", 
    "3p+p53", "WT", "p53", "p53", "p53", "p53", "3p+p53", "p53", 
    "p53", "p53", "p53", "3p+p53", "p53", "p53", "p53", "p53", "3p+p53", 
    "p53", "3p+p53", "p53", "WT", "p53", "3p+p53", "WT", "3p+p53", 
    "3p+p53", "p53", "p53", "p53", "3p+p53", "WT", "p53", "3p", "WT", 
    "p53", "p53", "3p+p53", "p53", "p53", "p53", "3p+p53", "p53", 
    "p53", "p53", "p53", "p53", "3p+p53", "p53", "WT", "p53", "p53", 
    "3p+p53", "3p+p53", "p53", "p53", "p53", "p53", "p53", "p53", 
    "p53", "p53", "p53", "3p", "3p+p53", "p53", "p53", "p53", "p53", 
    "p53", "3p+p53", "p53", "p53", "WT", "3p+p53", "p53", "p53", 
    "p53", "3p+p53", "p53", "p53", "p53", "3p+p53", "p53", "3p+p53", 
    "3p+p53", "3p+p53", "WT", "p53", "p53", "p53", "p53", "p53", 
    "p53", "p53", "p53", "p53", "p53", "p53", "WT", "3p+p53", "p53", 
    "p53", "p53", "p53", "p53", "p53", "p53", "WT", "WT", "p53", 
    "p53", "3p+p53", "p53", "p53", "p53", "WT", "p53", "3p+p53", 
    "p53", "p53", "3p+p53", "3p", "3p+p53", "p53", "WT", "3p+p53", 
    "p53", "p53", "3p+p53", "3p+p53", "3p+p53", "p53", "p53", "3p+p53", 
    "p53", "3p+p53", "p53", "p53", "p53", "p53", "3p+p53", "WT", 
    "p53", "WT", "3p+p53", "3p+p53", "3p", "p53", "WT", "p53", "p53", 
    "3p+p53", "p53", "WT", "p53", "3p+p53", "3p+p53", "p53", "p53", 
    "p53", "p53", "3p+p53", "p53", "p53", "p53", "p53", "WT", "3p+p53", 
    "3p", "p53", "p53", "3p+p53", "p53", "p53", "p53", "3p+p53", 
    "p53", "p53", "p53", "p53", "WT", "p53", "3p+p53", "3p+p53", 
    "3p+p53", "3p+p53", "p53", "p53", "3p+p53", "3p+p53", "p53", 
    "p53", "p53", "3p", "3p+p53", "p53", "p53", "p53", "WT", "p53", 
    "p53", "3p+p53", "3p+p53", "3p+p53", "p53", "3p+p53", "3p+p53"
    ), class = "AsIs"), days = structure(c(3.35342465753425, 3.41643835616438, 
    4.0958904109589, 0.167123287671233, 3.88493150684931, 1.54246575342466, 
    0.989041095890411, 5.0027397260274, 5.0027397260274, 4.06301369863014, 
    1.7972602739726, 5.0027397260274, 4.10684931506849, 2.70958904109589, 
    3.69315068493151, 2.80547945205479, 5.0027397260274, 5.0027397260274, 
    4.71232876712329, 1.66575342465753, 3.1972602739726, 5.0027397260274, 
    5.0027397260274, 4.81095890410959, 2.96164383561644, 0.83013698630137, 
    3.2027397260274, 5.0027397260274, 3.46575342465753, 4.7972602739726, 
    2.53972602739726, 0.50958904109589, 2.92876712328767, 5.0027397260274, 
    5.0027397260274, 5.0027397260274, 3.31232876712329, 1.03835616438356, 
    2.15890410958904, 5.0027397260274, 1.53972602739726, 2.69041095890411, 
    0.975342465753425, 0.323287671232877, 2.95068493150685, 1.26027397260274, 
    1.65479452054795, 0.517808219178082, 2.04931506849315, 3.71232876712329, 
    3.2986301369863, 0.227397260273973, 2.58904109589041, 1.26575342465753, 
    1.48219178082192, 4.5972602739726, 0.205479452054795, 2.83835616438356, 
    2.68767123287671, 3.24383561643836, 3.06027397260274, 2.39178082191781, 
    1.63013698630137, 0.676712328767123, 2.93972602739726, 1.13424657534247, 
    5.0027397260274, 5.0027397260274, 5.0027397260274, 5.0027397260274, 
    5.0027397260274, 5.0027397260274, 5.0027397260274, 5.0027397260274, 
    3.61369863013699, 4.78904109589041, 4.67945205479452, 4.75616438356164, 
    3.94794520547945, 4.02739726027397, 3.98630136986301, 3.74520547945205, 
    3.09041095890411, 2.30684931506849, 2.83835616438356, 2.4986301369863, 
    1.76164383561644, 1.24109589041096, 1.15068493150685, 0.534246575342466, 
    0.506849315068493, 0.906849315068493, 0.550684931506849, 0.613698630136986, 
    0.504109589041096, 0.493150684931507, 0.56986301369863, 0.391780821917808, 
    0.350684931506849, 0.443835616438356, 4.55068493150685, 5.0027397260274, 
    5.0027397260274, 2.44931506849315, 5.0027397260274, 1.86575342465753, 
    5.0027397260274, 5.0027397260274, 4.36986301369863, 5.0027397260274, 
    3.18082191780822, 4.1013698630137, 0.353424657534247, 4.52602739726027, 
    3.13424657534247, 3.1972602739726, 2.73698630136986, 1.38082191780822, 
    0.350684931506849, 0.367123287671233, 0.397260273972603, 0.345205479452055, 
    0.213698630136986, 0.257534246575342, 0.715068493150685, 1.92054794520548, 
    2.16164383561644, 2.18356164383562, 3.96164383561644, 3.37808219178082, 
    1.28219178082192, 2.23561643835616, 2.64383561643836, 4.11780821917808, 
    2.42465753424658, 1.57260273972603, 0.186301369863014, 4.54246575342466, 
    0.380821917808219, 2.77534246575342, 5.0027397260274, 3.25753424657534, 
    2.78904109589041, 4.84383561643836, 2.07945205479452, 5.0027397260274, 
    5.0027397260274, 5.0027397260274, 5.0027397260274, 0.635616438356164, 
    1.63561643835616, 0.96986301369863, 1.81643835616438, 3.46301369863014, 
    0.638356164383562, 1.86027397260274, 4.50684931506849, 5.0027397260274, 
    3.95068493150685, 0.526027397260274, 0.531506849315069, 0.652054794520548, 
    0.665753424657534, 0.654794520547945, 0.063013698630137, 0.520547945205479, 
    0.501369863013699, 0.495890410958904, 0.446575342465753, 0.4, 
    1.44657534246575, 1.5972602739726, 1.55616438356164, 3.62739726027397, 
    0.712328767123288, 5.0027397260274, 1.03835616438356, 3.76164383561644, 
    0.758904109589041, 0.284931506849315, 0.0986301369863014, 1.84931506849315, 
    2.24657534246575, 4.78082191780822, 1.35068493150685, 4.71506849315069, 
    4.32602739726027, 3.44931506849315, 4.84109589041096, 5.0027397260274, 
    5.0027397260274, 3.32328767123288, 1.62739726027397, 3.67397260273973, 
    3.79178082191781, 3.97534246575342, 2.15616438356164, 0.854794520547945, 
    1.43561643835616, 5.0027397260274, 5.0027397260274, 4.02739726027397, 
    3.18630136986301, 0.405479452054795, 4.84657534246575, 4.4958904109589, 
    2.34794520547945, 3.70684931506849, 4.75616438356164, 3.01369863013699, 
    0.0657534246575342, 0.0301369863013699, 3.96164383561644, 2.63561643835616, 
    5.0027397260274, 3.71506849315068, 0.0684931506849315, 5.0027397260274, 
    3.01917808219178, 1.38630136986301, 4.43835616438356, 2.24383561643836, 
    4.25479452054795, 3.08493150684932, 0.0849315068493151, 2.91506849315068, 
    5.0027397260274, 3.07945205479452, 2.83013698630137, 0.249315068493151, 
    0.246575342465753, 1.54794520547945, 1.5013698630137, 1.31506849315068, 
    0.0301369863013699, 1.71506849315068, 0.654794520547945, 0.0246575342465753, 
    4.74520547945205, 5.0027397260274, 2.96164383561644, 4.26849315068493, 
    1.08493150684932, 4.08767123287671, 0.0849315068493151, 3.16986301369863, 
    0.252054794520548, 3.75068493150685, 1.66575342465753, 3.5013698630137, 
    0.246575342465753, 2.41643835616438, 2.58082191780822, 2.24657534246575, 
    5.0027397260274, 3.02191780821918, 0.506849315068493, 2.5041095890411, 
    1.56438356164384, 0.558904109589041, 0.180821917808219, 2.39452054794521, 
    2.42191780821918, 1.79452054794521, 1.9041095890411, 1.87945205479452, 
    0.712328767123288, 2.5041095890411, 2.14520547945205, 2.32054794520548, 
    2.01917808219178, 2.0027397260274, 0.947945205479452, 1.76164383561644, 
    2.08767123287671, 5.0027397260274, 1.86027397260274, 5.0027397260274, 
    0.783561643835616, 2.86575342465753, 0.178082191780822, 0.449315068493151, 
    2.32328767123288, 4.06575342465753, 3.18082191780822, 0.164383561643836, 
    0.46027397260274, 1.20821917808219, 1.4986301369863, 0.331506849315069, 
    2.55068493150685, 3.32054794520548, 2.55342465753425, 0.0986301369863014, 
    0.397260273972603, 4.30958904109589, 5.0027397260274, 5.0027397260274, 
    5.0027397260274, 4.62465753424657, 0.53972602739726, 5.0027397260274, 
    1.72328767123288, 3.70958904109589, 5.0027397260274, 1.85205479452055
    ), class = "AsIs"), event = structure(c(1, 1, 0, 1, 0, 1, 1, 
    0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 
    0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 
    1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 
    0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 
    1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 
    0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 
    1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 
    0, 1, 1, 0, 1), class = "AsIs")), .Names = c("molecular_status", 
    "days", "event"), row.names = c("TCGA-04-1331", "TCGA-04-1332", 
    "TCGA-04-1336", "TCGA-04-1337", "TCGA-04-1338", "TCGA-04-1342", 
    "TCGA-04-1343", "TCGA-04-1346", "TCGA-04-1347", "TCGA-04-1348", 
    "TCGA-04-1349", "TCGA-04-1350", "TCGA-04-1356", "TCGA-04-1361", 
    "TCGA-04-1362", "TCGA-04-1364", "TCGA-04-1365", "TCGA-04-1367", 
    "TCGA-04-1514", "TCGA-04-1517", "TCGA-04-1525", "TCGA-04-1530", 
    "TCGA-04-1542", "TCGA-09-0366", "TCGA-09-0369", "TCGA-09-1659", 
    "TCGA-09-1661", "TCGA-09-1662", "TCGA-09-1665", "TCGA-09-1666", 
    "TCGA-09-1669", "TCGA-09-2044", "TCGA-09-2045", "TCGA-09-2049", 
    "TCGA-09-2050", "TCGA-09-2051", "TCGA-09-2053", "TCGA-09-2056", 
    "TCGA-10-0926", "TCGA-10-0927", "TCGA-10-0928", "TCGA-10-0931", 
    "TCGA-10-0933", "TCGA-10-0934", "TCGA-10-0935", "TCGA-10-0937", 
    "TCGA-10-0938", "TCGA-13-0714", "TCGA-13-0717", "TCGA-13-0720", 
    "TCGA-13-0723", "TCGA-13-0724", "TCGA-13-0726", "TCGA-13-0727", 
    "TCGA-13-0730", "TCGA-13-0751", "TCGA-13-0755", "TCGA-13-0761", 
    "TCGA-13-0762", "TCGA-13-0791", "TCGA-13-0792", "TCGA-13-0793", 
    "TCGA-13-0795", "TCGA-13-0800", "TCGA-13-0804", "TCGA-13-0807", 
    "TCGA-13-0883", "TCGA-13-0884", "TCGA-13-0885", "TCGA-13-0886", 
    "TCGA-13-0887", "TCGA-13-0889", "TCGA-13-0890", "TCGA-13-0891", 
    "TCGA-13-0893", "TCGA-13-0897", "TCGA-13-0899", "TCGA-13-0900", 
    "TCGA-13-0903", "TCGA-13-0904", "TCGA-13-0905", "TCGA-13-0906", 
    "TCGA-13-0910", "TCGA-13-0911", "TCGA-13-0912", "TCGA-13-0913", 
    "TCGA-13-0916", "TCGA-13-0919", "TCGA-13-0920", "TCGA-13-0923", 
    "TCGA-13-0924", "TCGA-13-1403", "TCGA-13-1404", "TCGA-13-1405", 
    "TCGA-13-1407", "TCGA-13-1408", "TCGA-13-1409", "TCGA-13-1410", 
    "TCGA-13-1411", "TCGA-13-1412", "TCGA-13-1477", "TCGA-13-1481", 
    "TCGA-13-1482", "TCGA-13-1483", "TCGA-13-1484", "TCGA-13-1487", 
    "TCGA-13-1488", "TCGA-13-1489", "TCGA-13-1491", "TCGA-13-1492", 
    "TCGA-13-1494", "TCGA-13-1495", "TCGA-13-1496", "TCGA-13-1497", 
    "TCGA-13-1498", "TCGA-13-1499", "TCGA-13-1501", "TCGA-13-1504", 
    "TCGA-13-1505", "TCGA-13-1506", "TCGA-13-1507", "TCGA-13-1509", 
    "TCGA-13-1510", "TCGA-13-1512", "TCGA-13-2060", "TCGA-20-0987", 
    "TCGA-20-0990", "TCGA-20-0991", "TCGA-23-1021", "TCGA-23-1023", 
    "TCGA-23-1024", "TCGA-23-1026", "TCGA-23-1027", "TCGA-23-1028", 
    "TCGA-23-1030", "TCGA-23-1031", "TCGA-23-1032", "TCGA-23-1110", 
    "TCGA-23-1116", "TCGA-23-1117", "TCGA-23-1118", "TCGA-23-1122", 
    "TCGA-23-1123", "TCGA-23-1124", "TCGA-23-2072", "TCGA-23-2077", 
    "TCGA-23-2078", "TCGA-23-2079", "TCGA-23-2081", "TCGA-24-0966", 
    "TCGA-24-0968", "TCGA-24-0970", "TCGA-24-0975", "TCGA-24-0979", 
    "TCGA-24-0980", "TCGA-24-0982", "TCGA-24-1103", "TCGA-24-1104", 
    "TCGA-24-1105", "TCGA-24-1413", "TCGA-24-1416", "TCGA-24-1417", 
    "TCGA-24-1418", "TCGA-24-1419", "TCGA-24-1422", "TCGA-24-1423", 
    "TCGA-24-1424", "TCGA-24-1425", "TCGA-24-1426", "TCGA-24-1427", 
    "TCGA-24-1428", "TCGA-24-1431", "TCGA-24-1434", "TCGA-24-1435", 
    "TCGA-24-1436", "TCGA-24-1463", "TCGA-24-1464", "TCGA-24-1466", 
    "TCGA-24-1469", "TCGA-24-1470", "TCGA-24-1471", "TCGA-24-1474", 
    "TCGA-24-1544", "TCGA-24-1545", "TCGA-24-1548", "TCGA-24-1549", 
    "TCGA-24-1551", "TCGA-24-1552", "TCGA-24-1553", "TCGA-24-1555", 
    "TCGA-24-1556", "TCGA-24-1557", "TCGA-24-1558", "TCGA-24-1560", 
    "TCGA-24-1562", "TCGA-24-1563", "TCGA-24-1564", "TCGA-24-1565", 
    "TCGA-24-1567", "TCGA-24-1603", "TCGA-24-1604", "TCGA-24-1614", 
    "TCGA-24-1616", "TCGA-24-2019", "TCGA-24-2024", "TCGA-24-2030", 
    "TCGA-24-2035", "TCGA-24-2038", "TCGA-24-2254", "TCGA-24-2260", 
    "TCGA-24-2261", "TCGA-24-2262", "TCGA-24-2267", "TCGA-24-2271", 
    "TCGA-24-2280", "TCGA-24-2281", "TCGA-24-2288", "TCGA-24-2289", 
    "TCGA-24-2290", "TCGA-24-2293", "TCGA-24-2298", "TCGA-25-1313", 
    "TCGA-25-1315", "TCGA-25-1316", "TCGA-25-1317", "TCGA-25-1318", 
    "TCGA-25-1319", "TCGA-25-1320", "TCGA-25-1321", "TCGA-25-1322", 
    "TCGA-25-1329", "TCGA-25-1623", "TCGA-25-1625", "TCGA-25-1626", 
    "TCGA-25-1627", "TCGA-25-1628", "TCGA-25-1630", "TCGA-25-1631", 
    "TCGA-25-1632", "TCGA-25-1633", "TCGA-25-1634", "TCGA-25-1635", 
    "TCGA-25-2042", "TCGA-25-2391", "TCGA-25-2392", "TCGA-25-2393", 
    "TCGA-25-2396", "TCGA-25-2398", "TCGA-25-2399", "TCGA-25-2400", 
    "TCGA-25-2401", "TCGA-25-2404", "TCGA-25-2408", "TCGA-25-2409", 
    "TCGA-29-2427", "TCGA-30-1853", "TCGA-30-1862", "TCGA-30-1891", 
    "TCGA-31-1950", "TCGA-31-1953", "TCGA-31-1959", "TCGA-36-1568", 
    "TCGA-36-1569", "TCGA-36-1570", "TCGA-36-1571", "TCGA-36-1574", 
    "TCGA-36-1575", "TCGA-36-1576", "TCGA-36-1577", "TCGA-36-1578", 
    "TCGA-36-1580", "TCGA-57-1582", "TCGA-57-1583", "TCGA-57-1584", 
    "TCGA-57-1993", "TCGA-59-2348", "TCGA-59-2350", "TCGA-59-2351", 
    "TCGA-59-2352", "TCGA-59-2354", "TCGA-59-2355", "TCGA-59-2363", 
    "TCGA-61-1728", "TCGA-61-1736", "TCGA-61-1919", "TCGA-61-1995", 
    "TCGA-61-1998", "TCGA-61-2000", "TCGA-61-2002", "TCGA-61-2003", 
    "TCGA-61-2008", "TCGA-61-2009", "TCGA-61-2012", "TCGA-61-2016", 
    "TCGA-61-2088", "TCGA-61-2092", "TCGA-61-2094", "TCGA-61-2095", 
    "TCGA-61-2097", "TCGA-61-2101", "TCGA-61-2102", "TCGA-61-2104", 
    "TCGA-61-2109", "TCGA-61-2110", "TCGA-61-2111", "TCGA-61-2113"
    ), class = "data.frame"))
    
    Null model
      log likelihood= -761.1305 
      n= 306 
    


.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_08.png

In[120]:

.. code:: python

    import Data.Firehose as FH

In[121]:

.. code:: python

    def is_disruptive(s):
        if 'fs' in s:
            return False
        if len(s) != 7:
            #print s
            return False
        if lo.Polarity[s[6]] == 'stop':
            return True
        aa = s[3:6]
        if int(aa) in range(163,196) + range(236, 252):
            if lo.Polarity[s[2]] == lo.Polarity[s[6]]:
                return True
        return False
    
    aa = pd.read_csv('/cellar/users/agross/Data/GeneSets/amino_acids.csv')
    lo = aa.set_index('Symbol').groupby(level=0).first()

In[123]:

.. code:: python

    reload(FH)

Out[123]:

.. parsed-literal::

    <module 'Data.Firehose' from 'Data/Firehose.pyc'>

In[124]:

.. code:: python

    p53 = FH.get_submaf(run.data_path, cancer.name, ['TP53'], fields='All').ix['TP53']
    status = pd.concat([combine(p53.Protein_Change.map(is_disruptive), p53.is_silent==0), p53.Tumor_Sample_Barcode], axis=1, keys=['status','barcode']).set_index('barcode')['status']
    status = (status == 'both').groupby(level=0).sum().clip_upper(1.)
    status = status.ix[mut.df.columns].fillna(-1).map({-1:'WT',0:'Non-Disruptive',1:'Disruptive'})

In[128]:

.. code:: python

    p53_mut.value_counts()

Out[128]:

.. parsed-literal::

    1    140
    0     37
    dtype: int64

In[127]:

.. code:: python

    draw_survival_curves(status, surv, show=True)

Out[127]:

.. parsed-literal::

    <Reports.Figures.Show at 0x10a29f90>

In[126]:

.. code:: python

    fig, axs = subplots(1,2, figsize=(8,3))
    
    cc = p53.set_index('Tumor_Sample_Barcode').Protein_Change
    cc = cc.groupby(level=0).agg(lambda s: s.ix[s.map(p53.Protein_Change.value_counts()).argmax()])
    cc = cc[cc.isin(true_index(cc.value_counts() > 4))]
    draw_survival_curve(cc, surv, ax=axs[0])
    axs[0].legend(loc='lower right', frameon=False)
    
    draw_survival_curve(status, surv, ax=axs[1])
    axs[1].legend(loc='lower left', frameon=False)
    fig.tight_layout()
    #plt.savefig('/cellar/users/agross/Desktop/Figures/TP53_characterization.pdf', transparent=True)

::

    ---------------------------------------------------------------------------
    RRuntimeError                             Traceback (most recent call last)
    <ipython-input-126-41b97dfd4343> in <module>()
          4 cc = cc.groupby(level=0).agg(lambda s: s.ix[s.map(p53.Protein_Change.value_counts()).argmax()])
          5 cc = cc[cc.isin(true_index(cc.value_counts() > 4))]
    ----> 6 draw_survival_curve(cc, surv, ax=axs[0])
          7 axs[0].legend(loc='lower right', frameon=False)
          8 
    
    /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/src/Figures/Survival.pyc in draw_survival_curve(feature, surv, q, std, **args)
         76         feature = to_quants(feature, q=q, std=std, labels=True)
         77     return feature
    ---> 78 
         79 def draw_survival_curve(feature, surv, q=.25, std=None, **args):
         80     feature = process_feature(feature, q, std)
    
    /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/src/Processing/Tests.pyc in get_cox_ph_ms(surv, feature, covariates, return_val, null_model, formula, get_model, interactions)
        109         fmla = robjects.Formula(formula)
        110 
    --> 111     if formula is not None:
        112         s = survival.coxph(fmla, df)
        113     else:
    
    /cellar/users/agross/semillon/epd/lib/python2.7/site-packages/rpy2/robjects/functions.pyc in __call__(self, *args, **kwargs)
         80                 v = kwargs.pop(k)
         81                 kwargs[r_k] = v
    ---> 82         return super(SignatureTranslatedFunction, self).__call__(*args, **kwargs)
    
    /cellar/users/agross/semillon/epd/lib/python2.7/site-packages/rpy2/robjects/functions.pyc in __call__(self, *args, **kwargs)
         32         for k, v in kwargs.iteritems():
         33             new_kwargs[k] = conversion.py2ri(v)
    ---> 34         res = super(Function, self).__call__(*new_args, **new_kwargs)
         35         res = conversion.ri2py(res)
         36         return res
    
    RRuntimeError: Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
      contrasts can be applied only to factors with 2 or more levels


.. image:: /cellar/users/agross/TCGA_Code/TCGA/Pathway_Merge/Notebooks/HNSC_validate_other_files/HNSC_validate_other_fig_09.png

In[117]:

.. code:: python

    draw_survival_curves(combine(p53_mut==1, del_3p), surv, show=True, ann='p')

Out[117]:

.. parsed-literal::

    <Reports.Figures.Show at 0x90950d0>

In[109]:

.. code:: python

    gistic

Out[109]:

.. parsed-literal::

    <class 'pandas.core.frame.DataFrame'>
    MultiIndex: 24174 entries, (1p36.33, 116983, ACAP3) to (Xp11.1, 6845, VAMP7)
    Columns: 358 entries, TCGA-18-3406 to TCGA-98-8023
    dtypes: float64(358)

In[124]:

.. code:: python

    def fisher_exact_test(hit_vec, response_vec):
        '''
        Wrapper to do a fischer's exact test on pandas Series
        ------------------------------------------------
        hit_vec: Series of labels (boolean, or (0,1))
        response_vec: Series of measurements (boolean, or (0,1))
        '''
        hit_vec.name = 'h' #crosstab can't handle multi_index
        response_vec.name = 'd' #so we use dummy names
        cont_table = pd.crosstab(hit_vec, response_vec)
        if (cont_table.shape != (2,2)):
            return pd.Series(index=['odds_ratio','p'])
        return pd.Series(fisher_exact(cont_table), index=['odds_ratio','p'])

In[98]:

.. code:: python

    draw_survival_curves(mut.features.ix['BIOCARTA_WNT_PATHWAY'] ,surv, two_hit, show=True)

::

    ---------------------------------------------------------------------------
    NameError                                 Traceback (most recent call last)
    <ipython-input-98-74bb12ebc162> in <module>()
    ----> 1 draw_survival_curves(mut.features.ix['BIOCARTA_WNT_PATHWAY'] ,surv, two_hit, show=True)
    
    NameError: name 'two_hit' is not defined

In[186]:

.. code:: python

    fisher_exact_test(cn.features.ix[g] == i, variable==c)

Out[186]:

.. parsed-literal::

    odds_ratio    1.60
    p             0.48
    dtype: float64

In[191]:

.. code:: python

    variable = p53_mut>0
    cn_association = pd.Series({(g, i,c, idx): val 
                                for c in variable.unique()
                                   for g in cn.features.index
                                   for i in ([1,2] if g[0] == 'Amplification' else [-1,-2])
                                   if (i in cn.features.ix[g].unique())
                                   for idx,val in fisher_exact_test(cn.features.ix[g] == i, variable==c).iteritems()
                                   })
    cn_association.index = pd.MultiIndex.from_tuples(cn_association.index)
    cn_association.index.names = ['feature','level','val','stat']

In[192]:

.. code:: python

    a = cn_association.unstack(level='stat')

In[201]:

.. code:: python

    a.sort(columns='p').head(20)

Out[201]:

.. parsed-literal::

                                                                                                                            odds_ratio     p
    feature                                                                                                    level val                    
    (Amplification, 14q13.3, (MBIP, MIPOL1, NKX2-8, PAX9, SLC25A21))                                            1    False        0.27  0.01
                                                                                                                     True         3.67  0.01
    (Deletion, 3p12.2, (GBE1,))                                                                                -1    False        0.35  0.01
                                                                                                                     True         2.85  0.01
    (Amplification, 18q12.1, (B4GALT6, DSC2, DSC3, DSG1, DSG2, DSG4, FAM59A, MEP1B, RNF138))                    1    False        0.26  0.01
                                                                                                                     True         3.78  0.01
    (Amplification, 18q11.2, Lesion)                                                                            1    False        0.30  0.02
                                                                                                                     True         3.34  0.02
    (Deletion, 5q11.2, Lesion)                                                                                 -1    False        0.39  0.02
                                                                                                                     True         2.57  0.02
    (Amplification, 14q21.1, (FOXA1,))                                                                          1    False        0.32  0.02
                                                                                                                     True         3.13  0.02
    (Amplification, 14q13.3, Lesion)                                                                            1    False        0.31  0.02
                                                                                                                     True         3.23  0.02
    (Amplification, 1q21.2, Lesion)                                                                             2    False        0.00  0.03
                                                                                                                     True          inf  0.03
    (Deletion, 3p12.1, Lesion)                                                                                 -1    False        0.41  0.03
                                                                                                                     True         2.43  0.03
    (Deletion, 9p21.3, (C9orf53, CDKN2A, CDKN2B, DMRTA1, ELAVL2, IFNE, KIAA1797, KLHL9, MLLT3, MTAP, PTPLAD2)) -2    False        2.94  0.04
                                                                                                                     True         0.34  0.04

In[203]:

.. code:: python

    cn_res = run_feature_matrix(cn.features, test, fp_cutoff=1.)

In[205]:

.. code:: python

    cn_res.head()

Out[205]:

.. parsed-literal::

                                                                                                                                                                                                                 Full  Full                                               Full  Univariate            
                                                                                                                                                                                                                   LR  LR_q                                               fmla     hazzard     p     q
    Amplification 3q26.32 (KCNMB2, KCNMB3, PIK3CA, TBL1XR1, ZMAT3)                                                                                                                                            0.00561  0.37  Surv(days, event) ~ feature + age + age:feature\n        0.69  0.01  0.66
                  3q27.1  (ABCC5, ABCF3, ALG3, AP2M1, B3GNT5, CAMK2N2, CLCN2, DVL3, ECE2, EIF2B5, EIF4G1, EPHB3, FAM131A, KLHL24, LAMP3, MAGEF1, MAP6D1, MCCC1, MCF2L2, PARL, POLR2H, PSMD2, VWA5B2, YEATS2)  0.00626  0.37  Surv(days, event) ~ feature + age + age:feature\n        0.70  0.01  0.66
                  3q26.33 (ACTL6A, ATP11B, CCDC39, DCUN1D1, DNAJC19, FXR1, GNB4, MFN1, MRPL47, NDUFB5, PEX5L, SOX2, TTC14, USP13, ZNF639)                                                                     0.00802  0.37  Surv(days, event) ~ feature + age + age:feature\n        0.70  0.01  0.66
                          Lesion                                                                                                                                                                              0.00932  0.37  Surv(days, event) ~ feature + age + age:feature\n        0.74  0.04  0.66
    Deletion      2q22.1  Lesion                                                                                                                                                                               0.0102  0.37  Surv(days, event) ~ feature + age + age:feature\n        1.05  0.85  0.98

In[206]:

.. code:: python

    cn_res = run_feature_matrix(rppa.features, test, fp_cutoff=1.)

In[207]:

.. code:: python

    cn_res.head()

Out[207]:

.. parsed-literal::

                                            Full  Full                                               Full  Univariate                
                                              LR  LR_q                                               fmla     hazzard         p     q
    protiens (MSH6, MSH6-R-C)           0.000858  0.21  Surv(days, event) ~ feature + age + age:feature\n        0.49  7.88e-04  0.26
             (ATM, ATM-R-C)              0.00151  0.21  Surv(days, event) ~ feature + age + age:feature\n        0.80  1.56e-01  0.78
             (MSH2, MSH2-M-C)            0.00187  0.21  Surv(days, event) ~ feature + age + age:feature\n        0.45  1.89e-03  0.31
             (STAT5A, STAT5-alpha-R-V)   0.00393  0.21  Surv(days, event) ~ feature + age + age:feature\n        0.80  2.15e-01  0.86
    phos_pc  CHEK2                        0.0043  0.21                Surv(days, event) ~ feature + age\n        0.00  8.23e-03  0.34

In[223]:

.. code:: python

    draw_survival_curves(rppa.features.ix[cn_res.index[0]], surv, show=True, show_legend=False)

Out[223]:

.. parsed-literal::

    <Reports.Figures.Show at 0xf2a6990>

In[225]:

.. code:: python

    draw_survival_curves(mut.features.ix['ANK2'], surv, p53_mut, show=True, show_legend=False)

Out[225]:

.. parsed-literal::

    <Reports.Figures.Show at 0xa055bd0>

In[226]:

.. code:: python

    p53_mut.name = 'p53'

In[227]:

.. code:: python

    survival_test = 'survival_5y'
    covariates = ['age','p53']
    cov_df = global_vars.join(clinical.clinical, how='outer').join(p53_mut)
    cov_df = cov_df[covariates]
    cov_df[('mutation', 'rate_non')] = global_vars[('mutation', 'rate_non')]
    remerge = lambda s: '__'.join(s) if type(s) != str else s
    cov_df = cov_df.rename(columns=remerge)
    surv = clinical.survival[survival_test]
    test = SurvivalTest(surv, cov_df)
    test.name = survival_test
    #test.check_feature = lambda s: (len(s.unique()) == 2) and ((pd.crosstab(s, two_hit).stack() > 9).sum() > 2)
    test.check_feature = lambda s: s.value_counts()[0] < (len(s) - 10)
        
    def fp(feature):
        return get_cox_ph_ms(test.surv, feature, covariates=cov_df, return_val='p_haz',
                             formula='Surv(days, event) ~ ' + ' + '.join(list(cov_df.columns) + ['feature']))
    test.first_pass = fp

In[229]:

.. code:: python

    mut_res = run_feature_matrix(mut.features, test, fp_cutoff=1.)

In[230]:

.. code:: python

    mut_res.head()

Out[230]:

.. parsed-literal::

                                                        Full  Full                                               Full  Univariate                
                                                          LR  LR_q                                               fmla     hazzard         p     q
    TEX15                                           0.000136  0.15  Surv(days, event) ~ feature + age + mutation__...        0.49  3.45e-01  1.00
    WHSC1L1                                         0.000316  0.18  Surv(days, event) ~ feature + age + p53 + feat...        1.00  9.98e-01  1.00
    COL5A2                                          0.000677  0.25  Surv(days, event) ~ feature + age + p53 + muta...        0.98  9.75e-01  1.00
    REACTOME_GLUTATHIONE_CONJUGATION                0.000997  0.28  Surv(days, event) ~ feature + mutation__rate_n...        4.13  3.14e-04  0.34
    REACTOME_METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS   0.00137  0.28  Surv(days, event) ~ feature + age + p53 + muta...        3.14  1.15e-01  1.00

In[ ]:

.. code:: python

    pd.crosstab(

In[240]:

.. code:: python

    draw_survival_curves(mut.features.ix['REACTOME_GLUTATHIONE_CONJUGATION'], surv, p53_mut, show=True)

Out[240]:

.. parsed-literal::

    <Reports.Figures.Show at 0xf2a1b90>

