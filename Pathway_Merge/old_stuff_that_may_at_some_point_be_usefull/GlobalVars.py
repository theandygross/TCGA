'''
Created on Jan 28, 2013

@author: agross
'''
pathway_file = '/cellar/users/agross/Data/GeneSets/c2.cp.v3.0.symbols_edit.csv'

meta_cols = ['bcrpatientbarcode', 'gender', 'daystoinitialpathologicdiagnosis',
             'daystobirth', 'race', 'icdo3histology', 'vitalstatus', 
             'patientid', 'tissuesourcesite','daystodeath',
             'ageatinitialpathologicdiagnosis', 'daystolastfollowup',
             'ethnicity']
meta_keep = ['vitalstatus','gender','tissuesourcesite','age','race','deceased']
event_cols = ['daystotumorrecurrence','daystotumorprogression', 
              'daystoadditionalsurgerylocoregionalprocedure',
              'daystoadditionalsurgerymetastaticprocedure', 
              'daystonewtumoreventafterinitialtreatment']

#event_cols = ['daystotumorrecurrence','daystotumorprogression', 
#              'daystonewtumoreventafterinitialtreatment']

feature_dict = {'anatomicorgansubdivision' : 'organ_subdivision',
               'neoplasmhistologicgrade' : 'histo_grade',
               'personneoplasmcancerstatus' : 'neo_status',
               'tumorresidualdisease' : 'tumor_residual',
               'tumorstage' : 'tumor_stage',
               'distantmetastasispathologicspread' : 'metastasis',
               'histologicaltype' : 'hist_type',
               'lymphnodepathologicspread' : 'lymphnode_stage',
               'primarytumorpathologicspread' : 'tumor_grade',
               'neoplasmanatomicsubdivision' : 'tumor_subdivision',
               'tobaccosmokinghistoryindicator' : 'smoking_status',
               'alcoholhistorydocumented' : 'drinker',
               'numberpackyearssmoked' : 'pack_years',
               'hemoglobinresult' : 'hemoglobin',
               'serumcalciumresult' : 'calcium_level',
               'whitecellcountresult' : 'white_cell_count',
               'venousinvasion' : 'venous_invasion',
               'acutemyeloidleukemiacalgbcytogeneticsriskcategory' : 'aml_cytogenics_risk',
               'cytogeneticabnormalities.cytogeneticabnormality' : 'aml_cytogenics',
               'fluorescenceinsituhybridcytogeneticsmetaphasenucleusresultcount' : 
                    'metaphase_nucleous_count',
               'leukemiafrenchamericanbritishmorphologycode' : 'mophology',
               'psaresultmostrecent' : 'psa',
               'clinicalcqcf.tumorfocality': 'tumor_focality',
               'neoplasmdimension.neoplasmdepth' : 'neo_depth',
               'neoplasmdimension.neoplasmlength' : 'neo_length',
               'neoplasmdimension.neoplasmwidth' : 'neo_width',
               'height' : 'height',
               'weight' : 'weight',
               'pcttumorinvasion' : 'pct_tumor_invasion',
               'breastcarcinomaestrogenreceptorstatus': 'ER',
               'breastcarcinomaprogesteronereceptorstatus': 'PR',
               'labprocedureher2neuinsituhybridoutcometype': 'her2_fish',
               'labprocher2neuimmunohistochemistryreceptorstatus': 'her2_chem'}

followup_dict = {'easterncanceroncologygroup' : 'eastern_cancer_oncology_group',
                 'karnofskyperformancescore' : 'karnofsky_score',
                 'siteoftumorfirstrecurrence' : 'recurrence_site',
                 'primarytherapyoutcomesuccess' : 'outcome'}
