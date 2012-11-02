'''
Created on Nov 1, 2012

@author: agross
'''

keeper_columns = ['bcrpatientbarcode', 'gender', 'daystoinitialpathologicdiagnosis',
                  'daystobirth', 'race', 'icdo3histology', 'vitalstatus', 
                  'patientid', 'tissuesourcesite','daystodeath',
                  'ageatinitialpathologicdiagnosis','daystolastfollowup',
                  'ethnicity']
keeper_columns = map(lambda s: 'patient.' + s, keeper_columns)