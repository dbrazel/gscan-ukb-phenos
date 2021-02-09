# Script to extract all of the GSCAN GWAS phenotypes to a CSV file, without
# the GWAS specific QC and filtering steps (eg, removing related individuals)

import pandas as pd

ukb_vars = pd.read_csv('ukb_gscan_vars.csv', index_col=False, quotechar='"')

out = pd.DataFrame(ukb_vars['f.eid'])

# CPD - cigarettes per day - by integrating across fields
# -10 -- Less than one a day, -1 -- Do not know, -3 -- Prefer not to answer
# Answers greater than 60 are considered invalid and are set to missing
# GWAS Version: Answers are binned - 1 = 1-5, 2 = 6-15, 3 = 16-25, 4 = 26-35, 5 = 36+
# Exome Version: Answers are binned - 1 = 1-10, 2 = 11-20, 3 = 21-30, 4 = >30

out['qcpd'] = ukb_vars['f.2887.0.0'].fillna(ukb_vars['f.3456.0.0']).fillna(
    ukb_vars['f.6183.0.0'])

out['qcpd'] = out['qcpd'].replace([-10, -1, -3], None)
out.loc[out.qcpd > 60, 'qcpd'] = None

out['cpd_gwas'] = out['qcpd']

out.loc[(out['cpd_gwas'] >= 1) & (out['cpd_gwas'] <= 5), 'cpd_gwas'] = 1
out.loc[(out['cpd_gwas'] >= 6) & (out['cpd_gwas'] <= 15), 'cpd_gwas'] = 2
out.loc[(out['cpd_gwas'] >= 16) & (out['cpd_gwas'] <= 25), 'cpd_gwas'] = 3
out.loc[(out['cpd_gwas'] >= 26) & (out['cpd_gwas'] <= 35), 'cpd_gwas'] = 4
out.loc[out['cpd_gwas'] > 35, 'cpd_gwas'] = 5

out['cpd_exome'] = out['qcpd']

out.loc[(out['cpd_exome'] >= 1) & (out['cpd_exome'] <= 10), 'cpd_exome'] = 1
out.loc[(out['cpd_exome'] >= 11) & (out['cpd_exome'] <= 20), 'cpd_exome'] = 2
out.loc[(out['cpd_exome'] >= 21) & (out['cpd_exome'] <= 30), 'cpd_exome'] = 3
out.loc[out['cpd_exome'] > 30, 'cpd_exome'] = 4

# SI - smoking initiation - 1 - Denies being a regular smoker,
# 2 - Regular smoker at some point
# Must have been a cigarette smoker

out['si'] = ukb_vars['f.2644.0.0'].replace(
    ['Yes', 'No', 'Do not know', 'Prefer not to answer'],
    ['2', '1', None, None])
out['si'] = out['si'].fillna(ukb_vars['f.2877.0.0'].replace(
    ['Hand-rolled cigarettes', 'Manufactured cigarettes',
        'Cigars or pipes', 'None of the above', 'Prefer not to answer'],
    ['2', '2', None, None, None]))
out['si'] = out['si'].fillna(ukb_vars['f.1249.0.0'].replace(
    ['Smoked occasionally', 'I have never smoked',
        'Smoked on most or all days', 'Just tried once or twice',
        'Prefer not to answer'],
    [None, '1', None, None, None]))

formerCig = ukb_vars['f.6183.0.0'].notnull()
formerCig.loc[formerCig == False] = None
formerCig = formerCig.astype(str)
formerCig = formerCig.replace(['nan', '1.0'], [None, '2'])
out['si'] = out['si'].fillna(formerCig)

out['si'] = out['si'].fillna(ukb_vars['f.3446.0.0'].replace(
    ['Manufactured cigarettes', 'Hand-rolled cigarettes', 'Cigars or pipes',
        'None of the above', 'Prefer not to answer'],
    ['2', '2', None, None, None]))

# SC - smoking cessation
# 2 - Current smoker, 1 - Former smoker

out['sc'] = ukb_vars['f.2644.0.0'].replace(
    ['Yes', 'No', 'Do not know', 'Prefer not to answer'],
    ['1', None, None, None])
out['sc'] = out['sc'].fillna(ukb_vars['f.2877.0.0'].replace(
    ['Hand-rolled cigarettes', 'Manufactured cigarettes', 'Cigars or pipes',
        'None of the above', 'Prefer not to answer'],
    ['1', '1', None, None, None]))
out['sc'] = out['sc'].fillna(ukb_vars['f.3446.0.0'].replace(
    ['Manufactured cigarettes', 'Hand-rolled cigarettes', 'Cigars or pipes',
        'None of the above', 'Prefer not to answer'],
    ['2', '2', None, None, None]))

# PY - pack years - qcpd divided by 20 to yield number of packs, multiplied by
# the period of regular smoking in years

age = ukb_vars['f.21003.0.0']
currSmokerAge = ukb_vars['f.3436.0.0']
prevSmokerAge = ukb_vars['f.2867.0.0']
stopSmokeAge = ukb_vars['f.2897.0.0']

currSmokerAge = currSmokerAge.replace([-3, -1], None)
prevSmokerAge = prevSmokerAge.replace([-3, -1], None)
stopSmokeAge = stopSmokeAge.replace([-3, -1], None)

smokerPeriod = age - currSmokerAge
smokerPeriod = smokerPeriod.fillna(stopSmokeAge - prevSmokerAge)
smokerPeriod.loc[smokerPeriod < 1] = None  # Some impossible intervals
out['py'] = (out['qcpd'] / 20) * smokerPeriod

# AI - age of initiation of regular smoking

out['ai'] = prevSmokerAge.fillna(currSmokerAge)
out.loc[out['ai'] < 10, 'ai'] = None
out.loc[out['ai'] > 35, 'ai'] = None

# DPW - alcoholic drinks per week

dpw = (ukb_vars['f.1568.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.1578.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.1588.0.0'].replace([-3, -1], 0) * 1.3) + \
    (ukb_vars['f.1598.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.1608.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.5364.0.0'].replace([-3, -1], 0).fillna(0) * 1.5)

dpm = (ukb_vars['f.4407.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.4418.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.4429.0.0'].replace([-3, -1], 0) * 1.3) + \
    (ukb_vars['f.4440.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.4451.0.0'].replace([-3, -1], 0)) + \
    (ukb_vars['f.4462.0.0'].replace([-3, -1], 0).fillna(0) * 1.5)

# Divide by the average number of days in a month times days in a week
dpw = dpw.fillna(dpm / 30.42 * 7)
dpw.loc[dpw > 168] = None  # Set > 24 drinks per day to missing
out['dpw'] = dpw

# Current or former smoker
# 1 - Formerly a regular smoker, 2 - Currently a regular smoker
out['currFormSmoker'] = ukb_vars['f.1239.0.0'].replace(
    ['No', 'Only occasionally', 'Yes, on most or all days',
        'Prefer not to answer'],
    [None, None, '2', None]).fillna(
    ukb_vars['f.1249.0.0'].replace(
        ['Smoked occasionally', 'I have never smoked',
            'Smoked on most or all days', 'Just tried once or twice',
            'Prefer not to answer'],
        [None, None, '1', None, None])).fillna(
    ukb_vars['f.2644.0.0'].replace(
        ['Yes', 'No', 'Do not know', 'Prefer not to answer'],
        ['1', None, None, None]))

# Current or former drinker - 1 - Formerly a drinker, 2 - Currently a drinker
out['currFormDrinker'] = ukb_vars['f.1558.0.0'].replace(
    ['Once or twice a week', 'Daily or almost daily', 'Never',
        'Three or four times a week', 'One to three times a month',
        'Special occasions only', 'Prefer not to answer'],
    ['2', '2', None, '2', '2', '2', None]).fillna(
    ukb_vars['f.3731.0.0'].replace(
        ['Yes', 'No', 'Prefer not to answer'],
        ['1', None, None]))

# Drinker/Non-Drinker - 1 - Not currently a drinker, 2 - Currently a drinker
out['dnd'] = ukb_vars['f.1558.0.0'].replace(
    ['Once or twice a week', 'Daily or almost daily', 'Never',
        'Three or four times a week', 'One to three times a month',
        'Special occasions only', 'Prefer not to answer'],
    ['2', '2', '1', '2', '2', '2', None])

out.to_csv('ukb_gscan_phen.csv', index=False)

