import pandas as pd

# Load the filtered VCF into a dataframe (you may need to use vcf parser like pyVCF or manually handle the fields)
vcf_file = "priority_snps.vcf"
df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)

# Define a scoring function based on criteria
def score_variant(row):
    score = 0
    if row['IMPACT'] in ['HIGH', 'MODERATE']:
        score += 2 if row['IMPACT'] == 'HIGH' else 1
    if row['CADD_PHRED'] >= 20:
        score += 1
    if row['ClinVar'] in ['Pathogenic', 'Likely_pathogenic']:
        score += 2
    if row['gnomAD_AF'] < 0.01:
        score += 1
    return score

# Apply the scoring function
df['score'] = df.apply(score_variant, axis=1)

# Sort variants by score
df_sorted = df.sort_values(by='score', ascending=False)

# Output the ranked variants to CSV
df_sorted.to_csv("ranked_snps.csv", index=False)
