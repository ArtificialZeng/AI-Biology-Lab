# Bounded claim and exact reference

For complete biallelic fixtures, positive explicit branch lengths, and finite nonconstant trait columns with positive exact sample variance, the pinned native baseline can lose the variation in a trait after adding a large common offset. Its variance formula subtracts the squared mean from the raw second moment. The candidate first subtracts each column's first sample, computes the mean of those differences, forms centered residuals, then sums their squares and scales by n-1. This is a conventional anchored two-pass sample-variance calculation shared by site, branch and node statistic modes.

The exact oracle uses Pearson squared correlation

`(n*sum(g*w)-sum(g)*sum(w))**2 / ((n*sum(g*g)-sum(g)**2)*(n*sum(w*w)-sum(w)**2))`.

All arithmetic forming expected values uses exact rational numbers. A constant allele indicator contributes zero. The trait columns have positive exact variance; the audit does not relax upstream's rejection of invalid or constant traits. For g=(0,0,1,1) and w=(0,1,2,3)+c, the exact site value is 4/5 for every tested offset c.

The fixed graph has sample node IDs 1 through n, deliberately excluding unused node 0 from the sample set. Explicit clades, two mutation positions, and known parent/child times yield exact expectations independently of native traversal. Site expectations average the two observed alleles; branch expectations sum branch length times the corresponding split correlation; node expectations retain every explicit node, including the unused node. Whole and split windows exercise span-normalized and raw outputs.

There are 396 API calls: n=2,4,6; zero offset plus both signs of 2^20, 2^26, 2^27, 2^30, and 2^40; two trait columns in every call; site/branch/node; whole/split windows; and both normalization options. Exact reference values are converted to binary floats only for comparison, with rtol=atol=2e-12. Thus the expected values are exact, while the API validation is tolerance-based floating-point evidence.

The packet establishes a concrete defect and its correction on these fixtures. It does not prove correctness for all networks, all finite doubles, missing genotypes, or extreme dynamic range. It does not demonstrate effects on published biological analyses. Independent AI review and script replay are computational checks, not human peer review.
