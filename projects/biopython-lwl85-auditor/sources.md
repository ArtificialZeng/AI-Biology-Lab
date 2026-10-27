# Primary sources and limits

- [Biopython source at the audited commit](https://github.com/biopython/biopython/blob/893f1543f87981927104fb7af9de43c4579dfe60/Bio/Align/analysis.py), plus its legacy `Bio/codonalign/codonseq.py` implementation and release tag `biopython-188`.
- [MEGA12 original Li–Wu–Luo method](https://www.megasoftware.net/web_help_12/Li-Wu-Luo_Method.htm), with the twofold/fourfold weighted formula.
- [PAML original LWL85 aggregation](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/yn00.c#L432); the same function separately implements LWL85m and LPB93.
- Li, Wu and Luo (1985), [DOI](https://doi.org/10.1093/oxfordjournals.molbev.a040343). Metadata and abstract inspected; full formula pages were inaccessible. Formula verification relies on the two named first-party implementations/documentation above.
- [NCBI standard genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1). The oracle contains the standard 64-entry TCAG-order amino-acid table and enumerates synonymous alternatives itself.

Current all-state keyword searches did not locate a competing weight correction. Existing issue #1611 concerns several dN/dS questions, including saturation; absence of a competing correction in bounded inspected sources is not proof of historical priority. Third-party source papers and private correspondence are not redistributed.
