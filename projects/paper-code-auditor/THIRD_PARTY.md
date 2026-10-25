# Third-party materials

`fixtures/upstream/` contains selected files from [Luiz Vieira's DMS-SiteEffect-PLM](https://github.com/ziul-bio/DMS-SiteEffect-PLM) at commit `e168a2c825be0c3e28c300b8f6277dda1d999c4d`. Its MIT copyright and permission notice is retained as `fixtures/upstream/LICENSE.txt`.

`fixtures/historical/` contains two original files from the same repository at commit `e0c2d6ae06903cb356abb7e1b579058115e5623a`, recovered from public Git history. The restoration patch reproduces those files with a documented logical-device adjustment in the trainer; the model implementation is unchanged. These files are included for inspectable regression testing and source restoration, not as newly authored models.

`fixtures/data/IF1_ECOLI_identifiers.csv` is a two-column projection (`ID`, `mutant`) of the upstream `data/cellular/metadata/IF1_ECOLI.csv`. It retains row order and identifiers for checking split behavior and contains neither sequences nor experimental phenotype values. Original assay and derived fixture hashes are recorded in `evidence/`.

Immutable file URLs, SHA-256 hashes, the complete pinned tree inventory, and restoration metadata are supplied in `evidence/provenance.json`, `evidence/upstream_inventory.json`, and `evidence/restoration_provenance.json`. No paper PDFs, model weights, private correspondence, or credentials are distributed.
