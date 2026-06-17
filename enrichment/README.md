# enrichment

GO-term enrichment of the DEG lists via DAVID and REVIGO: builds the gene/GO lists for DAVID,
visualizes DAVID functional annotation, and prepares REVIGO inputs.

Absorbed from Grace's repo. Scripts carry original paths; phase-4 rewrite pending.

## Layout

```
enrichment/
├── enrichment.Rproj
├── 01_code/
│   ├── 07-GOenrichment_listcreationforDAVID.Rmd
│   ├── 08-func_enrichment_DAVID_visual.Rmd
│   └── 09-GOenrichment_listcreation_REVIGO.Rmd
├── 02_data/                 no stored inputs; reads DEG lists from differential-expression
└── 03_analyses/
    ├── Func_annot_DAVID/    DAVID functional annotation results
    ├── Revigo_results/      REVIGO outputs
    └── uniprot_BG_DAVID.txt UniProt background gene list for DAVID
```

## Inputs (cross-folder)

These scripts read the DEG lists in `../differential-expression/03_analyses/DEG_lists/` and the
GO mapping in `../blast/03_analyses/genome-foot/`. In phase 4 those reads will use a `repo_root`
pointer. DAVID and REVIGO are web tools; the committed outputs are their results.

## Runnability (phase 4 pending)

Original paths, not yet repointed to this layout.
