`PSMFC-mytilus-byssus-pilot/data`

---

- `psmfc_mussel_rna_summary.csv`: Summary table of RNA isolations. All isolations performed using [RNAzol RT](https://www.mrcgene.com/product/rnazol-rt) (MRC) with [Direct-zol RNA Miniprep Kit](https://www.zymoresearch.com/collections/direct-zol-rna-kits#:~:text=The%20Direct%2Dzol%20RNA%20Kits,TRI%20Reagent%C2%AE%20or%20similar.&text=These%20RNA%20purification%20kits%20allow,complete%20in%20only%207%20minutes.) (ZymoResearch). Contains the following columns:

  - `species`: Mussel species of source tissue.

  - `sample`: Sample name.

  - `concentration(ng/uL)`: RNA sample concentration in ng/uL.

  - `volume(uL)`: Total elution volume used, in uL. Volume should _not_ be considered to be accurate, as it does not account for subsequent usage after elution.

  - `yield(ng)`: Total RNA yield, in ng. Calculation is based off of elution volume. See note regarding the `volume(uL)` accuracy.

  - `tissue`: Source tissue type.

  - `isolation_date`: Date RNA was isolated.

  - `storage_location`: RNA sample storage location. All boxes are stored in the [Roberts Lab -80<sup>o</sup>C freezer](https://docs.google.com/spreadsheets/d/1Qsvz3QTURlPF_hX05BQxjom3484WuMfqQ1ILl9LEljU/edit?usp=sharing) (Google Sheet).

  - `quant_method`: Method used for quantifying RNA on date of isolation.

  - `DNased(y/n)`: Indicates if RNA was treated with DNase. `y` = Yes. `n` = No.