# morphometrics

Shell and tissue measurements used to derive condition index (CI) and gonad index (GI).

```
morphometrics/
└── 02_data/
    └── morphometrics.xlsx   shell length/width/height and tissue dry masses
```

This folder is data-only. The CI and GI calculations and plots live in the cross-cutting
`summary-plots` analysis, which reads `morphometrics.xlsx` from here. There is no
`01_code/` because no script is specific to morphometrics.
