
# Brightway theory

Just for clarification, abbreviations used in the markdown cells of this notebook are:

| Abbreviation(s) | Written out |
| --- | --- |
| bw, BW2 | Brightway2 |
| FU | Functional unit |
| LCA | Life cycle assessment |
| LCI | Life cycle inventory |
| LCIA | Life cycle impact assessment |

---

## Recap on LCA

To make sure we are on the same page, let's recapitulate the basics of LCA, expressed in a single formula (the notation may differ to the one that you are familiar with):

$$h = CBA^{-1}f$$

where:

- $A$... technosphere matrix (dimensions *p x p*)
- $B$... biosphere matrix (dimensions *q x p*)
- $C$... characterisation matrix (dimensions *r x q*)
- $f$... final demand vector (dimensions *p x 1*)
- $h$... characterised inventory matrix (dimensions *r x 1*)

with the dimensions being:

- *p*... number of products/processes
- *q*... number of elementary flows
- *r*... number of impact categories

Specific parts of the above formula have distinct names. Since these names are also used in the Brightway2-framework, we will note them down here:

- $f$... demand array
- $A{^-1}f$... supply array
- $BA{^-1}f$... inventory referred sometimes as b
- $CBA{^-1}f$... characterised inventory
