# Natural Killer cells demonstrate distinct eQTL and transcriptome-wide disease associations, highlighting their role in autoimmunity.

**Abstract**
Natural Killer (NK) cells are innate lymphocytes with central roles in immunosurveillance and are implicated in autoimmune pathogenesis. The degree to which regulatory variants affect NK gene expression is poorly understood. We performed expression quantitative trait locus (eQTL) mapping of negatively selected NK cells from a population of healthy Europeans (n=245). We find a significant subset of genes demonstrate eQTL specific to NK cells and these are highly informative of human disease, in particular autoimmunity. An NK cell transcriptome-wide association study (TWAS) across five common autoimmune diseases identified further novel associations at 27 genes. In addition to these cis observations, we find novel master-regulatory regions impacting expression of trans gene networks at regions including 19q13.4, the Killer cell Immunoglobulin-like Receptor (KIR) Region, GNLY and MC1R. Our findings provide new insights into the unique biology of NK cells, demonstrating markedly different eQTL from other immune cells, with implications for disease mechanisms.

[Preprint](https://www.biorxiv.org/content/10.1101/2021.05.10.443088v1)

**Overview of repository**
* Figure folders (Figure1-4 & Supp_Figures): contains R script and source data to reproduce each main and supplementary Figure from the manuscript.
* QTL_mapping: contains example code for the eQTL mapping pipeline in cis and trans.
* colocalisation_examples: contains example code and data for colocalisation analysis with coloc (data: CD226 eQTL and inflammatory bowel disease colocalisation) and moloc (data: CD226 eQTL colocalisation across 5 immune cell subsets).
* FUSION_TWAS: contains example code for computation of functional weights for NK eQTL mapping data () and TWAS analysis (). 

***eQTL mapping data***

[cis mapping results Shiny app](https://jjgilchrist.shinyapps.io/nk_cis_eqtl/)

[trans mapping results Shiny app](https://jjgilchrist.shinyapps.io/nk_trans_eqtl/)
