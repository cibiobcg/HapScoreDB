# Filippo Gastaldello, Fabio Mazza, Alessandro Romanel - 13/05/25

library(shiny)
library(bslib)
library(tidyverse)
library(DT)
library(hrbrthemes)
library(plotly)
library(viridis)
library(ggpubr)
library(scales)
library(patchwork)
library(shinycssloaders)

# DATA PATHS
interface_data_path <- "../data/interface_data.tsv"
rsid_map_path <- "../data/rsid_map.tsv"
coords_map_path <- "../data/coords_map.tsv"
sequences_path <- "../data/all_sequences.tsv"
version_df_path <- "../data/version.csv"
api_df_path <- "../data/api.csv"
link_map_path<- "../data/link_map.rds"

# LOAD DATA ----
data <- read_tsv(interface_data_path)
rsid_map <- read_tsv(rsid_map_path)
rsid_map <- rsid_map %>% column_to_rownames("rsid") %>% mutate(haplo_ids = strsplit(haplo_ids, ","))
coords_map <- read_tsv(coords_map_path)
coords_map <- coords_map %>% column_to_rownames("variant_coord") %>% mutate(haplo_ids = strsplit(haplo_ids, ","))
version <- read_csv(version_df_path)
api <- read_csv(api_df_path)
link_map <- readRDS(link_map_path)
link_map <- link_map %>% column_to_rownames('rsid')
# Footer definition ----
app_footer <- tags$footer(
    class = "text-center text-white",
    style = "background-color: #004c5a; padding-top: 20px; padding-bottom: 15px;",
    div(class = "container-fluid",
        h5("CONTACTS", style = "font-weight: bold; text-transform: uppercase; letter-spacing: 1px;"),
        hr(style = "border-top: 1px solid rgba(255, 255, 255, 0.3); max-width: 60%; margin: 1.2rem auto;"),
        p(HTML(paste("&copy;", format(Sys.Date(), "%Y"), "Copyright: University of Trento"))),
        p("Department of Cellular, Computational and Integrative Biology CIBIO"),
        p("Laboratory of Bioinformatics and Computational Genomics", style = "margin-bottom: 0;")
    )
)

# Value box theme definition ----
valuebox_theme <- value_box_theme(
    bg = "#ffffff",
    fg = "#e39107"
)

# UI ----
ui <- tagList(
  
    withMathJax(), # enable mathjax
  
    page_navbar(
        ## CSS styling ----
        
        # Accordion title bold style
        tags$head(
            tags$style(HTML(".accordion-button {
                            font-weight: bold;
                            }")
                       )
        ),
        title = tagList("HapScoreDB", tags$sub("v1.0")),
        theme = bs_theme(bootswatch = "journal"),
        # Prevent the theme from changing title to uppercase
        tags$head(
            tags$style(HTML("
                            .navbar-brand {
                                text-transform: none !important;
                            }
                            ")
                       )
        ),
        nav_spacer(),
        ## Search ----
        nav_panel(
            title = "Search",
            layout_sidebar(
                fillable = FALSE,
                ### Sidebar ----
                sidebar = sidebar(
                            width = "20%",
                            h3("Search options"),
                            #### Search type selector ----
                            radioButtons(
                                inputId = "search_type",
                                label = h5("Search by"),
                                choices = list(
                                    "Ensembl Gene ID" = "ensg",
                                    "Gene Symbol" = "hgnc",
                                    "Ensembl Transcript ID" = "transcript",
                                    "Variant Coordinates" = "coordinates",
                                    "Variant rsid" = "rsid"
                                ),
                                selected = "hgnc"
                            ),
                            #### ID text input ----
                            textInput(
                                inputId = "search_id",
                                label = "ID",
                                placeholder = "Insert ID here"
                            ),
                            accordion(
                                accordion_panel(
                                    title = "Advanced search",
                                    #### Variant type filter ----
                                    checkboxGroupInput(
                                        inputId = "variant_type_filter",
                                        label = h5("Exclude haplotypes with:"),
                                        choices = list("Synonymous variants" = "synonymous_variant",
                                                       "Missense variants" = "missense_variant",
                                                       "Nonsense variants" = "stop_gained",
                                                       "Frameshift variants" = "framesift_variants",
                                                       "In frame deletions" = "inframe_deletion",
                                                       "In frame insertions" = "inframe_insertion",
                                                       "Start loss" = "start_lost",
                                                       "Start gain" = "initiator_codon_variant",
                                                       "Stop loss" = "stop_lost"
                                                       )
                                    ),
                                    #### Ensembl tsl filter ----
                                    sliderInput(
                                        inputId = "tsl_filter",
                                        label = h5("Select desired Ensembl transcript support level (tsl):"),
                                        min = 1,
                                        max = 5,
                                        value = 5
                                    ),
                                    input_switch("tsl_NA",
                                                 "Remove NAs",
                                                 value = FALSE
                                    ),
                                    input_switch("tsl_end_NF",
                                                 "Remove END_NFs",
                                                 value = FALSE
                                    )
                                ),
                                id = "advanced_search",
                                open = FALSE
                            )
                    ),
                    ### Main content ----
                    
                    #### Value box ----
                    layout_column_wrap(
                        value_box( 
                            title = "Selected haplotypes",
                            showcase = icon("user", "fa-regular"),
                            textOutput(outputId = "selected_haplotypes"), 
                            theme = valuebox_theme
                        ), 
                        value_box( 
                            title = "Affected transcripts",
                            showcase = icon("dna", "fa-regular"), 
                            textOutput(outputId = "affected_transcripts"), 
                            theme = valuebox_theme 
                        ), 
                        value_box( 
                            title = "Involved variants",
                            showcase = icon("disease", "fa-regular"),
                            textOutput(outputId = "involved_variants"), 
                            theme = valuebox_theme
                        ),
                        max_height = "10%"
                    ),
                    
                    #### datatable ----
                    card(
                        card_header("Haplotypes table"),
                        DTOutput("datatable"),
                        full_screen = TRUE
                    ),
                    
                    layout_column_wrap(
                        #### Score distributions ----
                        card(
                            card_header("Score distribution"),
                            layout_sidebar(
                                sidebar = sidebar(
                                    selectizeInput(
                                        inputId = "model",
                                        label = "Select the desired model:",
                                        choices = list("ESM2"="esm",
                                                       "PoET"="poet",
                                                       "proSST"="prosst"
                                                       )
                                    ),
                                    selectizeInput(
                                            inputId = "score_type",
                                            label = "Select the desired score type:",
                                            choices = list("PLL" = "_PLL",
                                                           "PLLR_wt"="_PLLR_wt",
                                                           "PLLR_mf"="_PLLR_mf"
                                            )
                                    ),
                                    selectizeInput(
                                        inputId = "score_distribution_filter",
                                        label = "Restrict to:",
                                        choice = NULL,
                                        options = list(placeholder = "Select an option")
                                    )
                                ),
                                withSpinner(plotlyOutput(outputId = "score_distribution"), type = 4, color = "#e39107")
                            ),
                            full_screen = TRUE,
                            fill = TRUE
                        ),
                        
                        #### Score deltas distributions ----
                        card(
                            card_header("Score deltas distribution"),
                            layout_sidebar(
                                sidebar = sidebar(
                                    selectizeInput(
                                        inputId = "delta",
                                        label = "Select the desired model's delta:",
                                        choices = list("ESM2 Delta" = "esm_PLL_delta",
                                                       "PoET Delta" = "poet_PLL_delta",
                                                       "proSST Delta" = "prosst_PLL_delta"
                                                       )
                                    ),
                                    selectizeInput(
                                        inputId = "delta_distribution_filter",
                                        label = "Restrict to:",
                                        choice = NULL,
                                        options = list(placeholder = "Select an option")
                                    )
                                ),
                                withSpinner(plotlyOutput(outputId = "score_delta_distribution"), type = 4, color = "#e49107")
                            ),
                            full_screen = TRUE,
                            fill = TRUE
                        )
                    ),
                    #### Ancestry ----
                    card(
                        card_header("Ancestry"),
                        layout_sidebar(
                            sidebar = sidebar(
                                radioButtons(
                                    inputId = "group_by",
                                    label = "Group by",
                                    choices = list(
                                        "Haplotype" = "haplotype",
                                        "Ancestry" = "ancestry"
                                    ),
                                    selected = "haplotype"
                                ),
                                selectizeInput(
                                    inputId = "ancestry_frequency_filter",
                                    label = "Restrict to:",
                                    choice = NULL,
                                    options = list(placeholder = "Select an option")
                                )
                            ),
                            withSpinner(plotlyOutput("population"), type = 4, color = "#e39107")
                        ),
                        full_screen = TRUE
                    ),
                    ### Models correlation selected scores ----
                    card(
                        card_header("Model scores correlation"),
                        layout_sidebar(
                            sidebar = sidebar(
                                selectizeInput(
                                    inputId = "model1_selected",
                                    label = "Select the desired model:",
                                    choices = list("ESM2" = "esm",
                                                   "PoET" = "poet",
                                                   "proSST" = "prosst")
                                ),
                                selectizeInput(
                                    inputId = "model2_selected",
                                    label = "Select the desired model:",
                                    choices = list("PoET" = "poet",
                                                   "ESM2" = "esm",
                                                   "proSST" = "prosst")
                                )
                            ),
                            withSpinner(plotOutput(outputId = "models_correlation_selected_score"), type = 4, color = "#e39107")
                        ),
                        full_screen = TRUE,
                        fill = TRUE 
                    )
                )
            ),
            ## About ----
            nav_panel(
                title = "About",
                page_fluid(
                    ### Value box ----
                    layout_column_wrap(
                        value_box( 
                            title = "Genes",
                            textOutput(outputId = "analyzed_genes"), 
                            theme = valuebox_theme
                        ), 
                        value_box( 
                            title = "Transcripts",
                            textOutput(outputId = "analyzed_transcripts"), 
                            theme = valuebox_theme
                        ), 
                        value_box(
                            title = "Haplotypes",
                            textOutput(outputId = "analyzed_haplotypes"), 
                            theme = valuebox_theme
                        ),
                        value_box( 
                            title = "Variants",
                            textOutput(outputId = "analyzed_variants"), 
                            theme = valuebox_theme
                        ),
                        max_height = "10%"
                    ),
                    layout_column_wrap(
                        ### Haplotypes per gene distributions ----
                        card(
                            card_header("Haplotypes per gene distribution"),
                            withSpinner(plotlyOutput(outputId = "haplotype_gene_distribution"), type = 4, color = "#e39107"),
                            full_screen = TRUE,
                            fill = TRUE
                        ),
                        ### Variants per haplotype distributions ----
                        card(
                            card_header("Variants per haplotype distribution"),
                            withSpinner(plotlyOutput(outputId = "variant_haplotype_distribution"), type = 4, color = "#e39107"),
                            full_screen = TRUE,
                            fill = TRUE
                        ),
                        max_height = "45%"
                    ),
                    ### Models correlation all scores ----
                    card(
                        card_header("Model scores correlation"),
                        layout_sidebar(
                            sidebar = sidebar(
                                selectizeInput(
                                    inputId = "model1_all",
                                    label = "Select the desired model:",
                                    choices = list("ESM2" = "esm",
                                                   "PoET" = "poet",
                                                   "proSST" = "prosst")
                                ),
                                selectizeInput(
                                    inputId = "model2_all",
                                    label = "Select the desired model:",
                                    choices = list("PoET" = "poet",
                                                   "ESM2" = "esm",
                                                   "proSST" = "prosst")
                                )
                            ),
                            withSpinner(plotOutput(outputId = "models_correlation_all_score"), type = 4, color = "#e39107")
                        ),
                        full_screen = TRUE,
                        fill = TRUE 
                    )
                )
            ),
            ## Download ----
            nav_panel(
                title = "Download",
                layout_sidebar(
                    fillable = FALSE,
                    ### Filters ----
                    sidebar = sidebar(
                        width = "20%",
                        h4("Select fields to download"),
                        accordion(
                            accordion_panel(
                                title = "IDs",
                                checkboxGroupInput(
                                    inputId = "IDs_selector",
                                    label = "",
                                    choices = list(
                                        "Ensembl Gene ID" = "gene_id",
                                        "Ensembl Transcript ID" = "transcript_id",
                                        "Haplotype ID" = "haplotype_id",
                                        "Gene Symbol" = "gene_symbol",
                                        "Uniprot ID" = "uniprot_id",
                                        "Variant rsid" = "rsid"
                                    ),
                                    selected = c("gene_id", "transcript_id", "haplotype_id")
                                )
                            ),
                            accordion_panel(
                                title = "Variant info",
                                checkboxGroupInput(
                                    inputId = "variant_info_selector",
                                    label = "",
                                    choices = list(
                                        "Variant coordinates" = "variant_coordinates_hg38",
                                        "Reference allele" = "ref_alleles",
                                        "Alternative allele" = "alt_alleles",
                                        "DNA changes" = "dna_changes",
                                        "Protein changes" = "protein_changes",
                                        "Variant type" = "variant_types"
                                    ),
                                    selected = "variant_coordinates_hg38"
                                )
                            ),
                            accordion_panel(
                                title = "Frequencies",
                                checkboxGroupInput(
                                    inputId = "frequency_selector",
                                    label = "",
                                    choices = list(
                                        "Global" = "frequency",
                                        "AFR" = "AFR_freq",
                                        "AMR" = "AMR_freq",
                                        "EAS" = "EAS_freq",
                                        "EUR" = "EUR_freq",
                                        "SAS" = "SAS_freq"
                                    ),
                                    selected = "frequency"
                                )
                            ),
                            accordion_panel(
                                title = "Transcript info",
                                checkboxGroupInput(
                                    inputId = "transcript_info_selector",
                                    label = "",
                                    choices = list(
                                        "Protein length" = "protein_length",
                                        "Number of variants" = "n_variants",
                                        "TSL - Ensembl Transcript Support Level" = "tsl"
                                    ),
                                    selected = "tsl"
                                )
                            ),
                            accordion_panel(
                                title = "Sequences",
                                checkboxGroupInput(
                                    inputId = "sequence_selector",
                                    label = "",
                                    choices = list(
                                        "DNA sequence" = "dna_sequence",
                                        "Protein sequence" = "protein_sequence"
                                    )
                                )
                            ),
                            accordion_panel(
                                title = "Scores",
                                checkboxGroupInput(
                                    inputId = "scores_selector",
                                    label = "",
                                    choices = list(
                                        "ESM2 PLL" = "esm_PLL",
                                        "ESM2 PLLR_mf" = "esm_PLLR_mf",
                                        "ESM2 PLLR_wt" = "esm_PLLR_wt",
                                        "ESM2 delta PLL" = "esm_PLL_delta",
                                        "PoET PLL" = "poet_PLL",
                                        "PoET PLLR_mf" = "poet_PLLR_mf",
                                        "PoET PLLR_wt" = "poet_PLLR_wt",
                                        "PoET delta PLL" = "poet_PLL_delta",
                                        "proSST PLL" = "prosst_PLL",
                                        "proSST PLLR_mf" = "prosst_PLLR_mf",
                                        "proSST PLLR_wt" = "prosst_PLLR_wt",
                                        "proSST delta PLL" = "prosst_PLL_delta",
                                        "AlphaMissense classes" = "am_classes",
                                        "AlphaMissense scores" = "am_scores",
                                        "AlphaMissense average score" = "am_avg_score",
                                        "AlphaMissense sum score" = "am_sum_score"
                                    ),
                                    selected = "esm_PLL"
                                )
                            ),
                            open = FALSE
                        )
                    ),
                    card(
                        card_header("Preview table"),
                        dataTableOutput(outputId = "preview_table"),
                        full_screen = TRUE
                    ),
                    downloadButton(
                        outputId = "download_button",
                        label = "Download")
                )
            ),
            ## FAQ ----
            nav_panel(
                title = "FAQ",
                h1("Any questions?", style = "text-align:center"),
                fluidRow(
                    column(
                        width = 8,
                        offset = 2,
                        accordion(
                            id = "FAQ",
                            accordion_panel(
                                title = "What is HapScoreDB?",
                                "HapScoreDB is a proteogenomic database providing a comprehensive collection of Protein Language
                                Model (PLM) scores for haplotype-resolved protein-coding sequences, encompassing all human
                                transcript isoforms. Our approach uniquely integrates GENCODE and Ensembl gene and transcript
                                models with phased variant data from the 1000 Genomes Project. For each protein-coding transcript,
                                we have meticulously reconstructed sequences of common protein haplotypes containing single or
                                multiple variants, including both SNPs (single nucleotide polymorphisms) and INDELs
                                (insertions/deletions). To quantify functional impact, we computed scores using state-of-the-art protein
                                language models such as ESM2, ProSST, and PoET, providing deep and contextualized
                                representations for each protein haplotype."
                            ),
                            accordion_panel(
                              title = "What is a Protein Language Model?",
                              "A protein language model is a class of probabilistic, computational
                                models engineered to capture the complex statistical patterns inherent
                                in protein sequences. The core purpose of these models is to
                                learn a high-dimensional representation of the sequence space
                                that corresponds to functional, evolutionarily viable proteins.
                                To achieve this, these models are trained on extensive biological
                                sequence databases with the same techniques used in Natural Language Processing. 
                                A predominant self-supervised strategy is masked language modeling, where the model learns to
                                predict masked amino acids based on their surrounding context.
                                Through this and similar tasks, the model implicitly learns the complex patterns 
                                underlying protein structure and function, including the co-evolutionary
                                constraints between residue positions—information that is otherwise
                                explicitly derived from the analysis of a Multiple Sequence Alignment (MSA).
                                Once trained, a key downstream capability of these models is to evaluate an arbitrary
                                protein sequence by assigning it a quantitative score. We compute this score as the 
                                pseudo-log-likelihood of the sequences, 
                                representing the average of conditional probabilities of each amino acid given
                                the context of the entire sequence. A higher pseudo-log-likelihood
                                score indicates that the sequence has a higher probability under the
                                learned statistical distribution, thereby providing a robust,
                                quantitative measure of its biological plausibility and consistency
                                with known evolutionary and structural constraints."
                            ),
                            accordion_panel(
                              title = "How are the scores calculated?",
                              tags$p("For ESM2 and ProSST, we compute the predicted pseudo-log-likelihoods of every protein sequence passed once to the models, without masking any amino acid, following Brandes et al 2023 (doi.org/10.1038/s41588-023-01465-0):"),
                              HTML("$$ P_{PLL}(S) = \\frac{1}{L} \\sum_{i=1}^{L} \\log(P(r_i = r_i^S | S_{mt})) $$"),
                              tags$p("Where $S$ is the input protein sequence, $L$ its length and $r_i^S$ the aminoacid at position $i$ in $S$. We find this approach for masked language models to correlate well with experimental data from the ProteinGym benchmark in the zero-shot setting. Albeit at a slighly lower performance for single variant effect prediction compared to the commonly used wild-type marginals method, this type of score can capture the effect of indels and of epistatic interactions, when multiple variants are present in the same sequence."),
                              tags$p("In the case of ProSST, the input structures are taken from the AlphaFoldDB (Varadi et al 2024, doi.org/10.1093/nar/gkad1011), and encoded using the 4096 long version of the structure sequence alphabet. For PoET, an autoregressive generative model, we use the same scoring function as in the original paper, averaging over multiple context lengths as in their work, and using as inputs MSAs generated from the Uniref100 database employing the ColabFold protocol.")                              
                            ),
                            accordion_panel(
                                title = "What do PLLR_wt and PLLR_mf mean?",
                                "The Pseudo Log-Likelihood (PLL) of each haplotype-transcript
                                couple should be compared to a reference to quantify how much a
                                certain haplotype can shift the score of the transcript. In
                                order to do so we adopted the Pseudo Log-Likelihood Ratio (PLLR) metric
                                already known in literature. Since we are handling log-likelihoods,
                                the ratio is the difference between two PLLs. We computed, for
                                each score, its ratio with respect to the score of the most
                                frequent haplotype for the specific transcript, namely PLLR-mf,
                                and with respect to what is annotated as wild type in Ensembl,
                                PLLR-wt. Most of the time these two ratios coincide."
                            ),
                            accordion_panel(
                                title = "What do score deltas represent?",
                                "Score deltas are transcript level metrics, meaning that all
                                entries in the table sharing the same transcript ID will share
                                the same score delta (one per model). The score delta is computed
                                as the difference between the maximum and the minimum score
                                obtained by any haplotype of that transcript. Comparing the
                                score delta of a transcript against the distribution of all
                                score deltas should give a quantitative measure of the functional
                                variability of a transcript in the population."
                            ),
                            accordion_panel(
                                title = "How were protein sequences generated?",
                                "For each transcript’s isoform, all haplotypes comprising 
                                common variants were identified. For each of these haplotypes,
                                a mutated sequence was generated by including the haplotype’s
                                variants in the wild type coding sequence. When a Transcription
                                Start Site is lost due to variants, the first next starting codon
                                is used as Transcription Start Site. When a nonsense mutation
                                occurs, the protein is truncated at that point."
                            ),
                            accordion_panel(
                                title = "What types of genes and transcript are taken into consideration?",
                                "Genes and transcripts annotations are downloaded from Ensembl.
                                In total, 18741 protein coding genes and all their isoforms leading
                                to 78282 transcripts were downloaded. Transcripts longer than 4000 amino acids were
                                removed due to model and computational limitations, while transcripts
                                shorter than 10 amino acids were removed due to "
                            ),
                            accordion_panel(
                                title = "Where do genotype data come from?",
                                "Genotypes data are taken from the phased 1000 Genomes release phase 3, aligned to reference genome hg38.
                                Lowy-Gallego E, Fairley S, Zheng-Bradley X et al. Variant calling on the GRCh38 assembly with the data 
                                from phase three of the 1000 Genomes Project . Wellcome Open Res 2019, 4:50 
                                (https://doi.org/10.12688/wellcomeopenres.15126.2)"
                            ),
                            accordion_panel(
                                title = "Which variants were taken into consideration?",
                                "To generate the database of mutated sequences all (SNPs and INDELs)
                                phased common (MAF > 0.5%) coding variants identified in 1000 Genomes were taken into consideration."
                            ),
                            open = FALSE
                        )
                    )
                )
            ),
            ## Access ----
            nav_panel(
                    title = "Access",
                    fluidRow(
                            column(
                                    width = 6,
                                    offset = 3,
                                    accordion(
                                            id = "ACCESS",
                                            accordion_panel(
                                                    title = "Version",
                                                    tableOutput("version")
                                            ),
                                            accordion_panel(
                                                    title = "License",
                                                    "All data and download files in HapScoreDB are freely available under a 'Creative Commons BY 4.0'
                                                    license. When using the data, please provide appropriate credit — and inform users of any
                                                    changes or additions that you might have made to the data.",
                                                    tags$br(),
                                                    tags$img(
                                                            src = "https://licensebuttons.net/l/by/4.0/88x31.png",
                                                            alt = "Creative Commons BY 4.0",
                                                            style = "margin-top:10px;"
                                                    )
                                            ),
                                            accordion_panel(
                                                    title = "API",
                                                    "HapScoreDB data are also available through rest API at the endpoint:",
                                                    tags$br(),
                                                    tags$br(),
                                                    h3("https://bcglab.cibio.unitn.it/hapscoredbAPI/data"),
                                                    tags$br(),
                                                    "with the following parameters:",
                                                    tags$br(),
                                                    tags$br(),
                                                    tableOutput("api"),
                                                    tags$br(),
                                                    "Example: ", tags$code("curl \"https://bcglab.cibio.unitn.it/hapscoredbAPI/data?genes=ENSG00000164002,ENSG00000065978\"")
                                            ),
                                            open = FALSE
                                    )
                            )
                    )
            ),
            ## Contacts ----
            nav_item(
                tags$a(
                    "Contacts",
                    href = "https://www.cibio.unitn.it/785/laboratory-of-bioinformatics-and-computational-genomics",
                    target = "_blank"
                )    
            ),
            nav_spacer()
    )
)

# SERVER LOGIC ----
server <- function(input, output, session){
    ## Search ----
    ### Update ID text input and search bar example ----
    # Observe changes in the radioButtons
    observeEvent(input$search_type, {
        selected_type <- input$search_type
        new_label <- ""
        
        if (selected_type == "ensg") {
            new_label <- "Ensembl Gene ID (e.g.):"
            example_search <- "ENSG00000164002"
        } else if (selected_type == "hgnc") {
            new_label <- "HGNC gene symbol (e.g.):"
            example_search <- "EXO5"
        }else if (selected_type == "transcript") {
            new_label <- "Ensembl Transcript ID (e.g.):"
            example_search <- "ENST00000682383"
        } else if (selected_type == "coordinates") {
            new_label <- "Variant coordinates (chr:pos.REF>ALT, e.g.):"
            example_search <- "chr1:40515059.G>T"
        }else if (selected_type == "rsid") {
            new_label <- "Variant rsid:"
            example_search <- "rs11208299" 
        }
        
        # Update the label of the textInput
        updateTextInput(
            inputId = "search_id",
            value = example_search,
            label = new_label
        )
    })
    
    ### Subset search dataframe ---- 
    df <- reactive({
        if (input$search_type == "ensg") {
            df <- data %>% filter(gene_id == input$search_id)
        }
        if (input$search_type == "hgnc") {
            df <- data %>% filter(gene_symbol == input$search_id)
        }
        if (input$search_type == "transcript") {
            df <- data %>% filter(transcript_id == input$search_id)
        }
        if (input$search_type == "coordinates") {
            haplotypes <- coords_map[input$search_id,"haplo_ids"][[1]]
            df <- data %>% filter(haplotype_id %in% haplotypes)
        }
        if (input$search_type == "rsid") {
            haplotypes <- rsid_map[input$search_id,"haplo_ids"][[1]]
            df <- data %>% filter(haplotype_id %in% haplotypes)
        }
        # Check advanced filters
        if (length(input$variant_type_filter)>0) {
            try(
                df <- df %>% filter(!grepl(paste0(input$variant_type_filter, collapse = "|"), df$variant_types))
            )
        }
        # Tsl filter is always set to 5 by default, no need to check if it is set
        try(
            df <- df %>% filter(tsl <= input$tsl_filter | is.na(tsl))
        )
        if (input$tsl_NA) {
            df <- df %>% filter(!grepl("NA",tx_notes))
        }
        if (input$tsl_end_NF) {
            df <- df %>% filter(!grepl("end_NF",tx_notes))
        }
        
        if (length(unique(df$gene_id))>1) {
                showModal(
                    modalDialog(
                        title = "Multiple genes",
                        "The selected variant affects more than one gene, to focus on one gene only,
                        copy its ID from the table and paste it into a new search!",
                        easyClose = TRUE
                    )
                )
        }
        
        # Create variants SNPs to clinvar
        df$rsid <- lapply(df$rsid, function(x){
                if (x=="wt") {
                        return("wt")
                }else{
                        link <- link_map[x,1][[1]]
                        return(ifelse(!is.null(link), link, sprintf('<a href="https://www.ncbi.nlm.nih.gov/snp/%s" target="_blank">%s</a>', x, x)))
                }
        })
        
        # Create links to uniprot
        df$uniprot_id <- lapply(df$uniprot_id, function(x){
                if (is.na(x)) {
                        return(x)
                }else{
                        return(sprintf('<a href="https://www.uniprot.org/uniprotkb/%s" target="_blank">%s</a>', x, x))
                }
        })
        
        return(df)
        
    })
    
    ### Update value boxes ----
    output$selected_haplotypes <- renderText({
        nrow(df() %>% select(haplotype_id) %>% unique())
    })
    output$affected_transcripts <- renderText({
        nrow(df() %>% select(transcript_id) %>% unique())
    })
        
    output$involved_variants <- renderText({
        subset <- df() %>% select(variant_coordinates_hg38) %>% filter(!variant_coordinates_hg38 == "wt")
        variants <- c()
        for(row in subset$variant_coordinates_hg38){
            variants <- c(variants, str_split_1(row, pattern = ","))
        }
        length(unique(variants))
    })
    
    ### Search datatable ---- 
    output$datatable <- DT::renderDT({
        df <- df()
        DT::datatable(df, escape=FALSE, options = list(pageLength = 5))
    })
    
    ### Update distribution filters, both for scores and deltas----
    observe({
        # Extract haplotype and trancsript ID from subset df 
        current_df_val <- df()
        haplotypes <- c("none", sort(unique(current_df_val$haplotype_id)))
        transcripts <- c("none", sort(unique(current_df_val$transcript_id)))
        updateSelectizeInput(session,
                             inputId = "score_distribution_filter",
                             label = "Restrict to:",
                             choices = list(
                                        "Haplotypes" = as.list(haplotypes),
                                        "transcripts" = as.list(transcripts)
                             ),
                             server = FALSE)
        updateSelectizeInput(session,
                             inputId = "delta_distribution_filter",
                             label = "Restrict to:",
                             choices = list(
                                    "Haplotypes" = as.list(haplotypes),
                                    "transcripts" = as.list(transcripts)
                             ),
                             server = FALSE)
        
    })
    
    ### Score distribution plots ----
    output$score_distribution <- renderPlotly({
        data_plot <- data
        if (nrow(df()>0)) {
            plot_df <- df()
            # Check if filtering options have been applied for the plot
            if (str_detect(input$score_distribution_filter, "ENST")) {
                plot_df <- plot_df %>% filter(transcript_id == input$score_distribution_filter)
            }
            if (str_detect(input$score_distribution_filter, "\\.")) {
                plot_df <- plot_df %>% filter(haplotype_id == input$score_distribution_filter)
            }
            # Check if advanced filtering was applied to df(), in that case, also apply it 
            # to the background distribution
            if (length(input$variant_type_filter)>0) {
                try(
                    data_plot <- data_plot %>% filter(!grepl(paste0(input$variant_type_filter, collapse = "|"), data_plot$variant_types))
                )
            }
            try(
                data_plot <- data_plot %>% filter(tsl <= input$tsl_filter | is.na(tsl))
            )
            if (input$tsl_NA) {
                data_plot <- data_plot %>% filter(!grepl("NA",tx_notes))
            }
            if (input$tsl_end_NF) {
                data_plot <- data_plot %>% filter(!grepl("end_NF",tx_notes))
            }
            if (str_detect(paste0(input$model, input$score_type), "PLLR")) {
                p <- data_plot %>% filter(!paste0(input$model, input$score_type)==0) %>% ggplot() +
                                                            geom_density(
                                                                aes(x = !!sym(paste0(input$model, input$score_type))),
                                                                fill = "#e6ab47",
                                                                color = "#e39107",
                                                                alpha = 0.8,
                                                                adjust = 0.5,
                                                                stat = "count") +
                                                            geom_vline(
                                                                data = plot_df, 
                                                                aes(xintercept = !!sym(paste0(input$model, input$score_type)),
                                                                    text = paste0(haplotype_id,
                                                                                  " on transcript ",
                                                                                  transcript_id,
                                                                                  "\n",
                                                                                  input$model,
                                                                                  ": ",
                                                                                  round(!!sym(paste0(input$model, input$score_type)),4)
                                                                                  )
                                                                    ),
                                                                color = "#e34907") +
                                                            theme_minimal() +
                                                            theme(legend.position = "none") +
                                                            scale_y_log10(name = "Density (log10 scale)",
                                                                          breaks = trans_breaks("log10", function(x) 10^x),
                                                                          labels = trans_format("log10", math_format(10^.x)
                                                                                                )
                                                                          )
            }else{
                p <- data_plot %>% ggplot() +
                                geom_density(
                                    aes(x = !!sym(paste0(input$model, input$score_type))),
                                    fill = "#e6ab47",
                                    color = "#e39107",
                                    alpha = 0.8,
                                    adjust = 0.1) +
                                geom_vline(
                                    data = plot_df, 
                                    aes(xintercept = !!sym(paste0(input$model, input$score_type)),
                                        text = paste0(haplotype_id,
                                                      " on transcript ",
                                                      transcript_id,
                                                      "\n",
                                                      input$model,
                                                      ": ",
                                                      round(!!sym(paste0(input$model, input$score_type)),4)
                                        )
                                    ),
                                    color = "#e34907") +
                                theme_minimal() +
                                theme(legend.position = "none")
            }
            ggplotly(p, tooltip = "text")
        }else{
            if (str_detect(paste0(input$model, input$score_type), "PLLR")) {
                p <- data_plot %>% filter(!paste0(input$model, input$score_type)==0) %>% ggplot(aes(x = !!sym(paste0(input$model, input$score_type)))) +
                                                            geom_density(fill = "#e6ab47",
                                                                         color = "#e39107",
                                                                         alpha = 0.8,
                                                                         adjust = 0.5,
                                                                         stat = "count"
                                                            ) +
                                                            theme_minimal() +
                                                            scale_y_log10(name = "Density (log10 scale)",
                                                                          breaks = trans_breaks("log10", function(x) 10^x),
                                                                          labels = trans_format("log10", math_format(10^.x)
                                                                                                )
                                                                          )
            }else{
                p <- data_plot %>% ggplot(aes(x = !!sym(paste0(input$model, input$score_type)))) +
                                geom_density(fill = "#e6ab47",
                                             color = "#e39107",
                                             alpha = 0.8,
                                             adjust = 0.1
                                ) +
                                theme_minimal()
            }
            ggplotly(p)
        }
    })
    
    ### Score deltas distribution plots ----
    output$score_delta_distribution <- renderPlotly({
        data_plot <- data
        if (nrow(df()>0)) {
            plot_df <- df()
            
            # Check if advanced filtering was applied to df(), in that case, also apply it 
            # to the background distribution
            if (length(input$variant_type_filter)>0) {
                try(
                    data_plot <- data_plot %>% filter(!grepl(paste0(input$variant_type_filter, collapse = "|"), data_plot$variant_types))
                )
            }
            try(
                data_plot <- data_plot %>% filter(tsl <= input$tsl_filter | is.na(tsl))
            )
            if (input$tsl_NA) {
                data_plot <- data_plot %>% filter(!grepl("NA",tx_notes))
            }
            if (input$tsl_end_NF) {
                data_plot <- data_plot %>% filter(!grepl("end_NF",tx_notes))
            }
            # Check if filtering options have been applied for the plot
            if (str_detect(input$delta_distribution_filter, "ENST")) {
                plot_df <- plot_df %>% filter(transcript_id == input$delta_distribution_filter)
            }
            if (str_detect(input$delta_distribution_filter, "ENSG")) {
                plot_df <- plot_df %>% filter(haplotype_id == input$delta_distribution_filter)
            }
            p <- data_plot %>%
                ggplot() +
                geom_density(aes(x = !!sym(input$delta)),
                             fill = "#e6ab47",
                             color = "#e39107",
                             alpha = 0.8,
                             adjust = 0.1) +
                theme_minimal() +
                geom_vline(data = plot_df, 
                           aes(xintercept = !!sym(input$delta),
                               text = paste0(haplotype_id, " on transcript ", transcript_id, "\n", input$delta, ": ", round(!!sym(input$delta),4))),
                           color = "#e34907") +
                theme(legend.position = "none")
            ggplotly(p, tooltip = "text")
        }else{
            p <- data_plot %>% ggplot(aes(x = !!sym(input$delta))) +
                geom_density(fill = "#e6ab47", color = "#e39107", alpha = 0.8, adjust = 0.1) +
                theme_minimal()
            ggplotly(p)
        }
    })
    
    ### Update ancestry plot filters ----
    observe({
        # Extract haplotype and trancsript ID from subset df 
        current_df_val <- df()
        transcripts <- c("none", sort(unique(current_df_val$transcript_id)))
        updateSelectizeInput(session,
                             inputId = "ancestry_frequency_filter",
                             label = "Restrict to:",
                             choices = list(
                                 "Transcripts" = as.list(transcripts)
                             ),
                             selected = as.list(transcripts)[2], # the first would be "none"
                             server = FALSE)
        
    })
    
    ### Population frequencies ----
    output$population <- renderPlotly({
        validate(need(try(nrow(df())>0),"Waiting for a subset of haplotypes to be selected."))
        plot_df <- df()
        # Check if filtering options have been applied for the plot
        if (str_detect(input$ancestry_frequency_filter, "ENST")) {
            plot_df <- plot_df %>% filter(transcript_id == input$ancestry_frequency_filter)
        }
        if (input$group_by == "ancestry") {
            # Group by ancestry
            p <- plot_df %>%
                    select(frequency, haplotype_id, AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq) %>%
                    dplyr::rename("Global" = "frequency") %>%
                    pivot_longer(cols = c("Global", "AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "SAS_freq"),
                                 names_to = "ancestry",
                                 values_to = "freq") %>% 
                    ggplot(aes(x = ancestry, y = freq, fill = haplotype_id)) +
                        geom_bar(position = "dodge", stat = "identity") +
                        scale_fill_viridis(discrete = TRUE, option = "F") +
                        theme_minimal() +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            ggplotly(p)
        }else{
            # Group by haplotype
            p <- plot_df %>%
                    select(haplotype_id, AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq) %>%
                    pivot_longer(cols = c("AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "SAS_freq"),
                                 names_to = "ancestry",
                                 values_to = "freq") %>% 
                    ggplot(aes(x = haplotype_id, y = freq, fill = ancestry)) +
                        geom_bar(position = "dodge", stat = "identity") +
                        scale_fill_viridis(discrete = TRUE, option = "F") +
                        theme_minimal() +
                        theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            ggplotly(p)
        }
    })
    
    ### Model correlation selected scores----
    output$models_correlation_selected_score <- renderPlot({
        validate(need(try(nrow(df())>0),"Waiting for a subset of haplotypes to be selected."))
        validate(need(!input$model1_selected=="none" && !input$model2_selected=="none", "Waiting for two models to be selected."))
        p1 <- df() %>% ggplot(aes(x = !!sym(paste0(input$model1_selected, "_PLL")), y = !!sym(paste0(input$model2_selected, "_PLL")))) +
            geom_hex() + 
            scale_fill_gradient(
                low = "gray90",
                high = "#e39107", 
                name = "Count"
            ) +
            theme_minimal()
        p2 <- df() %>% ggplot(aes(x = !!sym(paste0(input$model1_selected, "_PLLR_wt")), y = !!sym(paste0(input$model2_selected, "_PLLR_wt")))) +
            geom_hex() + 
            scale_fill_gradient(
                low = "gray90",
                high = "#e39107",
                name = "Count"
            ) +
            stat_cor(method = "pearson",
                     size = 5) +
            theme_minimal()
        p1 + p2
    })
    
    ## About ----
    
    ### Update value boxes ----
    output$analyzed_genes <- renderText({
        nrow(data %>% select(gene_id) %>% unique())
    })
    output$analyzed_haplotypes <- renderText({
        nrow(data %>% select(haplotype_id) %>% unique())
    })
    output$analyzed_transcripts <- renderText({
        nrow(data %>% select(transcript_id) %>% unique())
    })
    # Read from stats file
    output$analyzed_variants <- renderText({
        nrow(rsid_map)
    })
    
    ### haplotypes per gene distribution ----
    output$haplotype_gene_distribution <- renderPlotly({
        plot_df <- data %>%
                        select(haplotype_id, gene_id) %>%
                        unique() %>%
                        summarise(haplotype_count = n(), .by = gene_id)
        p <- plot_df %>% ggplot(aes(x = haplotype_count)) +
                            geom_histogram(stat = "count",
                                           binwidth = 2,
                                           fill = "#e6ab47",
                                           color = "#e39107",
                                           alpha = 0.8) +
                            theme_minimal()
        ggplotly(p)
    })
    ### variants per haplotype distribution ----
    output$variant_haplotype_distribution <- renderPlotly({
        p <- data %>% ggplot(aes(x = n_variants)) +
            theme_minimal() +
            geom_histogram(stat = "count",
                           binwidth = 2,
                           fill = "#e6ab47",
                           color = "#e39107",
                           alpha = 0.8)
        ggplotly(p)
    })
    
    ### Model correlation all scores----
    output$models_correlation_all_score <- renderPlot({
        validate(need(!input$model1_all=="none" && !input$model2_all=="none", "Waiting for two models to be selected."))
        p1 <- data %>% ggplot(aes(x = !!sym(paste0(input$model1_all, "_PLL")), y = !!sym(paste0(input$model2_all, "_PLL")))) +
                                geom_hex() + 
                                scale_fill_gradient(
                                    low = "gray90",
                                    high = "#e39107", 
                                    name = "Count"
                                ) +
                                theme_minimal()
        p2 <- data %>% ggplot(aes(x = !!sym(paste0(input$model1_all, "_PLLR_wt")), y = !!sym(paste0(input$model2_all, "_PLLR_wt")))) +
                                geom_hex() + 
                                scale_fill_gradient(
                                    low = "gray90",
                                    high = "#e39107",
                                    name = "Count"
                                ) +
                                stat_cor(method = "pearson",
                                         size = 5) +
                                theme_minimal()
        p1 + p2
    })
    
    ## Access ----
    
    output$version <- renderTable(version)
    output$api <- renderTable(api)
    ## Download ----
    
    ### Subset download datatable ----
    preview_df <- reactive({
        if (length(input$IDs_selector)>0) {
            preview_df <- data %>% select(input$IDs_selector)
        }
        if (length(input$variant_info_selector)>0) {
            preview_df <- cbind(preview_df, data %>% select(input$variant_info_selector))
        }
        if (length(input$frequency_selector)>0) {
            preview_df <- cbind(preview_df, data %>% select(input$frequency_selector))
        }
        if (length(input$transcript_info_selector)>0) {
            preview_df <- cbind(preview_df, data %>% select(input$transcript_info_selector))
        }
        if (length(input$sequence_selector)>0) {
            sequences_df <- read_tsv(sequences_path) %>% select(gene_id, transcript_id, haplotype_id, input$sequence_selector)
            preview_df <- preview_df %>% left_join(sequences_df)
        }
        if (length(input$scores_selector)>0) {
            preview_df <- cbind(preview_df, data %>% select(input$scores_selector))
        }
        
        return(preview_df)
    })
    
    ### Download datatable ----
    output$preview_table <- renderDataTable(
        preview_df() %>% head()
    )
    
    ### Download button ----
    output$download_button <- downloadHandler(
        filename = "scores.tsv",
        content = function(file){
            write_tsv(preview_df(), file)
        }
    )
    
}

# RUN APP ----
shinyApp(ui, server)
