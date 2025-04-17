library(dplyr)
library(tidyr)
library(ggplot2)
library(ggVennDiagram)
source("library/library.R")


metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA-NH_02102025.csv")
head(metadata)
metadata$BaselineType = metadata$dateOrBaselineType
metadata$BaselineType[grep("^[0-9].*", metadata$BaselineType)] = "sample"

metadata$date = paste0(metadata$dateOrBaselineType, "-2024")
metadata$date[grep("^BL.", metadata$date)] = NA
#metadata$mdy = mdy(metadata$date)



asv_tab.bl = read.csv("data/FEDRR_all_2024/asv_tab.baseline.csv", row.names = 1)
asv_tab.bl$SequenceID = rownames(asv_tab.bl)
asv_tab.samps = read.csv("data/FEDRR_all_2024/asv_tab.samps.csv", row.names = 1)
asv_tab.samps$SequenceID = rownames(asv_tab.samps)

#we can just join to trapID and then make a vec of the ASVs that are present
asv_tab.bl.meta = left_join(asv_tab.bl, metadata %>% select(SequenceID, State, BaselineType, trapID))
asv_tab.samps.meta = left_join(asv_tab.samps, metadata %>% select(SequenceID, State, BaselineType, trapID))

samps.asv_sum = asv_tab.samps.meta %>% 
    group_by(State, BaselineType, trapID) %>%
    summarize(across(where(is.numeric), .fns = sum)) %>%
    data.frame

bl.asv_sum = asv_tab.bl.meta %>% 
    group_by(State, BaselineType, trapID) %>%
    summarize(across(where(is.numeric), .fns = sum)) %>%
    data.frame

#list 
sample_list = list()
all_traps = samps.asv_sum$trapID
samps_asvs = colnames(samps.asv_sum)[4:ncol(samps.asv_sum)]
length(samps_asvs)
bl_asvs = colnames(bl.asv_sum)[4:ncol(bl.asv_sum)]
#sanity
length(samps_asvs) == ncol(samps.asv_sum)-3
length(bl_asvs) == ncol(bl.asv_sum)-3

for(i in 1:length(all_traps) ){
    pres_abs.samp = samps.asv_sum[samps.asv_sum$trapID == all_traps[i], 4:ncol(samps.asv_sum)]
    pres_abs.bl = bl.asv_sum[bl.asv_sum$trapID == all_traps[i], 4:ncol(bl.asv_sum)]
    temp_list = list()
    temp_list[["samples"]] = samps_asvs[pres_abs.samp > 0]
    temp_list[["baseline"]] = bl_asvs[pres_abs.bl > 0]
    sample_list[[all_traps[i]]] = temp_list
}

# how many ASVs per bl sample per trap
for(i in 1:length(all_traps)) {
    print(paste(all_traps[i], length( sample_list[[ all_traps[i] ]][["baseline"]] )) )
}

getOption('repr.plot.width') #default 7


p1 = ggVennDiagram(sample_list[[1]], label_alpha = 0.85, label_size = 6) +
    scale_fill_distiller(palette = "RdBu") +
    labs(title = paste("Trap", all_traps[1]), fill = "ASV count") +
    coord_cartesian(clip = "off") +
    theme(
        plot.title = element_text(size = 18),
        legend.position = "right",
        aspect.ratio = 1.35
    )
p1


#pdf("figures/FEDRR_all_2024/test_baseline_Venn.pdf", width = 5, height = 5)
#print(p1)
#dev.off()


# save plots list
plot_list = list()
for(i in 1:length(all_traps)) {
    
    plot_list[[ all_traps[i] ]] = ggVennDiagram(
            sample_list[[i]], 
            label_alpha = 0.85, 
            label_size = 5.5
        ) +
        scale_fill_distiller(palette = "RdBu") +
        labs(title = paste("Trap", all_traps[i]), fill = "ASV count") +
        coord_cartesian(clip = "off") +
        theme(
            plot.title = element_text(size = 18),
            legend.position = "right",
            aspect.ratio = 1.35
        )
}

pdf("figures/FEDRR_all_2024/baseline_Venns.NY-OH-IA-NH.pdf", width = 5, height = 5)
lapply(plot_list, print)
dev.off()


###############################################################
# Venn of all samples/traps
###############################################################
samps_asvs.sums = asv_tab.samps.meta %>% 
    select(where(is.numeric)) %>% 
    colSums()
bl_asvs.sums = asv_tab.bl.meta %>% 
    select(where(is.numeric)) %>% 
    colSums()
length(samps_asvs.sums[samps_asvs.sums > 0]) == length(samps_asvs.sums)
length(bl_asvs.sums[bl_asvs.sums > 0]) == length(bl_asvs.sums)

samps_asvs = names(samps_asvs.sums)
bl_asvs = names(bl_asvs.sums)

length(samps_asvs)
length(bl_asvs)
length(samps_asvs) - length(bl_asvs)

p1 = ggVennDiagram(
            list(samples = samps_asvs, baseline = bl_asvs), 
            label_alpha = 0.85, 
            label_size = 5.5
        ) +
        scale_fill_distiller(palette = "RdBu") +
        labs(fill = "ASV count") +
        coord_cartesian(clip = "off") +
        theme(
            plot.title = element_text(size = 18),
            legend.position = "right",
            aspect.ratio = 1.35
        )

pdf("figures/FEDRR_all_2024/baseline_all_samples_Venn.NY-OH-IA-NH.pdf", width = 5, height = 5)
p1
dev.off()


###############################################################
###############################################################
# also output lists of ASVs that are shared to look at taxonomy

shared_list = list()

for(i in 1:length(all_traps)) {
    temp_shared = sample_list[[ all_traps[i] ]][["baseline"]] %in% sample_list[[ all_traps[i] ]][["samples"]]
    shared_list[[ all_traps[i] ]] = sample_list[[ all_traps[i] ]][["baseline"]][temp_shared]
}
unlist(shared_list, use.names = F) %>% length
#3880
unlist(shared_list, use.names = F) %>% unique %>% length
#1664
shared_asvs = unlist(shared_list, use.names = F) %>% unique

head(asv_tab.bl)
######################
#reloading to refresh
asv_tab.bl = read.csv("data/FEDRR_all_2024/asv_tab.baseline.csv", row.names = 1)
asv_tab.bl.t = t(asv_tab.bl)
asv_tab.bl.t.shared = asv_tab.bl.t[rownames(asv_tab.bl.t) %in% shared_asvs,]
nrow(asv_tab.bl.t.shared)
nrow(asv_tab.bl.t) - nrow(asv_tab.bl.t.shared)
#dropped 2898 ASVs
nrow(asv_tab.bl.t[rowSums(asv_tab.bl.t) > 0,])
nrow(asv_tab.bl.t)
#all ASVs > 0 seqs

######################
#loading taxonomic data
asv_tax = read.table("data/FEDRR_all_2024/ASVs_taxonomy.tsv", header = T)
tax_bs = read.table("data/FEDRR_all_2024/ASVs_taxonomy_bootstrapVals.tsv", header = T)
asv_tax[tax_bs < 80] = NA
asv_tax.shared = asv_tax[rownames(asv_tax) %in% shared_asvs,]
nrow(asv_tax.shared)

asv_tax.shared$ASV = rownames(asv_tax.shared)
asv_tab.bl.t.shared = as.data.frame(asv_tab.bl.t.shared)
asv_tab.bl.t.shared$ASV = rownames(asv_tab.bl.t.shared)

asv_tab_tax = left_join(asv_tax.shared, asv_tab.bl.t.shared, by = join_by(ASV))
asv_tab_tax.sp = asv_tab_tax %>% 
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarize(across(where(is.numeric), sum) )
nrow(asv_tab_tax.sp)


################################
################################
# tax level summaries
# 
# First calculate rel abd
seqs_per_sample = colSums(asv_tab_tax %>% select(where(function(x) is.numeric(x) == T)))
seqs_per_sample > 0

asv_tab_tax.RA = cbind(
    asv_tab_tax %>% select(where(function(x) is.numeric(x) == F)),
    apply(
        asv_tab_tax %>% select(where(function(x) is.numeric(x) == T)), 
        1, 
        FUN = function(x) x/seqs_per_sample
    ) %>% t
)
colSums(asv_tab_tax.RA %>% select(where(function(x) is.numeric(x) == T)))

asv_tab_tax.order = asv_tab_tax.RA %>% 
    #filter(!is.na(Class)) %>%
    group_by(Class) %>%
    summarize(across(where(is.numeric), sum) ) %>%
    pivot_longer(cols = -Class, names_to = "SequenceID", values_to = "rel_abd") %>%
    left_join(., metadata %>% select(SequenceID, State, trapID, BaselineType))

asv_tab_tax.order %>%
    group_by(trapID, BaselineType) %>%
    summarize(sum(rel_abd)) %>%
    print(n = Inf)
# parsing most abundant taxa through dif ranks

#classes
p1 = ggplot(
    asv_tab_tax.RA %>% 
        #filter(!is.na(Class)) %>%
        group_by(Class) %>%
        summarize(across(where(is.numeric), sum) ) %>%
        pivot_longer(cols = -Class, names_to = "SequenceID", values_to = "rel_abd") %>%
        left_join(., metadata %>% select(SequenceID, State, trapID, BaselineType)), 
    aes(x = paste(trapID, BaselineType), y = rel_abd*100, fill = Class)
) +
    geom_col() +
    #scale_fill_manual(values = c(c25, "white", "grey", "black")) +
    labs(x = "Trap", y = "Relative abundance (%)", title = "Classes") +
    guides(fill=guide_legend(ncol=1)) +
    my_gg_theme +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
p1

#order of dothideos
p2 = ggplot(
    asv_tab_tax.RA %>% 
        filter(Class == "c__Dothideomycetes") %>%
        group_by(Order) %>%
        summarize(across(where(is.numeric), sum) ) %>%
        pivot_longer(cols = -Order, names_to = "SequenceID", values_to = "rel_abd") %>%
        left_join(., metadata %>% select(SequenceID, State, trapID)), 
    aes(x = trapID, y = rel_abd*100, fill = Order)
) +
    geom_col() +
    scale_fill_manual(values = c(c25, "white")) +
    labs(x = "Trap", y = "Relative abundance (%)", title = "Dothideomycetes orders") +
    guides(fill=guide_legend(ncol=1)) +
    my_gg_theme +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
p2
#fam of pleosporales
p3 = ggplot(
    asv_tab_tax.RA %>% 
        filter(Order == "o__Pleosporales") %>%
        group_by(Family) %>%
        summarize(across(where(is.numeric), sum) ) %>%
        pivot_longer(cols = -Family, names_to = "SequenceID", values_to = "rel_abd") %>%
        left_join(., metadata %>% select(SequenceID, State, trapID)), 
    aes(x = trapID, y = rel_abd*100, fill = Family)
) +
    geom_col() +
    scale_fill_manual(values = c(c25, "white")) +
    labs(x = "Trap", y = "Relative abundance (%)", title = "Pleosporales families") +
    guides(fill=guide_legend(ncol=1)) +
    my_gg_theme +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
p3

asv_tab_tax.RA %>% filter(Family == "f__Didymellaceae")
# g__Neoascochyta s__rosicola, g__Neoascochyta NA, NA
asv_tab_tax.RA %>% filter(Family == "f__Didymosphaeriaceae")
# g__Paracamarosporium and NA

#fam of capnodiales
p4 = ggplot(
    asv_tab_tax.RA %>% 
        filter(Order == "o__Capnodiales") %>%
        group_by(Family) %>%
        summarize(across(where(is.numeric), sum) ) %>%
        pivot_longer(cols = -Family, names_to = "SequenceID", values_to = "rel_abd") %>%
        left_join(., metadata %>% select(SequenceID, State, trapID)), 
    aes(x = trapID, y = rel_abd*100, fill = Family)
) +
    geom_col() +
    scale_fill_manual(values = c(c25, "white")) +
    labs(x = "Trap", y = "Relative abundance (%)", title = "Capnodiales families") +
    guides(fill=guide_legend(ncol=1)) +
    my_gg_theme +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
p4

#order of sordarios
p5 = ggplot(
    asv_tab_tax.RA %>% 
        filter(Class == "c__Sordariomycetes") %>%
        group_by(Order) %>%
        summarize(across(where(is.numeric), sum) ) %>%
        pivot_longer(cols = -Order, names_to = "SequenceID", values_to = "rel_abd") %>%
        left_join(., metadata %>% select(SequenceID, State, trapID)), 
    aes(x = trapID, y = rel_abd*100, fill = Order)
) +
    geom_col() +
    scale_fill_manual(values = c(c25, "white")) +
    labs(x = "Trap", y = "Relative abundance (%)", title = "Sordariomycetes orders") +
    guides(fill=guide_legend(ncol=1)) +
    my_gg_theme +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
p5
#fam of sordariales
p6 = ggplot(
    asv_tab_tax.RA %>% 
        filter(Order == "o__Sordariales") %>%
        group_by(Family) %>%
        summarize(across(where(is.numeric), sum) ) %>%
        pivot_longer(cols = -Family, names_to = "SequenceID", values_to = "rel_abd") %>%
        left_join(., metadata %>% select(SequenceID, State, trapID)), 
    aes(x = trapID, y = rel_abd*100, fill = Family)
) +
    geom_col() +
    scale_fill_manual(values = c(c25, "white")) +
    labs(x = "Trap", y = "Relative abundance (%)", title = "Sordariales families") +
    guides(fill=guide_legend(ncol=1)) +
    my_gg_theme +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
p6
asv_tab_tax.RA %>% filter(Family == "f__Chaetomiaceae")
# g__Collariella s__carteri, NA
asv_tab_tax.RA %>% filter(Family == "f__Sordariaceae")
# NA


#sort sp df to look at top taxa
# could just sum to get one column of totals in stead of using across for this ...


asv_tab_tax %>% 
    mutate(gen_sp = paste(
        sub("g__", "", Genus),
        sub("s__", "", Species)
        )
    ) %>% 
    group_by(gen_sp) %>%
    dplyr::summarize(total_seqs = sum(across(where(is.numeric))) ) %>%
    mutate(rel_abd = total_seqs/sum(total_seqs)) %>% # had some funny results here where . was needed and then not.
    slice_max(order_by = rel_abd, n = 25)


p7 = ggplot(
    asv_tab_tax %>% 
        mutate(gen_sp = paste(
            sub("g__", "", Genus),
            sub("s__", "", Species)
            )
        ) %>% 
        group_by(gen_sp) %>%
        dplyr::summarize(total_seqs = sum(across(where(is.numeric))) ) %>%
        mutate(rel_abd = total_seqs/sum(total_seqs)) %>% # had some funny results here where . was needed and then not.
        slice_max(order_by = rel_abd, n = 25), 
    aes(x = reorder(gen_sp, rel_abd), y = rel_abd*100)) +
    geom_col() +
    my_gg_theme +
    labs(y = "Relative abundance (%)", title = "Top 25 genus/species groups") +
    theme(
        axis.text.x = element_text(angle = 65, hjust = 1),
        axis.title.x = element_blank()
    )



pdf("figures/FEDRR_all_2024/baseline_taxa.NY-OH.pdf", width = 14, height = 8)
p1
p2
p3
p4
p5
p6
p7
dev.off()




