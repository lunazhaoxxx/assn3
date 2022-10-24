library(tidyverse)

#Read in the augustus prediction file
s_sup <- read_tsv(file = "augustus_prediction.csv") %>% separate(1, c("gene", "transcript"), sep = "\\.")

#Rename columns and select only columns we need
colnames(s_sup)[3:9] <- c("support",
                         "3'UTR exons and introns", 
                         "5'UTR exons and introns",
                         "CDS exons",
                         "CDS introns",
                         "obeyed",
                         "incompatible")

# predicted gene number
nrow(unique(s_sup[1]))

# predicted transcript number
nrow(s_sup[2])

# calculate the proportion of transcripts with support > 80
sub <- s_sup %>%
  filter(support > 80)
nrow(sub)/nrow(s_sup)

# calculate the number and proportion of transcripts with obeyed < 4
sub <- s_sup %>%
  filter(obeyed < 4)
nrow(sub)
nrow(sub)/nrow(s_sup)

# distribution of percentage of transcript supported by hints
ggplot(s_sup, aes(x=support)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=15)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("% of transcript supported by hints") + 
  ylab("Density")

# distribution of the number of hint groups fully obeyed
ggplot(s_sup, aes(x=obeyed)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=10)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("The number of hint groups fully obeyed") + 
  ylab("Density")

# distribution of the number of incompatible hint groups 
ggplot(s_sup, aes(x=incompatible)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=35)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("The number of incompatible hint groups") + 
  ylab("Density")

#Plot number of incompatible hint groups per transcript against % of support by hints
s_sup %>%
  ggplot(aes(x = support, y = incompatible, size = 2)) +
  geom_point() + 
  theme_minimal() +
  xlab("% of transcript supported by hints") +
  ylab("Number of incompatible hint groups") +
  guides(size = "none")

#Plot number of fully obeyed hint groups per transcript against % of support by hints
s_sup %>%
  ggplot(aes(x = support, y = obeyed, size = 2)) +
  geom_point() +
  theme_minimal() +
  xlab("% of transcript supported by hints") +
  ylab("Number of obeyed hint groups") +
  guides(size = "none")

#Plot number of incompatible hint groups vs the number of obeyed hint groups
s_sup %>%
  ggplot(aes(x = obeyed, y = incompatible, size = 2)) +
  geom_point() +
  theme_minimal() +
  xlab("Number of obeyed hint groups") +
  ylab("Number of incompatible hint groups") +
  guides(size = "none")

#Read in blast results
blast <- read_tsv("blast_species_descr.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "description", "species", "qlen"))

#Plot based on which related species the protein sequences have been annotated
blast %>%
  select(species) %>%
  group_by(species) %>%
  count() %>%
  arrange(desc(n)) %>% 
  head(10) %>% 
  ggplot(aes(x = species, y = n)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + xlab("Species") + 
  ylab("Count") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# keep blastq results with highest bitscore for each predicted transcript
# get top 10 species in blast alignment
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  select(species) %>%
  group_by(species) %>%
  count() %>%
  arrange(desc(n)) %>% 
  head(10) %>% 
  ggplot(aes(x = species, y = n)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + xlab("Species") + 
  ylab("Count") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# keep blastq results with highest bitscore for each predicted transcript
# distribution of  pident (percentage of identical matches)
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x=pident)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=10)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("Percentage of identical matches") + 
  ylab("Density")

# keep blastq results with highest bitscore for each predicted transcript
# distribution of alignment length
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x=length)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=50)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("Alignment length") + 
  ylab("Density")

# keep blastq results with highest bitscore for each predicted transcript
# distribution of the number of gap openings
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x=gapopen)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=1)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("The number of gap openings") + 
  ylab("Density")

# keep blastq results with highest bitscore for each predicted transcript
# Plot alignment length against % of identical matches
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x = pident, y = length)) +
  geom_point(size = 2) + 
  theme_minimal() +
  xlab("% of identical matches") +
  ylab("Alignment length") +
  guides(size = "none")+
  geom_density_2d(alpha=0.5)

# keep blastq results with highest bitscore for each predicted transcript
# Plot the number of gap openings against % of identical matches
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x = pident, y = gapopen)) +
  geom_point(size = 2) + 
  theme_minimal() +
  xlab("% of identical matches") +
  ylab("The number of gap openings") +
  guides(size = "none") + 
  geom_density_2d(alpha=0.5)

# keep blastq results with highest bitscore for each predicted transcript
# Plot bitscore against % of identical matches
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x = pident, y = bitscore)) +
  geom_point(size = 2) + 
  theme_minimal() +
  xlab("% of identical matches") +
  ylab("bitscore") +
  guides(size = "none") + 
  geom_density_2d(alpha=0.5)

# distribution of  query sequence length (predicted transcript length)
blast %>%
  group_by(qseqid) %>%
  top_n(1, bitscore) %>%
  ggplot(aes(x=qlen)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=70)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("The length of predicted transcripts") + 
  ylab("Density")

# Read in PFAM HMM search results
pfam <- read_tsv("hmmer_out.tsv", col_names = c("targetName", "accession", "tlen", "queryName", "qlen", "eval", "ndom", "totalDom", "hmmStart", "hmmEnd", "aliStart", "aliEnd", "accuracyPP", "description"))

# plot the distribution of accuracyPP (reliability of alignment) of all annotation results
pfam %>%
  ggplot(aes(x=accuracyPP)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=0.1)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("Reliability of alignment (posterior probability)") + 
  ylab("Density")

#Count occurrences of individual PFAM accession numbers
accession_freq <- pfam %>%
  group_by(accession) %>%
  count() %>%
  arrange(desc(n)) %>% 
  head(10) 

colnames(accession_freq) <- c('accession', "count")

accession_freq <- unique(left_join(accession_freq, pfam[,c(2,14)], by='accession', all.x = T, no.dups = T))

#plot description of top 10 of the mostly commonidentified protein
accession_freq %>%
  ggplot(aes(x = reorder(description.x, -n), y = n)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + xlab("Description") + 
  ylab("Count") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

#Calculate the percentage of the protein sequence that is covered by the HMM domain model hit
cov <- pfam %>% select(queryName, qlen, aliStart, aliEnd) %>%
  rowwise() %>%
  mutate(cov = list(c(rep(0,aliStart-1), rep(1, aliEnd-(aliStart-1)), rep(0, qlen-aliEnd)))) %>%
  group_by(queryName,qlen) %>%
  summarise(collapse = list(reduce(.x = cov, .f = `+`))) %>%
  unnest() %>% 
  filter(collapse > 0) %>%
  mutate(collapse = 1) %>% 
  group_by(queryName,qlen) %>%
  summarise(cov = sum(collapse)) %>% 
  mutate(cov = (cov/qlen)*100)

sub <- cov %>%
  filter(cov > 60)
nrow(sub)
nrow(cov)
nrow(sub)/89

# plot the distribution of coverage from pfam for each transcript
ggplot(cov, aes(x=cov)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=10)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("% of sequence covered by Pfam annotation") + 
  ylab("Density")

# Plot the largest 10 coverage from pfam and their transcript ids.
cov %>%   
  arrange(desc(cov)) %>% 
  head(20) %>%
  ggplot(aes(x = reorder(queryName, desc(cov)), y = cov)) +
  geom_histogram(stat = "identity") +
  theme_classic() +
  xlab("Protein identifier") +
  ylab("% of sequence covered by Pfam annotation") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

#Calculate the percentage of the protein sequence that is covered by the BLAST hit
cov_blast <- blast %>% select(qseqid, qlen, qstart, qend) %>%
  rowwise() %>%
  mutate(cov = list(c(rep(0,qstart-1), rep(1, qend-(qstart-1)), rep(0, qlen-qend)))) %>%
  group_by(qseqid, qlen) %>%
  summarise(collapse = list(reduce(.x = cov, .f = `+`))) %>%
  unnest() %>% 
  filter(collapse > 0) %>%
  mutate(collapse = 1) %>% 
  group_by(qseqid, qlen) %>%
  summarise(cov_blast = sum(collapse)) %>% 
  mutate(cov_blast = (cov_blast/qlen)*100)

sub <- cov_blast %>% 
  filter(cov_blast > 99)
nrow(sub)
nrow(cov_blast)
nrow(sub)/89

# plot the distribution of coverage from BLAST for each transcript
ggplot(cov_blast, aes(x=cov_blast)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=1) +
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("% of sequence covered by BLAST") + 
  ylab("Density")

# Plot the largest 10 coverage from BLAST and their transcript ids.
cov_blast %>% 
  arrange(desc(cov_blast)) %>% 
  head(20) %>%
  ggplot(aes(x = reorder(qseqid, desc(cov_blast)), y = cov_blast)) +
  geom_histogram(stat = "identity") +
  theme_minimal() +
  xlab("Protein identifier") +
  ylab("% of sequence covered by BLAST") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# Combine the results from BLAST and PFAM and see if coverage from the two methods correlates
full_join(cov, cov_blast, by=c("queryName"="qseqid")) %>%
  mutate(cov = replace_na(cov, 0)) %>% 
  mutate(cov_blast = replace_na(cov_blast, 0)) %>% 
  select(-qlen.x,-qlen.y) %>%
  rename(cov_ipro = cov) %>%
  mutate(cov_blast = if_else(is.na(cov_blast), 0, cov_blast)) %>%
  ggplot(aes(x = cov_ipro, y = cov_blast)) +
  geom_point(size = 5) +
  coord_fixed(1) +
  xlim(0,100) +
  ylim(0,100) +
  theme_minimal() +
  xlab("% sequence coverage by Pfam domains") +
  ylab("% sequence coverage by blastp")

tmp_table <- full_join(cov, cov_blast, by=c("queryName"="qseqid")) %>%
  mutate(cov = replace_na(cov, 0)) %>% 
  mutate(cov_blast = replace_na(cov_blast, 0)) %>% 
  select(-qlen.x,-qlen.y) %>%
  rename(cov_ipro = cov) %>%
  mutate(cov_blast = if_else(is.na(cov_blast), 0, cov_blast))

# the number of transcripts with results from blast and pfam
both <- tmp_table %>%
  filter(cov_blast > 0) %>%
  filter(cov_ipro > 0)
nrow(both)
69/89

# the number of transcripts with results from pfam
pfam_table <- tmp_table %>%
  filter(cov_ipro > 0)
nrow(pfam_table)

# the number of transcripts with results from pfam
blast_table <- tmp_table %>%
  filter(cov_blast > 0)

nrow(blast_table) - nrow(both)
(nrow(blast_table) - nrow(both))/89

# join blast and pfam results and get all information
merge_all <- full_join(blast, pfam, by=c("qseqid"="queryName"))
# plot the correlation between blast bitscore and pfam accuracyPP
merge_all %>%
  ggplot(aes(x = bitscore, y = accuracyPP, size = 2)) +
  geom_point(alpha = 0.3, size = 0.5) +
  theme_minimal() +
  xlab("Bitscore from BLAST") +
  ylab("Reliability (accuracyPP) from pfam") +
  guides(size = "none") +
  geom_density_2d(alpha = 0.5)

# plot distribution of length 
qlen <- read_tsv("qlen.tsv", col_names = c("queryName", "qlen"))
nrow(qlen)
mean(qlen$qlen)
median(qlen$qlen)
ggplot(qlen, aes(x=qlen)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=70)+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_classic() + 
  xlab("Transcript_length") + 
  ylab("Density")


